# Central reader for environmental metadata across datasets.
#
# Usage (examples):
#   Rscript code/read_data/read_envdata.R --config code/read_data/config/ocean.R
#   Rscript code/read_data/read_envdata.R --config code/read_data/config/wwtp.R
#   Rscript code/read_data/read_envdata.R --config code/read_data/config/soil.R
#
# Outputs:
# - Writes `envdata.csv` and `site_locations.csv` (sample_id, latitude, longitude) into the dataset-specific `data/processed_data/` directory
# - Returns a list (when sourced/used as functions) containing:
#   - envdata (tibble): `sample_id` + selected environmental variables
#   - var_catalog (tibble): variable documentation (type/units/description/source column)
#   - extras (list): dataset-specific extras (e.g., Tara `sample_id_lookup`)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(readxl)
})

or_default <- function(x, y) if (!is.null(x)) x else y

stopf <- function(...) stop(sprintf(...), call. = FALSE)

safe_as_numeric <- function(x) {
  # Convert common non-numeric encodings to NA without warnings.
  # This avoids "NAs introduced by coercion" from as.numeric().
  s <- trimws(as.character(x))
  s[s %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  s <- gsub(",", "", s, fixed = TRUE)

  # Accept plain numbers and scientific notation.
  ok <- !is.na(s) & grepl("^[+-]?(\\d+\\.?\\d*|\\d*\\.\\d+)([eE][+-]?\\d+)?$", s)
  out <- rep(NA_real_, length(s))
  out[ok] <- as.numeric(s[ok])
  out
}

as_logical_flag <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(x != 0)
  if (is.character(x)) {
    v <- tolower(trimws(x))
    if (v %in% c("true", "t", "1", "yes", "y")) return(TRUE)
    if (v %in% c("false", "f", "0", "no", "n")) return(FALSE)
  }
  stopf("Could not parse logical flag value: '%s'", x)
}

parse_cli_args <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (!requireNamespace("optparse", quietly = TRUE)) {
    stopf("R package 'optparse' is required. Install it with: install.packages('optparse')")
  }

  option_list <- list(
    optparse::make_option(
      c("-c", "--config"),
      type = "character",
      help = "Path to an R config file that evaluates to a list()"
    ),
    optparse::make_option(
      "--root",
      type = "character",
      default = NA_character_,
      help = "Repo root path (optional; if omitted, uses repo root from code/setup.R)"
    ),
    optparse::make_option(
      "--dry_run",
      action = "store_true",
      default = FALSE,
      help = "If set, do not write any outputs"
    )
  )

  p <- optparse::OptionParser(option_list = option_list)
  a <- optparse::parse_args(p, args = argv, positional_arguments = FALSE)

  # Normalize to the keys used elsewhere in the script.
  out <- list(
    config = a$config,
    root = a$root,
    dry_run = a$dry_run
  )

  out
}

load_config <- function(config_path) {
  if (is.null(config_path) || !nzchar(config_path)) stopf("Missing required --config <path>")
  if (!file.exists(config_path)) stopf("Config file not found: %s", config_path)

  cfg <- source(config_path, local = TRUE)$value
  if (!is.list(cfg)) stopf("Config must evaluate to an R list: %s", config_path)
  cfg
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

print_var_catalog <- function(var_catalog) {
  vc <- var_catalog %>%
    mutate(type = as.character(.data$type)) %>%
    arrange(.data$type, .data$var)
  cat("\n=== Variable catalog (by type) ===\n")
  for (ty in unique(vc$type)) {
    vars <- vc %>% filter(.data$type == ty) %>% pull(.data$var)
    cat(sprintf("- %s (%d): %s\n", ty, length(vars), paste(vars, collapse = ", ")))
  }
}

apply_rename_map <- function(df, rename_map = c()) {
  # Helper to apply a config-specified rename map consistently across datasets.
  # Convention: rename_map is a named character vector: old_name = "new_name"
  # (i.e., names(rename_map) are the source columns; values are the canonical names)
  rename_map <- or_default(rename_map, c())
  if (length(rename_map) == 0) return(df)

  old_names <- intersect(names(rename_map), names(df))
  if (length(old_names) == 0) return(df)

  # dplyr::rename expects: new_name = old_name
  df %>% rename(!!!setNames(old_names, rename_map[old_names]))
}

get_var_catalog <- function(cfg, dataset, fallback = NULL) {
  vc <- or_default(cfg$var_catalog, fallback)
  if (is.null(vc)) {
    stopf("Missing var_catalog: add cfg$var_catalog to the config for dataset '%s'", dataset)
  }

  vc <- as_tibble(vc)

  # If dataset column is missing, fill it
  if (!"dataset" %in% names(vc)) vc$dataset <- dataset

  required <- c("dataset", "var", "type", "units", "description", "source_column")
  missing <- setdiff(required, names(vc))
  if (length(missing) > 0) {
    stopf("cfg$var_catalog is missing required columns: %s", paste(missing, collapse = ", "))
  }

  vc
}

apply_range_filters <- function(df, ranges = list()) {
  # Apply per-variable numeric range filters.
  #
  # `ranges` is a named list like:
  #   ranges = list(
  #     ph = list(min = 3, max = 10, keep_na = FALSE),
  #     soil_c = list(max = 15, keep_na = TRUE)
  #   )
  #
  # Key detail: if keep_na is FALSE, rows with NA in that variable are *dropped*.
  # (This used to be buggy because NA values in the logical mask can propagate and
  # data.frame subsetting with NA keeps "NA rows".)

  if (length(ranges) == 0) return(list(df = as_tibble(df), removed = tibble()))
  removed_log <- list()

  for (v in names(ranges)) {
    if (!v %in% names(df)) {
      removed_log[[v]] <- tibble(var = v, rule = "range", removed_n = NA_integer_)
      next
    }

    r <- ranges[[v]]
    minv <- or_default(r$min, -Inf)
    maxv <- or_default(r$max, Inf)
    keep_na <- as_logical_flag(or_default(r$keep_na, TRUE))

    before_n <- nrow(df)

    # Ensure numeric comparisons behave consistently, even if the column is read as text.
    x <- df[[v]]
    if (!is.numeric(x)) x <- safe_as_numeric(x)

    # IMPORTANT: make the keep mask strictly TRUE/FALSE (no NA).
    keep <- (!is.na(x) & x >= minv & x <= maxv) | (keep_na & is.na(x))

    df <- df[keep, , drop = FALSE]
    after_n <- nrow(df)

    removed_log[[v]] <- tibble(
      var = v,
      rule = sprintf(
        "[%s,%s]%s",
        as.character(minv),
        as.character(maxv),
        if (keep_na) " (NA kept)" else " (NA dropped)"
      ),
      removed_n = before_n - after_n
    )
  }

  list(df = as_tibble(df), removed = bind_rows(removed_log))
}

choose_variables <- function(envdata, base_vars, extra_vars = character()) {
  # Manual-only selection:
  # - Keep base vars (ordered) + extra vars (ordered) if present

  base_vars <- unique(base_vars)
  extra_vars <- unique(extra_vars)

  missing_base <- setdiff(base_vars, names(envdata))
  if (length(missing_base) > 0) {
    cat("\nWARNING: Missing base variables in envdata:\n")
    cat(paste0("  - ", missing_base, collapse = "\n"), "\n")
  }

  missing_extra <- setdiff(extra_vars, names(envdata))
  if (length(missing_extra) > 0) {
    cat("\nNOTE: Missing extra variables in envdata:\n")
    cat(paste0("  - ", missing_extra, collapse = "\n"), "\n")
  }

  kept <- c(intersect(base_vars, names(envdata)), intersect(extra_vars, names(envdata)))
  list(kept_vars = kept)
}

write_outputs <- function(output_dir, envdata, var_catalog, extras = list(), site_locations = NULL, dry_run = FALSE) {
  if (dry_run) {
    cat("\n(dry run) Not writing outputs.\n")
    return(invisible(NULL))
  }
  ensure_dir(output_dir)
  write.csv(envdata, file.path(output_dir, "envdata.csv"), row.names = FALSE)
  write.csv(var_catalog, file.path(output_dir, "var_catalog.csv"), row.names = FALSE)

  if (!is.null(site_locations) && nrow(site_locations) > 0) {
    write.csv(site_locations, file.path(output_dir, "site_locations.csv"), row.names = FALSE)
  }

  if (!is.null(extras$sample_id_lookup)) {
    write.csv(extras$sample_id_lookup, file.path(output_dir, "sample_id_lookup.csv"), row.names = FALSE)
  }

  cat(sprintf("\nWrote outputs to: %s\n", output_dir))
}

read_env_soil <- function(root, cfg) {
  in_path <- file.path(root, cfg$input$path)
  if (!file.exists(in_path)) {
    stopf(
      "soil input not found: %s\nSet cfg$input$path to the correct location (default expects soil/data/downloaded_data/env_metadata.xlsx).",
      in_path
    )
  }

  raw <- readxl::read_xlsx(in_path)
  # Standardize column names using the config rename_map
  raw <- apply_rename_map(raw, cfg$rename_map)

  # Canonical required id
  if (!"sample_id" %in% names(raw)) stopf("soil: expected a sample id column mapped to 'sample_id'")

  # Coerce types
  raw <- raw %>% mutate(sample_id = as.character(.data$sample_id))

  # Extract the variable catalog from config
  var_catalog <- get_var_catalog(cfg, dataset = "soil")

  list(
    envdata_raw = as_tibble(raw),
    var_catalog = var_catalog,
    extras = list()
  )
}

read_tara_carbonate_xlsx <- function(root, cfg, sample_id_lookup) {
  # Read TARA carbonate chemistry XLSX (PARAMETER-row layout, COMMENT row for min/q1/q2/q3/max).
  # Returns tibble with sample_id + carbonate columns (q2/median only), canonical names via rename_map.
  carb_path <- file.path(root, cfg$input$carbonate_xlsx)
  if (!file.exists(carb_path)) {
    warning("ocean: carbonate XLSX not found, skipping: ", carb_path)
    return(NULL)
  }

  target_measurements <- tolower(or_default(cfg$carbonate$target_measurements, character()))
  if (length(target_measurements) == 0) return(NULL)

  raw_sheet <- readxl::read_xlsx(carb_path, col_names = FALSE, .name_repair = "minimal")
  header_row <- which(apply(raw_sheet, 1, function(row) {
    any(tolower(as.character(row)) == "parameter")
  }))
  if (length(header_row) == 0) stopf("ocean: no header row containing 'PARAMETER' in %s", carb_path)
  header_row <- header_row[1]

  col_names <- raw_sheet[header_row, 2:ncol(raw_sheet)] %>%
    unlist(use.names = FALSE) %>%
    as.character() %>%
    tolower()
  comments <- raw_sheet[header_row + 3, 2:ncol(raw_sheet)] %>%
    unlist(use.names = FALSE) %>%
    as.character() %>%
    tolower()

  pangea_col <- 1
  col_names[pangea_col] <- "pangea_id"
  sample_col <- which(col_names == "sample material")
  if (length(sample_col) == 0) stopf("ocean: no 'sample material' column in %s", carb_path)
  col_names[sample_col] <- "sample_label"

  comment_suffix <- dplyr::case_when(
    stringr::str_detect(comments, "minimum value") ~ "min",
    stringr::str_detect(comments, "lower quartile") ~ "q1",
    stringr::str_detect(comments, "median value") ~ "q2",
    stringr::str_detect(comments, "upper quartile") ~ "q3",
    stringr::str_detect(comments, "maximum value") ~ "max",
    TRUE ~ NA_character_
  )
  col_names <- ifelse(
    !is.na(comment_suffix),
    paste0(col_names, "_", comment_suffix),
    col_names
  )

  measurement_cols <- which(col_names %in% paste0(target_measurements, "_q2"))
  cols_to_keep <- c(pangea_col, sample_col, measurement_cols)

  data <- raw_sheet[(header_row + 5):nrow(raw_sheet), cols_to_keep + 1]
  colnames(data) <- col_names[cols_to_keep]

  data <- data %>%
    dplyr::inner_join(sample_id_lookup, by = c("pangea_id", "sample_label")) %>%
    dplyr::mutate(sample_id = as.character(.data$analysis_id)) %>%
    dplyr::select("sample_id", dplyr::all_of(col_names[cols_to_keep][!col_names[cols_to_keep] %in% c("pangea_id", "sample_label")]))

  data <- apply_rename_map(data, cfg$rename_map)
  num_cols <- setdiff(names(data), "sample_id")
  data <- data %>% dplyr::mutate(dplyr::across(dplyr::all_of(num_cols), safe_as_numeric))
  as_tibble(data)
}

read_env_ocean <- function(root, cfg) {
  mitag_path <- file.path(root, cfg$input$mitag_tax_profiles_tsv)
  comp_path <- file.path(root, cfg$input$companion_tables_xlsx)
  if (!file.exists(mitag_path)) stopf("ocean: miTAG table not found: %s", mitag_path)
  if (!file.exists(comp_path)) stopf("ocean: CompanionTables not found: %s", comp_path)

  # Build sample_id_lookup from miTAG header only (fast).
  hdr <- names(utils::read.delim(mitag_path, nrows = 0, check.names = FALSE))
  sample_labels <- hdr[grepl("^TARA_", hdr)]
  sample_id_lookup <- tibble(sample_label = sample_labels, analysis_id = seq_along(sample_labels))

  # Join pangea_id from Table W1
  sample_metadata <- readxl::read_xlsx(comp_path, sheet = "Table W1") %>%
    select(
      `Sample label [TARA_station#_environmental-feature_size-fraction]`, # nolint
      `PANGAEA sample identifier` # nolint
    ) %>%
    rename(
      sample_label = `Sample label [TARA_station#_environmental-feature_size-fraction]`,
      pangea_id = `PANGAEA sample identifier`
    )

  sample_id_lookup <- sample_id_lookup %>%
    left_join(sample_metadata, by = "sample_label")

  # Environmental data from Table W8
  w8 <- readxl::read_xlsx(comp_path, sheet = "Table W8")
  # Standardize column names using the config rename_map
  w8 <- apply_rename_map(w8, cfg$rename_map)

  if (!"pangea_id" %in% names(w8)) stopf("ocean: expected pangea_id after renaming")

  envdata <- w8 %>%
    filter(str_starts(as.character(.data$pangea_id), "TARA_")) %>%
    left_join(sample_id_lookup, by = "pangea_id") %>%
    filter(!is.na(.data$analysis_id)) %>%
    mutate(sample_id = as.character(.data$analysis_id)) %>%
    select(-any_of("analysis_id"))

  # Convert numeric columns (everything except ids/labels)
  num_cols <- setdiff(names(envdata), c("sample_id", "pangea_id", "sample_label"))
  envdata <- envdata %>%
    mutate(across(all_of(num_cols), safe_as_numeric))

  # Optional: merge carbonate chemistry from CARB XLSX (q2/median only)
  if (!is.null(cfg$input$carbonate_xlsx) && nzchar(trimws(cfg$input$carbonate_xlsx))) {
    envdata_carbonate <- read_tara_carbonate_xlsx(root, cfg, sample_id_lookup)
    if (!is.null(envdata_carbonate) && nrow(envdata_carbonate) > 0) {
      envdata <- envdata %>%
        left_join(envdata_carbonate, by = "sample_id")
    }
  }

  # Extract the variable catalog from config
  var_catalog <- get_var_catalog(cfg, dataset = "ocean")

  list(
    envdata_raw = envdata,
    var_catalog = var_catalog,
    extras = list(sample_id_lookup = sample_id_lookup)
  )
}

read_env_wwtp <- function(root, cfg) {
  in_path <- file.path(root, cfg$input$path)
  if (!file.exists(in_path)) stopf("wwtp input not found: %s", in_path)

  raw <- readxl::read_excel(in_path, sheet = or_default(cfg$input$sheet, "samples"))
  # Standardize column names using the config rename_map
  raw <- apply_rename_map(raw, cfg$rename_map)
  if (!"sample_id" %in% names(raw)) stopf("wwtp: expected a sample id column mapped to 'sample_id'")

  raw <- raw %>%
    mutate(sample_id = as.character(.data$sample_id)) %>%
    mutate(across(where(is.character), ~ ifelse(.x == "", NA_character_, .x)))

  # Attempt numeric coercion for known numeric env vars (and any columns that already are numeric remain numeric)
  candidate_num <- setdiff(names(raw), c("sample_id"))
  raw <- raw %>% mutate(across(all_of(candidate_num), ~ {
    if (is.numeric(.x)) return(.x)
    safe_as_numeric(.x)
  }))

  # Extract the variable catalog from config
  var_catalog <- get_var_catalog(cfg, dataset = "wwtp")

  list(
    envdata_raw = as_tibble(raw),
    var_catalog = var_catalog,
    extras = list()
  )
}

process_envdata <- function(dataset, root, cfg) {
  if (dataset == "soil") {
    res <- read_env_soil(root, cfg)
  } else if (dataset == "ocean") {
    res <- read_env_ocean(root, cfg)
  } else if (dataset == "wwtp") {
    res <- read_env_wwtp(root, cfg)
  } else {
    stopf("Unknown dataset: %s", dataset)
  }

  env <- res$envdata_raw
  vc <- res$var_catalog
  extras <- or_default(res$extras, list())

  # Filtering
  ranges <- or_default(cfg$filters$ranges, list())

  range_res <- apply_range_filters(env, ranges = ranges)
  env2 <- range_res$df

  base_vars <- or_default(cfg$base_vars, character())

  # Variable selection / pruning
  sel <- choose_variables(
    envdata = env2,
    base_vars = base_vars,
    extra_vars = or_default(cfg$extra_vars, character())
  )

  kept_vars <- sel$kept_vars
  out <- env2 %>% select(all_of(c("sample_id", kept_vars)))

  # Site locations: sample_id, latitude, longitude for same sample set as envdata
  site_locations <- env2 %>%
    filter(.data$sample_id %in% out$sample_id) %>%
    select("sample_id", any_of(c("latitude", "longitude")))
  if (!"latitude" %in% names(site_locations)) site_locations$latitude <- NA_real_
  if (!"longitude" %in% names(site_locations)) site_locations$longitude <- NA_real_
  site_locations <- site_locations %>% select("sample_id", "latitude", "longitude")

  # Reporting
  print_var_catalog(vc)
  cat("\n=== Selection summary ===\n")
  cat(sprintf("- dataset: %s\n", dataset))
  cat(sprintf("- samples: %d\n", nrow(out)))
  if (nrow(range_res$removed) > 0) {
    cat("\nOutlier/range filters removed:\n")
    print(range_res$removed)
  }
  cat(sprintf("\nBase vars (ordered): %s\n", paste(base_vars, collapse = ", ")))
  cat(sprintf("Kept vars: %s\n", paste(kept_vars, collapse = ", ")))

  list(envdata = out, var_catalog = vc, extras = extras, site_locations = site_locations)
}

read_envdata <- function(config_path,
                         root = getwd(),
                         dry_run = FALSE) {
  cfg <- load_config(config_path)

  dataset <- cfg$dataset
  if (is.null(dataset) || !nzchar(dataset)) stopf("Missing dataset: set cfg$dataset in the config file")
  cfg$dataset <- dataset

  res <- process_envdata(dataset = dataset, root = root, cfg = cfg)

  out_dir <- file.path(root, cfg$output$dir)
  write_outputs(out_dir, res$envdata, res$var_catalog, extras = res$extras, site_locations = res$site_locations, dry_run = dry_run)
  res
}

main <- function() {
  args <- parse_cli_args()
  config_path <- or_default(args$config, NULL)
  root <- args$root
  if (is.na(root) || !nzchar(trimws(root))) {
    args_inner <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args_inner, value = TRUE)
    if (length(file_arg) != 1) stop("This script must be run with Rscript so --file= is available.")
    script_path <- normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)
    script_dir  <- dirname(script_path)
    setup_path <- normalizePath(file.path(script_dir, "..", "setup.R"), winslash = "/", mustWork = TRUE)
    source(setup_path)
    root <- REPO_ROOT
  }
  dry_run <- or_default(args$dry_run, FALSE)

  invisible(read_envdata(
    config_path = config_path,
    root = root,
    dry_run = dry_run
  ))
}

if (!interactive()) main()