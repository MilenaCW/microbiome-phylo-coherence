# SI figures S1, S2, S3 — environmental variable diagnostics for manuscript.
#
#   S1: Correlogram of filtered soil environmental variables
#   S2: Missingness structure of full ocean environmental variable set (UpSet)
#   S3: Correlogram of core ocean environmental variables
#
# Usage (from repo root):
#   Rscript code/manuscript_plotting/SI_env_diagnostics.R --figure S1
#   Rscript code/manuscript_plotting/SI_env_diagnostics.R --figure S2
#   Rscript code/manuscript_plotting/SI_env_diagnostics.R --figure S3
#
# Output: manuscript/SI/SI_{soil_correlogram,ocean_upset_missingness,ocean_correlogram}.pdf
#
# Reads from already-processed envdata.csv files; run code/read_data/read_envdata.R first
# if those CSVs are missing.

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

# ---------------------------------------------------------------------------
# Figure configs
# ---------------------------------------------------------------------------

# Column labels for correlograms use plotmath strings so GGally can render
# chemical notation via labeller = "label_parsed". Plain-text labels (depth,
# temperature, pH, etc.) are just strings and parse harmlessly as symbols.
# The UpSet (S2) uses plain ASCII because ComplexUpset's labeller does not
# parse plotmath expressions.

FIGURE_CONFIGS <- list(

  S1 = list(
    env_csv   = "soil/data/processed_data/environmental/filtered/envdata.csv",
    plot_type = "correlogram",
    out_name  = "SI_soil_correlogram.pdf",
    col_labels = c(
      ph        = "pH",
      soil_c    = "C",
      soil_n    = "N",
      soil_p    = "P",
      clay_silt = "clay+silt"
    )
  ),

  S2 = list(
    env_csv   = "ocean/data/processed_data/environmental/full/envdata.csv",
    plot_type = "upset",
    out_name  = "SI_ocean_upset_missingness.pdf",
    # Plotmath strings: ComplexUpset's labeller can return an expression object,
    # which ggplot2 renders as plotmath on the y-axis (scale_y_discrete).
    label_map = c(
      depth                      = "'depth'",
      temperature                = "'temperature'",
      salinity                   = "'salinity'",
      oxygen                     = "'oxygen'",
      no2                        = "NO[2]^{'-'}",
      po4                        = "PO[4]^{'3-'}",
      no2no3                     = "NO[2]^{'-'}*' + '*NO[3]^{'-'}",
      si                         = "'Si'",
      pH                         = "'pH'",
      co2                        = "CO[2]",
      co2_pp                     = "'pCO'[2]",
      co2_f                      = "'fCO'[2]",
      hco3                       = "HCO[3]^{'-'}",
      co3                        = "CO[3]^{'2-'}",
      carbon_total               = "'DIC'",
      alkalinity_total           = "'TA'",
      calcite_saturation_state   = "Omega~'calcite'",
      aragonite_saturation_state = "Omega~'aragonite'"
    )
  ),

  S3 = list(
    env_csv   = "ocean/data/processed_data/environmental/partial-filtered/envdata.csv",
    plot_type = "correlogram",
    out_name  = "SI_ocean_correlogram.pdf",
    col_labels = c(
      depth       = "depth",
      temperature = "temperature",
      salinity    = "salinity",
      oxygen      = "oxygen",
      no2         = "NO[2]^{'-'}",
      po4         = "PO[4]^{'3-'}",
      no2no3      = "NO[2]^{'-'}*' + '*NO[3]^{'-'}",
      si          = "Si"
    )
  )

)

# ---------------------------------------------------------------------------
# Plot functions
# ---------------------------------------------------------------------------

plot_correlogram <- function(envdata, vars, col_labels, out_path) {
  if (!requireNamespace("GGally", quietly = TRUE)) {
    stop("Package 'GGally' is required: install.packages('GGally')", call. = FALSE)
  }
  suppressPackageStartupMessages(library(GGally))

  df <- envdata %>% select(all_of(vars))

  display_labels <- vapply(vars, function(v) {
    if (!is.null(col_labels) && v %in% names(col_labels)) unname(col_labels[[v]]) else v
  }, character(1))

  my_fn <- function(data, mapping, method = "p", use = "pairwise", ...) {
    x    <- GGally::eval_data_col(data, mapping$x)
    y    <- GGally::eval_data_col(data, mapping$y)
    corr <- stats::cor(x, y, method = method, use = use)
    colFn <- grDevices::colorRampPalette(c("blue", "white", "red"), interpolate = "spline")
    fill  <- colFn(100)[findInterval(corr, seq(-1, 1, length = 100))]
    GGally::ggally_cor(data = data, mapping = mapping, color = "black", stars = FALSE, title = "corr", ...) +
      theme_void() +
      theme(
        panel.background = element_rect(fill = fill),
        text = element_text(color = "black", face = "bold", size = 10)
      )
  }

  p <- GGally::ggpairs(
    df,
    columnLabels = display_labels,
    labeller     = "label_parsed",
    upper = list(continuous = my_fn),
    lower = list(continuous = GGally::wrap("points", alpha = 0.5, size = 0.5)),
    diag  = list(continuous = GGally::wrap("densityDiag", alpha = 0.5))
  ) +
    theme(
      text        = element_text(size = 10),
      axis.text   = element_text(size = 10),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.title  = element_text(size = 10),
      strip.text  = element_text(size = 10)
    )

  num_vars <- length(vars)
  ggsave(out_path, p, width = num_vars, height = num_vars, dpi = 300)
  message("Saved: ", out_path)
  invisible(TRUE)
}

plot_upset_missingness <- function(envdata, vars, label_map, out_path) {
  if (!requireNamespace("ComplexUpset", quietly = TRUE)) {
    stop("Package 'ComplexUpset' is required: install.packages('ComplexUpset')", call. = FALSE)
  }
  suppressPackageStartupMessages(library(ComplexUpset))

  membership <- envdata %>%
    select(any_of("sample_id"), all_of(vars)) %>%
    mutate(across(all_of(vars), ~ !is.na(.x))) %>%
    as.data.frame()

  # ComplexUpset calls labeller(non_sanitized_labels[sets]) and passes the result
  # to scale_y_discrete(labels = ...). Returning an expression object causes
  # ggplot2 to render the labels as plotmath, giving proper chemical notation.
  expr_labeller <- function(x) {
    parts <- lapply(x, function(lbl) {
      if (!is.null(label_map) && lbl %in% names(label_map)) {
        parse(text = label_map[[lbl]])[[1]]
      } else {
        lbl
      }
    })
    do.call(expression, parts)
  }

  p <- ComplexUpset::upset(
    membership,
    intersect             = vars,
    labeller              = expr_labeller,
    mode                  = "inclusive_intersection",
    min_size              = 1,
    sort_sets             = FALSE,
    sort_intersections_by = "cardinality",
    sort_intersections    = "descending",
    height_ratio          = 1
  ) +
    theme(
      text       = element_text(size = 10),
      axis.text  = element_text(size = 10),
      axis.title = element_text(size = 10),
      strip.text = element_text(size = 10)
    )

  ggsave(out_path, p, width = 8, height = 8, dpi = 300)
  message("Saved: ", out_path)
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

get_repo_root <- function() {
  args     <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) != 1) stop("This script must be run with Rscript.", call. = FALSE)
  script_path <- normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)
  setup_path  <- normalizePath(
    file.path(dirname(script_path), "..", "setup.R"), winslash = "/", mustWork = TRUE
  )
  source(setup_path, local = TRUE)
  get("REPO_ROOT", envir = .GlobalEnv)
}

parse_cli <- function(argv = commandArgs(trailingOnly = TRUE)) {
  option_list <- list(
    make_option(c("-f", "--figure"), type = "character", default = NA_character_,
                help = "Figure to produce: S1, S2, or S3")
  )
  parse_args(OptionParser(option_list = option_list), args = argv, positional_arguments = FALSE)
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

main <- function() {
  args <- parse_cli()
  fig  <- trimws(toupper(or_null(args$figure, "")))
  if (!nzchar(fig) || !fig %in% c("S1", "S2", "S3")) {
    stop("Required flag --figure must be one of: S1, S2, S3", call. = FALSE)
  }

  root    <- get_repo_root()
  cfg     <- FIGURE_CONFIGS[[fig]]
  env_csv <- file.path(root, cfg$env_csv)

  if (!file.exists(env_csv)) {
    stop(sprintf(
      "Processed envdata not found: %s\nRun code/read_data/read_envdata.R first.",
      env_csv
    ), call. = FALSE)
  }

  envdata <- read.csv(env_csv, check.names = FALSE) %>% as_tibble()
  if (!"sample_id" %in% names(envdata)) {
    stop(sprintf("envdata.csv must contain a 'sample_id' column: %s", env_csv), call. = FALSE)
  }

  out_dir <- file.path(root, "manuscript", "SI")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_path <- file.path(out_dir, cfg$out_name)

  cat(sprintf("\n=== SI figure %s ===\n", fig))
  cat(sprintf("- envdata : %s  (N=%d)\n", env_csv, nrow(envdata)))
  cat(sprintf("- output  : %s\n", out_path))

  if (cfg$plot_type == "correlogram") {
    vars <- intersect(names(cfg$col_labels), setdiff(names(envdata), "sample_id"))
    if (length(vars) == 0) stop("No matching variables found in envdata.", call. = FALSE)
    cat(sprintf("- vars    : %s\n", paste(vars, collapse = ", ")))
    plot_correlogram(envdata, vars = vars, col_labels = cfg$col_labels, out_path = out_path)

  } else if (cfg$plot_type == "upset") {
    vars <- intersect(names(cfg$label_map), setdiff(names(envdata), "sample_id"))
    if (length(vars) == 0) stop("No matching variables found in envdata.", call. = FALSE)
    cat(sprintf("- vars    : %s\n", paste(vars, collapse = ", ")))
    plot_upset_missingness(envdata, vars = vars, label_map = cfg$label_map, out_path = out_path)
  }
}

or_null <- function(x, default) if (!is.null(x) && !is.na(x)) x else default

if (!interactive()) main()
