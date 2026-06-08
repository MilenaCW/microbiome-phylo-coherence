# Plot cross-taxonomic environmental loadings (comparative) for manuscript.
# Follows the env loadings comparative logic from CCA/scripts/05_crosstaxonomic_compare.R.
#
# Usage (from repo root):
#   Rscript code/manuscript_plotting/plot_crosstax_env_loadings.R --dataset soil --perc_identity 0.90
#   Rscript code/manuscript_plotting/plot_crosstax_env_loadings.R --dataset ocean --perc_identity 0.99 --n_cds 3
#
# Output: manuscript/<dataset_short>_env_loadings_crosstax_<perc_short>.pdf

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
})

stopf <- function(...) stop(sprintf(...), call. = FALSE)

# Taxonomic level order (match CCA pipeline)
tax_levels_order <- c("OTU", "Species", "Genus", "Family", "Order", "Class", "Phylum")

# Manuscript color palette for taxonomic groups
tax_palette <- c(
  Phylum  = "#decfb2",
  Class   = "#ddc68e",
  Order   = "#dcbd6c",
  Family  = "#dab449",
  Genus   = "#bd863e",
  Species = "#a05532",
  OTU     = "#812727"
)

# Converts env variable display labels to plotmath expressions where needed.
# Chemical species (NO2, PO4, HCO3, CO3) are rendered with proper subscripts and
# superscripts; all other labels pass through as plain upright text.
parse_env_label <- function(labels) {
  chem_map <- c(
    "NO2"  = "NO[2]^{'-'}",
    "PO4"  = "PO[4]^{'3-'}",
    "HCO3" = "HCO[3]^{'-'}",
    "CO3"  = "CO[3]^{'2-'}"
  )
  parts <- lapply(labels, function(lbl) {
    if (lbl %in% names(chem_map)) parse(text = chem_map[[lbl]])[[1]] else lbl
  })
  do.call(expression, parts)
}

parse_cli_args <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (!requireNamespace("optparse", quietly = TRUE)) {
    stopf("R package 'optparse' is required. Install it with: install.packages('optparse')")
  }

  option_list <- list(
    optparse::make_option(
      "--dataset",
      type = "character",
      default = NA_character_,
      help = "Dataset: soil, ocean, or wwtp"
    ),
    optparse::make_option(
      "--perc_identity",
      type = "character",
      default = "0.90",
      help = "Perc identity (e.g. 0.99) [default %default]"
    ),
    optparse::make_option(
      "--n_cds",
      type    = "integer",
      default = NA_integer_,
      help    = "Number of canonical directions to include (plots 1:n_cds); default all"
    )
  )

  p <- optparse::OptionParser(option_list = option_list)
  a <- optparse::parse_args(p, args = argv, positional_arguments = FALSE)

  list(
    dataset       = a$dataset,
    perc_identity = a$perc_identity,
    n_cds         = a$n_cds
  )
}

get_repo_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) != 1) {
    stop("This script must be run with Rscript so --file= is available.", call. = FALSE)
  }
  script_path <- normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)
  script_dir <- dirname(script_path)
  setup_path <- normalizePath(file.path(script_dir, "..", "setup.R"), winslash = "/", mustWork = TRUE)
  source(setup_path, local = TRUE)
  get("REPO_ROOT", envir = .GlobalEnv)
}

main <- function() {
  args <- parse_cli_args()
  if (is.na(args$dataset) || !nzchar(trimws(args$dataset))) {
    stop("Required flag --dataset must be one of: soil, ocean, wwtp", call. = FALSE)
  }
  ds <- trimws(tolower(args$dataset))
  if (!ds %in% c("soil", "ocean", "wwtp")) {
    stopf("Unknown --dataset '%s'. Must be one of: soil, ocean, wwtp", args$dataset)
  }
  perc <- args$perc_identity
  root <- get_repo_root()

  # Load env var order/labels from CCA plot config if available
  env_var_labels <- NULL
  config_path <- file.path(root, "code", "CCA", "config", paste0(ds, ".R"))
  if (file.exists(config_path)) {
    cfg <- source(config_path, local = TRUE)$value
    if (is.list(cfg) && !is.null(cfg$plotting_params$env_var_labels)) {
      env_var_labels <- cfg$plotting_params$env_var_labels
    }
  }

  results_path <- file.path(root, ds, "results", "CCA")
  compare_base <- file.path(results_path, perc)

  # Discover tax levels: dirs under compare_base that have step2_loadings/env_loadings.csv
  all_dirs <- list.dirs(compare_base, full.names = FALSE, recursive = FALSE)
  tax_dirs <- character(0)
  for (d in all_dirs) {
    if (d == "" || d == "crosstax_compare" || d == "Domain") next
    env_file <- file.path(compare_base, d, "step2_loadings", "env_loadings.csv")
    if (file.exists(env_file)) tax_dirs <- c(tax_dirs, d)
  }

  if (length(tax_dirs) == 0) {
    stopf("No tax-level directories with step2_loadings/env_loadings.csv found under %s", compare_base)
  }

  plot_tax_levels <- unique(c(tax_levels_order, tax_dirs))
  plot_tax_levels <- plot_tax_levels[plot_tax_levels %in% tax_dirs]

  # Collect env loadings per tax level
  env_loadings_list <- list()
  for (tax_level in tax_dirs) {
    env_file <- file.path(compare_base, tax_level, "step2_loadings", "env_loadings.csv")
    env_raw <- read_csv(env_file, show_col_types = FALSE)
    env_summary <- env_raw %>%
      group_by(canonical_direction, var) %>%
      summarise(mean = mean(value, na.rm = TRUE), se = sd(value, na.rm = TRUE) / sqrt(dplyr::n()), .groups = "drop") %>%
      mutate(tax_level = tax_level)
    env_summary <- env_summary %>% select(tax_level, canonical_direction, var, mean, se)
    env_loadings_list[[tax_level]] <- env_summary
  }

  all_env_loadings <- bind_rows(env_loadings_list)

  include_cd <- if (!is.na(args$n_cds)) seq_len(args$n_cds) else integer(0)

  # Align env loadings across tax levels (sign flip per canonical_direction), relative to OTU
  aligned_env_loadings <- all_env_loadings
  tax_levels_vec <- unique(all_env_loadings$tax_level)
  directions <- unique(all_env_loadings$canonical_direction)
  ref_tax_level <- if ("OTU" %in% tax_levels_vec) "OTU" else tax_levels_vec[1]

  for (CD in directions) {
    cd_data <- all_env_loadings %>% filter(canonical_direction == CD)
    ref_data <- cd_data %>%
      filter(tax_level == ref_tax_level) %>%
      select(var, mean) %>%
      arrange(var)

    for (tax in tax_levels_vec) {
      if (tax == ref_tax_level) next
      tax_data <- cd_data %>%
        filter(tax_level == tax) %>%
        arrange(var)
      if (!setequal(ref_data$var, tax_data$var)) {
        stopf("Variables do not match between %s and %s for canonical_direction %s", ref_tax_level, tax, CD)
      }
      correlation_og <- cor(ref_data$mean, tax_data$mean)
      correlation_flipped <- cor(ref_data$mean, -tax_data$mean)
      if (correlation_og < correlation_flipped) {
        aligned_env_loadings <- aligned_env_loadings %>%
          mutate(mean = ifelse(canonical_direction == CD & tax_level == tax, -mean, mean))
      }
    }
  }

  # Prepare plot data
  env_plot_data <- aligned_env_loadings
  if (length(include_cd) > 0) {
    env_plot_data <- env_plot_data %>% filter(canonical_direction %in% include_cd)
  }
  env_plot_data <- env_plot_data %>%
    mutate(tax_level = factor(tax_level, levels = plot_tax_levels))

  if (!is.null(env_var_labels) && length(env_var_labels) > 0 && !is.null(names(env_var_labels))) {
    env_plot_data <- env_plot_data %>%
      mutate(var = factor(var, levels = names(env_var_labels), labels = unname(env_var_labels)))
  } else {
    env_plot_data <- env_plot_data %>%
      mutate(var = factor(var, levels = sort(unique(var))))
  }

  # Palette: manuscript colors by taxonomic level name; error if any level has no color
  missing <- plot_tax_levels[!plot_tax_levels %in% names(tax_palette)]
  if (length(missing) > 0) {
    stopf("No color specified for taxonomic level(s): %s. Add them to tax_palette.", paste(missing, collapse = ", "))
  }
  tax_palette_named <- setNames(unname(tax_palette[plot_tax_levels]), plot_tax_levels)

  n_facets <- dplyr::n_distinct(env_plot_data$canonical_direction)
  # n_facet_cols <- min(3L, max(2L, ceiling(sqrt(n_facets))))
  n_facet_cols <- n_facets
  n_facet_rows <- 1
  # n_facet_rows <- ceiling(n_facets / n_facet_cols)
  inch_per_facet_height <- 1.5
  inch_per_facet_width <- 1.5
  env_plot_height <- 0.5 + n_facet_rows * inch_per_facet_height
  env_plot_width <- min(0.25 + n_facet_cols * inch_per_facet_width, 5.25)

  env_plot <- env_plot_data %>%
    ggplot(aes(x = var, y = mean, color = tax_level)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey60") +
    geom_line(aes(group = tax_level), linetype = "solid", linewidth = 0.5, alpha = 0.5) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0) +
    geom_point(size = 1.5) +
    scale_color_manual(values = tax_palette_named) +
    scale_x_discrete(labels = parse_env_label) +
    facet_wrap(
      vars(canonical_direction),
      ncol = n_facet_cols,
      # labeller = labeller(canonical_direction = function(x) paste0("canonical direction ", x))
      labeller = labeller(canonical_direction = function(x) paste0(""))
    ) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      legend.position = "none"
    )

  out_dir <- file.path(root, "manuscript", ds)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  perc_short <- paste0("p", sub("^[^.]*\\.?", "", perc))
  out_path <- file.path(out_dir, sprintf("%s_env_loadings_crosstax_%s.pdf", ds, perc_short))
  ggsave(out_path, plot = env_plot, height = env_plot_height, width = env_plot_width)
  message("Saved: ", out_path)
}

main()
