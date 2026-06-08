# Plot test correlation vs null for one dataset at a given perc_identity and tax_level.
#
# Usage (from repo root):
#   Rscript code/manuscript_plotting/plot_corr_performance.R --dataset soil --perc_identity 0.90 --tax_level OTU
#   Rscript code/manuscript_plotting/plot_corr_performance.R --dataset ocean --perc_identity 0.90 --tax_level Genus --root /path/to/repo
#
# Output: manuscript/<dataset_short>_corr_performance_<perc_short>_<tax_level>.pdf

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
})

stopf <- function(...) stop(sprintf(...), call. = FALSE)

# Dataset colors (null is always grey in the plot)
dataset_colors <- c(
  "soil"  = "#8f723d",
  "ocean" = "#90afa7",
  "wwtp"  = "#606161"
)

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
      "--tax_level",
      type = "character",
      default = "OTU",
      help = "Taxonomic level (e.g. OTU, Genus) [default %default]"
    ),
    optparse::make_option(
      "--root",
      type = "character",
      default = NA_character_,
      help = "Repo root path (optional; if omitted, uses code/setup.R)"
    )
  )

  p <- optparse::OptionParser(option_list = option_list)
  a <- optparse::parse_args(p, args = argv, positional_arguments = FALSE)

  list(
    dataset = a$dataset,
    perc_identity = a$perc_identity,
    tax_level = a$tax_level,
    root = a$root
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
  tax  <- args$tax_level

  root <- args$root
  if (is.na(root) || !nzchar(trimws(root))) {
    root <- get_repo_root()
  } else {
    root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  }

  # CCA results path: <dataset>/results/CCA (same layout as CCA config results_path)
  results_path <- file.path(root, ds, "results", "CCA")
  step2_dir <- file.path(results_path, perc, tax, "step2_loadings")
  step3_dir <- file.path(results_path, perc, tax, "step3_null")

  corr_file <- file.path(step2_dir, "correlations_per_fold.csv")
  null_file <- file.path(step3_dir, "null_correlations_per_fold.csv")
  if (!file.exists(corr_file)) {
    stopf("Correlations file not found: %s\nRun CCA step 2 for this dataset/perc_identity/tax_level first.", corr_file)
  }
  if (!file.exists(null_file)) {
    stopf("Null correlations file not found: %s\nRun CCA step 3 for this dataset/perc_identity/tax_level first.", null_file)
  }

  corr_per_fold <- read.csv(corr_file)
  true_stats <- corr_per_fold %>%
    group_by(.data$canonical_direction) %>%
    summarise(mean = mean(.data$test, na.rm = TRUE), sd = sd(.data$test, na.rm = TRUE), .groups = "drop")

  null_corr <- read.csv(null_file)
  null_cd1 <- null_corr %>% filter(.data$canonical_direction == 1)
  null_per_seed <- null_cd1 %>% group_by(.data$seed) %>% summarise(mean_test = mean(.data$test, na.rm = TRUE), .groups = "drop")
  null_mean <- mean(null_per_seed$mean_test)
  null_se <- sd(null_per_seed$mean_test) / sqrt(nrow(null_per_seed))

  point_color <- dataset_colors[ds]
  if (is.na(point_color)) stopf("Could not find color for dataset: %s", ds)

  n_cd <- max(true_stats$canonical_direction)
  y_raw_min <- min(true_stats$mean - true_stats$sd, na.rm = TRUE)
  y_raw_max <- max(true_stats$mean + true_stats$sd, na.rm = TRUE)

  y_min <- ceiling(min(0, y_raw_min) / 0.25) * 0.25
  y_max <- floor(max(0, y_raw_max) / 0.25) * 0.25

  p <- ggplot(true_stats, aes(x = factor(.data$canonical_direction), y = .data$mean)) +
    geom_hline(yintercept = null_mean, linetype = "dashed", color = "grey", linewidth = 0.7) +
    annotate("rect", xmin = 0.5, xmax = n_cd + 0.5, ymin = null_mean - null_se, ymax = null_mean + null_se,
             fill = "grey", alpha = 0.15) +
    geom_point(color = point_color) +
    geom_errorbar(aes(ymin = .data$mean - .data$sd, ymax = .data$mean + .data$sd), width = 0, linewidth = 0.8, color = point_color) +
    labs(
      x = "",
      y = ""
      # title = sprintf("%s: true vs null (CD1)", dataset_label),
      # subtitle = sprintf("perc_identity = %s, tax_level = %s", perc, tax)
    ) +
    scale_y_continuous(
      limits = c(min(0,y_raw_min), max(1.0,y_raw_max)),
      breaks = seq(y_min, y_max, by = 0.25)
    ) +
    theme_minimal(base_size = 10)

  out_dir <- file.path(root, "manuscript", ds)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  perc_short <- paste0("p", sub("^[^.]*\\.?", "", perc))
  out_path <- file.path(out_dir, sprintf("%s_corr_performance_%s_%s.pdf", ds, perc_short, tax))
  ggsave(out_path, plot = p, width = 2.75, height = 1.75)
  message("Saved: ", out_path)
}

main()
