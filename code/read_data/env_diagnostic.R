# Dataset-agnostic diagnostic plots for environmental variables.
#
# Produces (per dataset):
#  1) Missingness UpSet plot (inclusive intersections)
#  2) Correlogram (GGally ggpairs with colored correlation tiles)
#
# Typical usage:
#   Rscript code/read_data/env_diagnostic.R \
#     --config code/read_data/config/ocean.R
#   Rscript code/read_data/env_diagnostic.R \
#     --config code/read_data/config/soil.R
#   Rscript code/read_data/env_diagnostic.R \
#     --config code/read_data/config/wwtp.R
#
# Notes:
# - Prefers reading processed env tables from cfg$output$dir/{envdata.csv,var_catalog.csv}
# - If envdata.csv is missing, will invoke the env reader via an Rscript subprocess:
#     code/read_data/read_envdata.R

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
})

stopf <- function(...) stop(sprintf(...), call. = FALSE)

or_default <- function(x, y) if (!is.null(x)) x else y

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

is_numericish <- function(x) is.numeric(x) || is.integer(x)

run_env_reader_if_needed <- function(root, config_path, env_csv_path) {
  if (file.exists(env_csv_path)) return(invisible(TRUE))

  reader_path <- file.path(root, "code/read_data/read_envdata.R")
  if (!file.exists(reader_path)) {
    stopf(
      "Processed envdata not found and env reader script missing: %s",
      reader_path
    )
  }

  cmd <- file.path(R.home("bin"), "Rscript")
  args <- c(reader_path, "--config", config_path, "--root", root)
  cat(sprintf("\nProcessed envdata not found; generating it via:\n  %s %s\n", cmd, paste(shQuote(args), collapse = " ")))

  res <- system2(cmd, args = args, stdout = TRUE, stderr = TRUE)
  cat(paste(res, collapse = "\n"), "\n")

  if (!file.exists(env_csv_path)) {
    stopf("Env reader ran but envdata.csv still not found at: %s", env_csv_path)
  }
  invisible(TRUE)
}

plot_upset_missingness <- function(envdata, vars, out_path) {
  if (!requireNamespace("ComplexUpset", quietly = TRUE)) {
    cat("\nNOTE: Skipping UpSet plot because package 'ComplexUpset' is not installed.\n")
    cat("Install with: install.packages('ComplexUpset')\n")
    return(invisible(FALSE))
  }
  suppressPackageStartupMessages({
    library(ComplexUpset)
  })

  membership <- envdata %>%
    select(any_of("sample_id"), all_of(vars)) %>%
    mutate(across(all_of(vars), ~ !is.na(.x))) %>%
    as.data.frame()

  p <- ComplexUpset::upset(
    membership,
    intersect = vars,
    mode = "inclusive_intersection",
    min_size = 1,
    sort_sets = FALSE,
    sort_intersections_by = "cardinality",
    sort_intersections = "descending",
    # set_sizes = FALSE,
    # set_sizes=(
    #     upset_set_size()
    #     + geom_text(aes(label=..count..), hjust=1.1, stat='count')
    # ),
    height_ratio = 1
  ) +
    theme(
      text       = element_text(size = 10),
      axis.text  = element_text(size = 10),
      axis.title = element_text(size = 10),
      strip.text = element_text(size = 10)
    )

  ggsave(out_path, p, width = 6, height = 6, dpi = 300)
  invisible(TRUE)
}

plot_correlogram <- function(envdata, vars, out_path) {
  if (!requireNamespace("GGally", quietly = TRUE)) {
    cat("\nNOTE: Skipping correlogram because package 'GGally' is not installed.\n")
    cat("Install with: install.packages('GGally')\n")
    return(invisible(FALSE))
  }
  suppressPackageStartupMessages({
    library(GGally)
  })

  df <- envdata %>% select(all_of(vars))

  # Ported from tara/code/diagnostic.R
  my_fn <- function(data, mapping, method = "p", use = "pairwise", ...) {
    x <- GGally::eval_data_col(data, mapping$x)
    y <- GGally::eval_data_col(data, mapping$y)
    corr <- stats::cor(x, y, method = method, use = use)

    colFn <- grDevices::colorRampPalette(c("blue", "white", "red"), interpolate = "spline")
    fill <- colFn(100)[findInterval(corr, seq(-1, 1, length = 100))]

    GGally::ggally_cor(data = data, mapping = mapping, color = "black", stars = FALSE, title = "corr", ...) +
      theme_void() +
      theme(
        panel.background = element_rect(fill = fill),
        text = element_text(color = "black", face = "bold", size = 10)
      )
  }

  p <- GGally::ggpairs(
    df,
    upper = list(continuous = my_fn),
    lower = list(continuous = GGally::wrap("points", alpha = 0.5, size = 0.5)),
    diag = list(continuous = GGally::wrap("densityDiag", alpha = 0.5))
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
  invisible(TRUE)
}

parse_cli <- function(argv = commandArgs(trailingOnly = TRUE)) {
  option_list <- list(
    make_option(c("-c", "--config"), type = "character",
                help = "Path to an R config file (list()) from code/read_data/config/"),
    make_option(c("--root"), type = "character",
                default = NA_character_,
                help = "Repo root path (optional; if omitted, uses repo root from code/setup.R)"),
    make_option(c("--outdir"), type = "character",
                default = NA_character_,
                help = "Override output directory for figures/CSVs (default: cfg$output$dir, alongside envdata.csv)"),
    make_option(c("--cutoff"), type = "double",
                default = 3.5,
                help = "Default modified z-score cutoff to report in CSV/plot reference (default: 3.5)")
  )

  p <- OptionParser(option_list = option_list)
  parse_args(p, args = argv, positional_arguments = FALSE)
}

main <- function() {
  args <- parse_cli()
  config_path <- args$config
  root <- args$root
  if (is.na(root) || !nzchar(trimws(root))) {
    args_inner <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args_inner, value = TRUE)
    if (length(file_arg) != 1) {
      stop("This script must be run with Rscript so --file= is available.")
    }
    script_path <- normalizePath(
      sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE
    )
    script_dir  <- dirname(script_path)
    setup_path <- normalizePath(
      file.path(script_dir, "..", "setup.R"), winslash = "/", mustWork = TRUE
    )
    source(setup_path)
    root <- REPO_ROOT
  }

  if (is.null(config_path) || !nzchar(config_path)) stopf("Missing --config")

  # Allow config path relative to root for convenience
  config_abs <- if (file.exists(config_path)) config_path else file.path(root, config_path)
  cfg <- load_config(config_abs)
  dataset <- cfg$dataset
  if (is.null(dataset) || !nzchar(dataset)) {
    stopf("Missing dataset: set cfg$dataset in the config file")
  }

  processed_dir <- file.path(root, cfg$output$dir)
  env_csv <- file.path(processed_dir, "envdata.csv")
  var_catalog_csv <- file.path(processed_dir, "var_catalog.csv")

  vars_mode <- or_default(cfg$vars_mode, basename(cfg$output$dir))

  run_env_reader_if_needed(root = root, config_path = config_abs, env_csv_path = env_csv)

  envdata <- read.csv(env_csv, check.names = FALSE) %>% as_tibble()
  if (!"sample_id" %in% names(envdata)) {
    stopf("envdata.csv must contain a 'sample_id' column: %s", env_csv)
  }
  vars <- setdiff(colnames(envdata), "sample_id")
  if (length(vars) == 0) {
    stopf("No usable numeric variables found (mode=%s) in: %s", vars_mode, env_csv)
  }

  out_dir <- if (!is.na(args$outdir) && nzchar(args$outdir)) {
    if (file.exists(args$outdir) || dir.exists(args$outdir)) args$outdir else file.path(root, args$outdir)
  } else {
    processed_dir
  }
  ensure_dir(out_dir)

  cat("\n=== Env diagnostic ===\n")
  cat(sprintf("- dataset: %s\n", dataset))
  cat(sprintf("- envdata: %s\n", env_csv))
  if (file.exists(var_catalog_csv)) cat(sprintf("- var_catalog: %s\n", var_catalog_csv))
  cat(sprintf("- vars_mode: %s\n", vars_mode))
  cat(sprintf("- n_samples: %d\n", nrow(envdata)))
  cat(sprintf("- n_vars: %d\n", length(vars)))
  cat(sprintf("- out_dir: %s\n", out_dir))

  # Plot 1: UpSet missingness
  plot_upset_missingness(
    envdata = envdata,
    vars = vars,
    out_path = file.path(out_dir, "env_upset_missingness.pdf")
  )

  # Plot 2: Correlogram
  plot_correlogram(
    envdata = envdata,
    vars = vars,
    out_path = file.path(out_dir, "env_correlogram.pdf")
  )

  cat("\nDone.\n")
}

if (!interactive()) main()
