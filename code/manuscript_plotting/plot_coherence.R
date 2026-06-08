# Plot coherence: Pagel's lambda or consenTRAIT tau_D.
#
# Usage (from repo root):
#   Rscript code/manuscript_plotting/plot_coherence.R --dataset soil
#   Rscript code/manuscript_plotting/plot_coherence.R --dataset soil --method consentrait
#   Rscript code/manuscript_plotting/plot_coherence.R --dataset ocean --n_cds 3
#
# Output: manuscript/<ds>/<ds>_coherence_<method>_<perc_short>.pdf (1.75" x 2")

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
})

stopf <- function(...) stop(sprintf(...), call. = FALSE)

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
      help = "Perc identity (e.g. 0.90) [default %default]"
    ),
    optparse::make_option(
      "--n_cds",
      type    = "integer",
      default = NA_integer_,
      help    = "Number of canonical directions to include (plots 1:n_cds); default all"
    ),
    optparse::make_option(
      "--method",
      type    = "character",
      default = "pagel",
      help    = "Coherence method: pagel or consentrait [default %default]"
    )
  )

  p <- optparse::OptionParser(option_list = option_list)
  a <- optparse::parse_args(p, args = argv, positional_arguments = FALSE)

  method <- trimws(tolower(a$method))
  if (!method %in% c("pagel", "consentrait")) {
    stopf("--method must be 'pagel' or 'consentrait', got '%s'", a$method)
  }

  list(
    dataset       = a$dataset,
    perc_identity = a$perc_identity,
    n_cds         = a$n_cds,
    method        = method
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

make_jitter_lookup <- function(keys_df, seed = 1, width = 0.2) {
  keys <- keys_df %>% distinct(.data$canonical_direction, .data$fold)
  set.seed(seed)
  keys %>%
    mutate(jitter = runif(n(), -width, width))
}

load_lambda_true <- function(stats_dir, include_cd) {
  if (length(include_cd) > 0) {
    required_files <- file.path(stats_dir, paste0("true_lambda_cd", include_cd, ".csv"))
    missing <- required_files[!file.exists(required_files)]
    if (length(missing) > 0) {
      stopf("Missing required Pagel lambda files for --include_cd %s: %s. Run 07_pagelsLambda.R for those --direction first.",
            paste(include_cd, collapse = ","), paste(basename(missing), collapse = ", "))
    }
    out <- bind_rows(lapply(required_files, function(f) read_csv(f, show_col_types = FALSE)))
  } else {
    true_files <- list.files(stats_dir, pattern = "^true_lambda_cd[0-9]+\\.csv$", full.names = TRUE)
    if (length(true_files) == 0) {
      stopf("No true_lambda_cd*.csv files found in %s. Run 07_pagelsLambda.R for one or more --direction first.", stats_dir)
    }
    out <- bind_rows(lapply(true_files, function(f) read_csv(f, show_col_types = FALSE)))
  }
  out
}

load_tauD_data <- function(base_dir, include_cd) {
  tau_D_path <- file.path(base_dir, "tau_D.csv")
  null_path  <- file.path(base_dir, "null_tau_D.csv")
  if (!file.exists(tau_D_path)) {
    stopf("tau_D.csv not found in %s. Run 06_consenTRAIT.R first.", base_dir)
  }
  if (!file.exists(null_path)) {
    stopf("null_tau_D.csv not found in %s. Run 06_consenTRAIT.R with null shuffles first.", base_dir)
  }

  obs      <- read_csv(tau_D_path, show_col_types = FALSE)
  null_raw <- read_csv(null_path,  show_col_types = FALSE)

  if (length(include_cd) > 0) {
    obs      <- filter(obs,      canonical_direction %in% include_cd)
    null_raw <- filter(null_raw, canonical_direction %in% include_cd)
    if (nrow(obs) == 0) stopf("No tau_D rows found for --include_cd %s", paste(include_cd, collapse = ","))
  }

  null_summary <- null_raw %>%
    group_by(canonical_direction, fold) %>%
    summarise(null_mean = mean(tau_D), null_sd = sd(tau_D), .groups = "drop")

  list(obs = obs, null = null_summary)
}

median_segments <- function(plot_df, y_col, width = 0.5) {
  plot_df %>%
    group_by(canonical_direction) %>%
    summarise(y_med = median(.data[[y_col]]), .groups = "drop") %>%
    mutate(x_min = canonical_direction - width,
           x_max = canonical_direction + width)
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
  perc       <- args$perc_identity
  include_cd <- if (!is.na(args$n_cds)) seq_len(args$n_cds) else integer(0)
  method     <- args$method

  root        <- get_repo_root()
  point_color <- dataset_colors[ds]
  if (is.na(point_color)) stopf("Could not find color for dataset: %s", ds)

  perc_short <- paste0("p", sub("^[^.]*\\.?", "", perc))
  out_dir    <- file.path(root, "manuscript", ds)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  out_path <- file.path(out_dir, sprintf("%s_coherence_%s_%s.pdf", ds, method, perc_short))

  results_path <- file.path(root, ds, "results", "CCA")

  if (method == "pagel") {
    stats_dir <- file.path(results_path, perc, "OTU", "step7_pagelsLambda", "stats")
    if (!dir.exists(stats_dir)) {
      stopf("Stats directory not found: %s. Run 07_pagelsLambda.R for one or more --direction first.", stats_dir)
    }

    lambda_true <- load_lambda_true(stats_dir, include_cd)
    cd_breaks   <- sort(unique(lambda_true$canonical_direction))

    jitter_lut <- make_jitter_lookup(lambda_true, seed = 1, width = 0.2)
    plot_df <- lambda_true %>%
      left_join(jitter_lut, by = c("canonical_direction", "fold")) %>%
      mutate(x_pos = .data$canonical_direction + .data$jitter)

    med_df <- median_segments(plot_df, "lambda")

    p <- ggplot(plot_df, aes(x = .data$x_pos, y = .data$lambda)) +
      geom_hline(yintercept = 0, linetype = "dotted", color = "grey") +
      geom_segment(
        data = med_df,
        aes(x = x_min, xend = x_max, y = y_med, yend = y_med),
        color = point_color, linewidth = 0.6, inherit.aes = FALSE
      ) +
      geom_point(color = point_color, alpha = 0.5) +
      scale_x_continuous(breaks = cd_breaks) +
      ylim(0, 1) +
      theme_minimal(base_size = 10) +
      labs(x = "", y = "")

  } else {
    base_dir <- file.path(results_path, perc, "OTU", "step6_consentrait")
    if (!dir.exists(base_dir)) {
      stopf("Directory not found: %s. Run 06_consenTRAIT.R first.", base_dir)
    }

    tauD     <- load_tauD_data(base_dir, include_cd)
    obs      <- tauD$obs
    null_sum <- tauD$null
    cd_breaks <- sort(unique(obs$canonical_direction))

    all_keys   <- bind_rows(
      select(obs,      canonical_direction, fold),
      select(null_sum, canonical_direction, fold)
    )
    jitter_lut <- make_jitter_lookup(all_keys, seed = 1, width = 0.2)

    obs_j <- obs %>%
      left_join(jitter_lut, by = c("canonical_direction", "fold")) %>%
      mutate(x_pos = canonical_direction + jitter)

    null_j <- null_sum %>%
      left_join(jitter_lut, by = c("canonical_direction", "fold")) %>%
      mutate(x_pos = canonical_direction + jitter)

    med_df <- median_segments(obs_j, "tau_D")

    p <- ggplot() +
      geom_errorbar(
        data = null_j,
        aes(x = x_pos, ymin = null_mean - null_sd, ymax = null_mean + null_sd),
        color = "grey50", linewidth = 0.4, alpha = 0.5, width = 0
      ) +
      geom_point(
        data = null_j,
        aes(x = x_pos, y = null_mean),
        color = "grey50", alpha = 0.5
      ) +
      geom_segment(
        data = med_df,
        aes(x = x_min, xend = x_max, y = y_med, yend = y_med),
        color = point_color, linewidth = 0.6, inherit.aes = FALSE
      ) +
      geom_point(
        data = obs_j,
        aes(x = x_pos, y = tau_D),
        color = point_color, alpha = 0.5
      ) +
      scale_x_continuous(breaks = cd_breaks) +
      theme_minimal(base_size = 10) +
      labs(x = "", y = "")
  }

  ggsave(out_path, plot = p, width = 1.75, height = 2.00)
  message("Saved: ", out_path)
}

main()
