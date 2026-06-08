# SI_coherence.R — SI figure: phylogenetic coherence for soil and ocean.
#
# Reads from (consentrait):
#   {dataset}/results/CCA/{perc}/OTU/step6_consentrait/tau_D.csv
#   {dataset}/results/CCA/{perc}/OTU/step6_consentrait/null_tau_D.csv
#
# Reads from (pagel):
#   {dataset}/results/CCA/{perc}/OTU/step7_pagelsLambda/stats/true_lambda_cd{N}.csv
#
# Output (one PDF):
#   manuscript/SI/SI_coherence_{method}_p{perc_short}.pdf
#
# Usage (from repo root):
#   Rscript code/manuscript_plotting/SI_coherence.R
#   Rscript code/manuscript_plotting/SI_coherence.R --method pagel
#   Rscript code/manuscript_plotting/SI_coherence.R --perc_identity 0.90

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(patchwork)
  library(optparse)
  library(withr)
})

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
DATASET_COLORS <- c(soil = "#8f723d", ocean = "#90afa7")
NULL_COLOR     <- "grey50"
SIG_CDS        <- list(soil = 1:4, ocean = 1:3)
BASE_SIZE      <- 10

# ---------------------------------------------------------------------------
# Root detection
# ---------------------------------------------------------------------------
get_repo_root <- function() {
  args     <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) != 1)
    stop("Run with Rscript so --file= is available.", call. = FALSE)
  script_path <- normalizePath(sub("^--file=", "", file_arg),
                               winslash = "/", mustWork = TRUE)
  setup_path  <- normalizePath(
    file.path(dirname(script_path), "..", "setup.R"),
    winslash = "/", mustWork = TRUE
  )
  source(setup_path, local = TRUE)
  get("REPO_ROOT", envir = .GlobalEnv)
}

# ---------------------------------------------------------------------------
# CLI args
# ---------------------------------------------------------------------------
parse_cli_args <- function(argv = commandArgs(trailingOnly = TRUE)) {
  option_list <- list(
    make_option("--perc_identity", type = "character", default = "0.90",
                help = "Percent identity [default %default]"),
    make_option("--method", type = "character", default = "consentrait",
                help = "Coherence method: consentrait or pagel [default %default]"),
    make_option("--root", type = "character", default = NA,
                help = "Repo root (inferred from script location if omitted)")
  )
  a <- parse_args(OptionParser(option_list = option_list), args = argv)
  method <- trimws(tolower(a$method))
  if (!method %in% c("consentrait", "pagel"))
    stop(sprintf("--method must be 'consentrait' or 'pagel', got '%s'", a$method), call. = FALSE)
  list(perc = a$perc_identity, method = method, root = trimws(a$root))
}

# ---------------------------------------------------------------------------
# Data loaders
# ---------------------------------------------------------------------------

load_tau_D <- function(root, dataset, perc) {
  base <- file.path(root, dataset, "results", "CCA", perc,
                    "OTU", "step6_consentrait")
  sig  <- SIG_CDS[[dataset]]

  obs <- read_csv(file.path(base, "tau_D.csv"), show_col_types = FALSE) %>%
    filter(canonical_direction %in% sig)

  null_raw <- read_csv(file.path(base, "null_tau_D.csv"), show_col_types = FALSE) %>%
    filter(canonical_direction %in% sig)

  null_summary <- null_raw %>%
    group_by(canonical_direction, fold) %>%
    summarise(null_mean = mean(tau_D), null_sd = sd(tau_D), .groups = "drop")

  list(obs = obs, null = null_summary)
}

load_lambda_true <- function(stats_dir, sig_cds) {
  required_files <- file.path(stats_dir, paste0("true_lambda_cd", sig_cds, ".csv"))
  missing <- required_files[!file.exists(required_files)]
  if (length(missing) > 0) {
    stop(sprintf("Missing Pagel lambda files for CDs %s: %s. Run 07_pagelsLambda.R first.",
                 paste(sig_cds, collapse = ","), paste(basename(missing), collapse = ", ")),
         call. = FALSE)
  }
  bind_rows(lapply(required_files, function(f) read_csv(f, show_col_types = FALSE)))
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

make_jitter_lookup <- function(cd_vals, fold_vals, seed = 1L, width = 0.2) {
  keys <- expand.grid(
    cd   = sort(unique(cd_vals)),
    fold = sort(unique(fold_vals))
  )
  keys$jitter <- withr::with_seed(seed, runif(nrow(keys), -width, width))
  keys
}

median_segments <- function(df, y_col, width = 0.5) {
  df %>%
    group_by(canonical_direction) %>%
    summarise(y_med = median(.data[[y_col]]), .groups = "drop") %>%
    mutate(x_min = canonical_direction - width,
           x_max = canonical_direction + width)
}

# ---------------------------------------------------------------------------
# Panel builders
# ---------------------------------------------------------------------------

make_panel1 <- function(tau_D_obs, null_summary, dataset_color, sig_cds) {
  jitter_lut <- make_jitter_lookup(
    cd_vals   = c(tau_D_obs$canonical_direction, null_summary$canonical_direction),
    fold_vals = c(tau_D_obs$fold,               null_summary$fold)
  )

  obs <- tau_D_obs %>%
    left_join(jitter_lut, by = c("canonical_direction" = "cd", "fold")) %>%
    mutate(x_plot = canonical_direction + jitter)

  null <- null_summary %>%
    left_join(jitter_lut, by = c("canonical_direction" = "cd", "fold")) %>%
    mutate(x_plot = canonical_direction + jitter)

  med_df <- median_segments(obs, "tau_D")

  ggplot() +
    geom_point(
      data = null,
      aes(x = x_plot, y = null_mean),
      color = NULL_COLOR,
      alpha = 0.5
    ) +
    geom_errorbar(
      data = null,
      aes(x = x_plot, y = null_mean,
          ymin = null_mean - null_sd, ymax = null_mean + null_sd),
      color = NULL_COLOR,
      linewidth = 0.5,
      alpha = 0.5
    ) +
    geom_segment(
      data = med_df,
      aes(x = x_min, xend = x_max, y = y_med, yend = y_med),
      color = dataset_color, linewidth = 0.6, inherit.aes = FALSE
    ) +
    geom_point(
      data  = obs,
      aes(x = x_plot, y = tau_D),
      color = dataset_color,
      alpha = 0.5
    ) +
    scale_x_continuous(breaks = sig_cds, labels = sig_cds) +
    labs(
      x = "canonical direction (CD)",
      y = expression("weighted mean genetic depth (" * tau[D] * ")")
    ) +
    theme_minimal(base_size = BASE_SIZE) +
    theme(panel.grid.minor = element_blank())
}

make_pagel_panel <- function(lambda_df, dataset_color, sig_cds) {
  jitter_lut <- make_jitter_lookup(
    cd_vals   = lambda_df$canonical_direction,
    fold_vals = lambda_df$fold
  )

  plot_df <- lambda_df %>%
    left_join(jitter_lut, by = c("canonical_direction" = "cd", "fold")) %>%
    mutate(x_plot = canonical_direction + jitter)

  med_df <- median_segments(plot_df, "lambda")

  ggplot(plot_df, aes(x = x_plot, y = lambda)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey") +
    geom_segment(
      data = med_df,
      aes(x = x_min, xend = x_max, y = y_med, yend = y_med),
      color = dataset_color, linewidth = 0.6, inherit.aes = FALSE
    ) +
    geom_point(color = dataset_color, alpha = 0.5) +
    scale_x_continuous(breaks = sig_cds, labels = sig_cds) +
    ylim(0, 1) +
    labs(
      x = "canonical direction (CD)",
      y = expression("Pagel's " * lambda)
    ) +
    theme_minimal(base_size = BASE_SIZE) +
    theme(panel.grid.minor = element_blank())
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
main <- function() {
  cli    <- parse_cli_args()
  root   <- if (!is.na(cli$root)) cli$root else get_repo_root()
  perc   <- cli$perc
  method <- cli$method

  perc_short <- gsub("\\.", "", perc)
  out_dir    <- file.path(root, "manuscript", "SI")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_path   <- file.path(out_dir,
                           paste0("SI_coherence_", method, "_p", perc_short, ".pdf"))

  message("Loading data...")

  if (method == "consentrait") {
    soil_tauD  <- load_tau_D(root, "soil",  perc)
    ocean_tauD <- load_tau_D(root, "ocean", perc)

    message("Building panels...")
    p_soil  <- make_panel1(soil_tauD$obs,  soil_tauD$null,
                            DATASET_COLORS["soil"],  SIG_CDS[["soil"]])
    p_ocean <- make_panel1(ocean_tauD$obs, ocean_tauD$null,
                            DATASET_COLORS["ocean"], SIG_CDS[["ocean"]])

    fig <- (p_soil + p_ocean) &
      coord_cartesian(ylim = c(
        min(c(soil_tauD$null$null_mean, ocean_tauD$null$null_mean), na.rm = TRUE) * 0.9,
        max(c(soil_tauD$obs$tau_D,     ocean_tauD$obs$tau_D),      na.rm = TRUE) * 1.1
      ))
  } else {
    soil_stats_dir  <- file.path(root, "soil",  "results", "CCA", perc,
                                  "OTU", "step7_pagelsLambda", "stats")
    ocean_stats_dir <- file.path(root, "ocean", "results", "CCA", perc,
                                  "OTU", "step7_pagelsLambda", "stats")

    lambda_soil  <- load_lambda_true(soil_stats_dir,  SIG_CDS[["soil"]])
    lambda_ocean <- load_lambda_true(ocean_stats_dir, SIG_CDS[["ocean"]])

    message("Building panels...")
    p_soil  <- make_pagel_panel(lambda_soil,  DATASET_COLORS["soil"],  SIG_CDS[["soil"]])
    p_ocean <- make_pagel_panel(lambda_ocean, DATASET_COLORS["ocean"], SIG_CDS[["ocean"]])

    fig <- p_soil + p_ocean
  }

  message("Saving: ", out_path)
  ggsave(out_path, fig, width = 6.5, height = 3, device = "pdf")
  message("Done.")
}

if (!interactive()) main()
