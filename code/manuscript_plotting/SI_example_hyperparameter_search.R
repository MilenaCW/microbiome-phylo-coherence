# SI_example_hyperparameter_search.R — Recreate the cross-validation lambda plot
# (canonical direction 1 only) from saved cv_results.csv / best_hyperparam.csv.
# No re-analysis needed. Saves a 3x3 in PDF to manuscript/.
#
# Usage (from repo root):
#   Rscript code/manuscript_plotting/SI_example_hyperparameter_search.R
#   Rscript code/manuscript_plotting/SI_example_hyperparameter_search.R \
#     --dataset soil --tax_level OTU --perc_identity 0.90

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

stopf <- function(...) stop(sprintf(...), call. = FALSE)

parse_cli_args <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (!requireNamespace("optparse", quietly = TRUE)) {
    stopf("R package 'optparse' is required. Install it with: install.packages('optparse')")
  }
  option_list <- list(
    optparse::make_option("--dataset",       type = "character", default = "soil",
                          help = "Dataset: soil, ocean, or wwtp [default %default]"),
    optparse::make_option("--tax_level",     type = "character", default = "OTU",
                          help = "Taxonomic level [default %default]"),
    optparse::make_option("--perc_identity", type = "character", default = "0.90",
                          help = "Perc identity (e.g. 0.90) [default %default]"),
    optparse::make_option("--root",          type = "character", default = NA_character_,
                          help = "Repo root path (optional; if omitted, uses code/setup.R)")
  )
  p <- optparse::OptionParser(option_list = option_list)
  a <- optparse::parse_args(p, args = argv, positional_arguments = FALSE)
  list(dataset = a$dataset, tax_level = a$tax_level,
       perc_identity = a$perc_identity, root = a$root)
}

get_repo_root <- function() {
  args     <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) != 1) {
    stop("This script must be run with Rscript so --file= is available.", call. = FALSE)
  }
  script_path <- normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)
  script_dir  <- dirname(script_path)
  setup_path  <- normalizePath(file.path(script_dir, "..", "setup.R"), winslash = "/", mustWork = TRUE)
  source(setup_path, local = TRUE)
  get("REPO_ROOT", envir = .GlobalEnv)
}

main <- function() {
  args <- parse_cli_args()
  ds   <- trimws(tolower(args$dataset))
  perc <- args$perc_identity
  tax  <- args$tax_level

  root <- args$root
  if (is.na(root) || !nzchar(trimws(root))) {
    root <- get_repo_root()
  } else {
    root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  }

  step1_dir <- file.path(root, ds, "results", "CCA", perc, tax, "step1_hyperparam")
  cv_file   <- file.path(step1_dir, "cv_results.csv")
  bp_file   <- file.path(step1_dir, "best_hyperparam.csv")

  if (!file.exists(cv_file)) stopf("CV results not found: %s", cv_file)
  if (!file.exists(bp_file)) stopf("Best hyperparam not found: %s", bp_file)

  results <- read.csv(cv_file)
  best    <- read.csv(bp_file)

  # Canonical direction 1 is what drove lambda selection
  results <- results[results$canonical_direction == 1, ]

  summ <- results %>%
    group_by(k, lambda1, lambda2, canonical_direction) %>%
    summarise(
      mean_train = mean(cor.train, na.rm = TRUE),
      sd_train   = sd(cor.train,   na.rm = TRUE),
      mean_test  = mean(cor.test,  na.rm = TRUE),
      sd_test    = sd(cor.test,    na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols      = c(mean_train, sd_train, mean_test, sd_test),
      names_to  = c("stat", "type"),
      names_sep = "_",
      values_to = "value"
    ) %>%
    pivot_wider(names_from = stat, values_from = value) %>%
    mutate(type = factor(type, levels = c("test", "train")))

  ymin <- min(summ$mean - summ$sd, 0, na.rm = TRUE)
  ymax <- max(summ$mean + summ$sd, 1, na.rm = TRUE)

  p <- ggplot(summ, aes(x = log10(lambda1), y = mean)) +
    geom_vline(xintercept = log10(best$lambda1), linetype = "dashed",
               color = "red", linewidth = 0.8) +
    geom_line(color = "black") +
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd),
                alpha = 0.15, fill = "darkblue") +
    facet_grid(
      rows = vars(type),
      cols = vars(canonical_direction),
      labeller = labeller(
        type                = as_labeller(c(test = "Test", train = "Train")),
        canonical_direction = as_labeller(c("1" = "CD 1 (used for selection)"))
      )
    ) +
    labs(x = expression(log[10](lambda[1])), y = NULL) +
    # scale_y_continuous(limits = c(ymin, ymax), breaks = c(0.0, 0.5, 1.0)) +
    theme_minimal() +
    theme(strip.text.y = element_blank())

  out_dir    <- file.path(root, "manuscript", "SI")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  perc_short <- sub("^0\\.", "", perc)   # "0.90" -> "90"
  out_file   <- file.path(out_dir, sprintf("SI_%s_cv_lambda_%s_%s.pdf", ds, perc_short, tax))
  ggsave(out_file, plot = p, width = 2.75, height = 3, units = "in")
  message("Saved: ", out_file)
}

main()
