#!/usr/bin/env Rscript
# 07_pagelsLambda.R — Pagel's lambda (phylolm) per fold for canonical directions 1..N; optional null (shuffle traits on tree).
# Usage: Rscript 07_pagelsLambda.R --config <path> --perc_identity <p> --n_cds <N> [--n_shuffles 0] [--seed 0] [--verbose]

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) != 1) stop("Run with Rscript so --file= is available.")
script_path <- normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)
script_dir <- dirname(script_path)
setup_path <- normalizePath(file.path(script_dir, "..", "..", "setup.R"), winslash = "/", mustWork = TRUE)
source(setup_path)

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(readr)
  library(ape)
  library(phylolm)
})
source(file.path(script_dir, "..", "functions", "CCA_functions.R"))
source(file.path(script_dir, "..", "..", "utility_functions.R"))

option_list <- list(
  make_option(c("--config"), type = "character", default = NULL, help = "Path to CCA config R file.", metavar = "FILE"),
  make_option(c("--perc_identity"), type = "character", default = "0.90", help = "Perc identity [default %default]."),
  make_option(c("--n_cds"), type = "integer", default = NULL, help = "Number of canonical directions to process (1..N); required.", metavar = "N"),
  make_option(c("--verbose", "-v"), action = "store_true", default = FALSE, help = "Verbose output."),
  make_option(c("--n_shuffles"), type = "integer", default = 0L, help = "Number of null shuffles; 0 = skip null [default %default]."),
  make_option(c("--seed"), type = "integer", default = 0L, help = "Seed before null shuffles [default %default].")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$config) || !nzchar(opt$config)) stop("Missing required --config", call. = FALSE)
if (is.null(opt$perc_identity) || !nzchar(opt$perc_identity)) stop("Missing required --perc_identity", call. = FALSE)
if (is.null(opt$n_cds) || is.na(opt$n_cds) || opt$n_cds < 1L) stop("Missing required --n_cds (positive integer)", call. = FALSE)
if (!file.exists(opt$config)) stop("Config file not found: ", opt$config, call. = FALSE)
cfg <- source(opt$config, local = TRUE)$value
if (!is.list(cfg)) stop("Config must evaluate to an R list.", call. = FALSE)

n_shuffles <- if (is.null(opt$n_shuffles) || is.na(opt$n_shuffles) || opt$n_shuffles <= 0L) 0L else as.integer(opt$n_shuffles)

root <- if (exists("REPO_ROOT")) REPO_ROOT else getwd()
results_path <- cfg$results_path
if (!grepl("^/", results_path)) results_path <- file.path(root, results_path)
data_path <- cfg$data_path
if (!grepl("^/", data_path)) data_path <- file.path(root, data_path)

step2_dir <- file.path(results_path, opt$perc_identity, "OTU", "step2_loadings")
abundance_loadings_file <- file.path(step2_dir, "abundance_loadings.csv")
composition_final <- file.path(data_path, "16S", "GG2", opt$perc_identity, "final")
tree_file <- file.path(composition_final, "tree.nwk")
out_base <- file.path(results_path, opt$perc_identity, "OTU", "step7_pagelsLambda")
stats_dir <- file.path(out_base, "stats")

# --- Compute: true lambda (and optionally null) for directions 1..n_cds ---
verbose_print(paste("Starting Pagel's lambda (step 7) for", opt$n_cds, "direction(s)..."), verbose = opt$verbose)
verbose_print(paste("Parameters: n_cds =", opt$n_cds, ", n_shuffles =", n_shuffles, ", seed =", opt$seed, ", output =", out_base), verbose = opt$verbose)
if (!file.exists(abundance_loadings_file)) stop("Abundance loadings file not found at: ", abundance_loadings_file, call. = FALSE)
if (!file.exists(tree_file)) stop("Tree file not found at: ", tree_file, call. = FALSE)
loadings <- read_csv(abundance_loadings_file, show_col_types = FALSE)
tree <- read.tree(tree_file)
if (n_shuffles > 0L) set.seed(opt$seed)

safe_phylolm_lambda <- function(tr, template, y_vec, verbose_flag = FALSE) {
  template$y <- y_vec
  out <- tryCatch(
    phylolm::phylolm(y ~ 1, data = template, phy = tr, model = "lambda"),
    error = function(e) {
      if (verbose_flag) warning("phylolm failed: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(out)) return(NA_real_)
  as.numeric(out$optpar)
}

fold_to_y <- function(cell_data, tip_order) {
  x <- setNames(cell_data$value, cell_data$var)
  y <- as.numeric(x[tip_order])
  if (any(!is.finite(y))) stop("Non-finite or missing values after alignment to tree tips.", call. = FALSE)
  y
}

for (direction in seq_len(opt$n_cds)) {
  verbose_print(paste("Processing direction", direction, "of", opt$n_cds), verbose = opt$verbose)
  loadings_dir <- loadings %>% filter(canonical_direction == direction)
  if (nrow(loadings_dir) == 0) {
    warning("No loadings for canonical_direction ", direction, "; skipping.", call. = FALSE)
    next
  }
  otus <- sort(unique(loadings_dir$var))
  missing_in_tree <- setdiff(otus, tree$tip.label)
  if (length(missing_in_tree) > 0) stop("Some OTUs not in tree: ", length(missing_in_tree), " missing.", call. = FALSE)
  tree_trimmed <- keep.tip(tree, otus)
  tip_order <- tree_trimmed$tip.label
  n_tips <- length(tip_order)
  df_template <- data.frame(y = rep(NA_real_, n_tips), row.names = tip_order)
  loadings_dir <- loadings_dir %>% filter(var %in% tip_order)
  loadings_by_fold <- split(loadings_dir, loadings_dir$fold)
  folds <- sort(as.integer(names(loadings_by_fold)))
  verbose_print(paste("Trimmed tree to", n_tips, "tips; folds:", paste(folds, collapse = ", ")), verbose = opt$verbose)

  true_rows <- vector("list", length(folds))
  null_rows <- list()

  for (i in seq_along(folds)) {
    f <- folds[i]
    cell_data <- loadings_by_fold[[as.character(f)]]
    y <- fold_to_y(cell_data, tip_order)
    lambda_true <- safe_phylolm_lambda(tree_trimmed, df_template, y, opt$verbose)
    true_rows[[i]] <- tibble(canonical_direction = direction, fold = f, lambda = lambda_true)
    if (n_shuffles > 0L) {
      fold_null <- vector("list", n_shuffles)
      for (s in seq_len(n_shuffles)) {
        y_shuf <- sample(y)
        lambda_null <- safe_phylolm_lambda(tree_trimmed, df_template, y_shuf, opt$verbose)
        fold_null[[s]] <- tibble(canonical_direction = direction, fold = f, lambda = lambda_null, shuffle = s)
      }
      null_rows[[length(null_rows) + 1L]] <- bind_rows(fold_null)
    }
    verbose_print(paste("Completed lambda calculation for fold", f), verbose = opt$verbose)
  }
  true_df <- bind_rows(true_rows)
  dir.create(stats_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv(true_df, file.path(stats_dir, paste0("true_lambda_cd", direction, ".csv")))
  verbose_print(paste("Wrote", file.path(stats_dir, paste0("true_lambda_cd", direction, ".csv"))), verbose = opt$verbose)
  if (n_shuffles > 0L && length(null_rows) > 0L) {
    null_df <- bind_rows(null_rows)
    write_csv(null_df, file.path(stats_dir, paste0("null_lambda_cd", direction, ".csv")))
    verbose_print(paste("Wrote", file.path(stats_dir, paste0("null_lambda_cd", direction, ".csv"))), verbose = opt$verbose)
  }
  verbose_print(paste("Direction", direction, "complete."), verbose = opt$verbose)
}

# --- Plot: aggregate across all processed directions ---
verbose_print("=== Plotting Pagel's lambda across canonical directions ===", verbose = opt$verbose)
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(stats_dir)) stop("Stats directory not found: ", stats_dir, ". No directions were processed.", call. = FALSE)
true_files <- list.files(stats_dir, pattern = "^true_lambda_cd[0-9]+\\.csv$", full.names = TRUE)
if (length(true_files) == 0) stop("No true_lambda_cd*.csv files found in ", stats_dir, call. = FALSE)
true_all <- bind_rows(lapply(true_files, function(f) read_csv(f, show_col_types = FALSE)))
null_files <- list.files(stats_dir, pattern = "^null_lambda_cd[0-9]+\\.csv$", full.names = TRUE)
null_all <- if (length(null_files) > 0) bind_rows(lapply(null_files, function(f) read_csv(f, show_col_types = FALSE))) else tibble(canonical_direction = integer(), fold = integer(), lambda = numeric(), shuffle = integer())
verbose_print(paste("Read", nrow(true_all), "true lambda rows,", nrow(null_all), "null lambda rows"), verbose = opt$verbose)
cd_breaks <- sort(unique(true_all$canonical_direction))
p <- ggplot()
if (nrow(null_all) > 0) {
  null_summ <- null_all %>%
    group_by(canonical_direction, fold) %>%
    summarise(mean_lambda = mean(lambda, na.rm = TRUE), sd_lambda = sd(lambda, na.rm = TRUE), .groups = "drop")
  p <- p +
    geom_pointrange(
      data = null_summ,
      aes(x = canonical_direction, y = mean_lambda, ymin = mean_lambda - sd_lambda, ymax = mean_lambda + sd_lambda),
      color = "grey", alpha = 0.5, position = position_jitter(width = 0.2, height = 0, seed = 1)
    )
}
p <- p +
  geom_point(
    data = true_all,
    aes(x = canonical_direction, y = lambda),
    color = "black", alpha = 0.5, position = position_jitter(width = 0.2, height = 0, seed = 1)
  ) +
  scale_x_continuous(breaks = cd_breaks) +
  ylim(0,1) +
  theme_minimal() +
  labs(x = "canonical direction", y = expression(lambda))
ggsave(file.path(out_base, "pagel_lambda.jpeg"), plot = p, width = 5, height = 5)
verbose_print(paste("Saved", file.path(out_base, "pagel_lambda.jpeg")), verbose = opt$verbose)
verbose_print("Pagel's lambda (step 7) complete.", verbose = opt$verbose)
