#!/usr/bin/env Rscript
# 02_loadings.R — Step 2: Run rcca_loadings with best hyperparams; normalize and align loadings; save correlations_per_fold and loadings (long).
# Usage: Rscript 02_loadings.R --config <path> --tax_level <level> --perc_identity <p> [--verbose] [--debug]

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
})
source(file.path(script_dir, "..", "functions", "CCA_functions.R"))
source(file.path(script_dir, "..", "..", "utility_functions.R"))

option_list <- list(
  make_option(c("--config"), type = "character", default = NULL, help = "Path to CCA config R file.", metavar = "FILE"),
  make_option(c("--tax_level"), type = "character", default = "OTU", help = "Taxonomic level [default %default]"),
  make_option(c("--perc_identity"), type = "character", default = "0.90", help = "Perc identity [default %default]"),
  make_option(c("--verbose", "-v"), action = "store_true", default = FALSE, help = "Verbose output."),
  make_option(c("--debug"), action = "store_true", default = FALSE, help = "Save intermediate files (untouched loadings, fold_indices, projections, flips) [default %default].")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$config) || !nzchar(opt$config)) stop("Missing required --config", call. = FALSE)
if (!file.exists(opt$config)) stop("Config file not found: ", opt$config, call. = FALSE)
cfg <- source(opt$config, local = TRUE)$value
if (!is.list(cfg)) stop("Config must evaluate to an R list.", call. = FALSE)

root <- if (exists("REPO_ROOT")) REPO_ROOT else getwd()
results_path <- cfg$results_path
if (!grepl("^/", results_path)) results_path <- file.path(root, results_path)
step0_dir <- file.path(results_path, opt$perc_identity, opt$tax_level, "step0_data")
step1_dir <- file.path(results_path, opt$perc_identity, opt$tax_level, "step1_hyperparam")
step2_dir <- file.path(results_path, opt$perc_identity, opt$tax_level, "step2_loadings")
dir.create(step2_dir, recursive = TRUE, showWarnings = FALSE)

X <- as.matrix(read_csv(file.path(step0_dir, "X_matrix.csv"), show_col_types = FALSE))
Y <- as.matrix(read_csv(file.path(step0_dir, "Y_matrix.csv"), show_col_types = FALSE))
best <- read.csv(file.path(step1_dir, "best_hyperparam.csv"))
k <- best$k[1]
lambda1 <- best$lambda1[1]
lambda2 <- best$lambda2[1]
seed <- cfg$hyperparam$seed

verbose_print("Extracting loadings...", verbose = opt$verbose)
loadings <- rcca_loadings(X, Y, k = k, lambda1 = lambda1, lambda2 = lambda2, seed = seed)
x_loadings_df <- loadings$x_loadings
y_loadings_df <- loadings$y_loadings
if (opt$debug) {
  write.csv(x_loadings_df, file.path(step2_dir, "abundance_loadings_untouched.csv"), row.names = FALSE)
  write.csv(y_loadings_df, file.path(step2_dir, "env_loadings_untouched.csv"), row.names = FALSE)
  fold_indices <- loadings$fold_indices
  write.table(tibble(fold = fold_indices), file.path(step2_dir, "fold_indices.csv"), row.names = FALSE)
  projections <- loadings$projections
  write.csv(projections, file.path(step2_dir, "projections.csv"), row.names = FALSE)
}
corr_df <- loadings$corr
write.csv(corr_df, file.path(step2_dir, "correlations_per_fold.csv"), row.names = FALSE)

verbose_print("Normalizing and aligning the canonical directions...", verbose = opt$verbose)
aligned <- normalize_and_align_loadings(x_loadings_df, y_loadings_df)
x_loadings_normalized <- aligned$x_loadings
y_loadings_normalized <- aligned$y_loadings
y_loadings_original   <- aligned$y_prealign  # normalised, pre-alignment; used by plot_alignment_assessment
if (opt$debug) {
  flips_long <- tibble::enframe(
    aligned$flipped_folds,
    name = "canonical_direction",
    value = "fold"
  ) %>%
    tidyr::unnest_longer(fold) %>%
    dplyr::mutate(canonical_direction = as.integer(canonical_direction))
  readr::write_csv(flips_long, file.path(step2_dir, "flips.csv"))
}

write.csv(x_loadings_normalized, file.path(step2_dir, "abundance_loadings.csv"), row.names = FALSE)
write.csv(y_loadings_normalized, file.path(step2_dir, "env_loadings.csv"), row.names = FALSE)

# Alignment assessment: original vs aligned loadings (env only; colored by fold)
plot_alignment_assessment(y_loadings_original, y_loadings_normalized, dir_path = step2_dir)

verbose_print(paste("Wrote", step2_dir), verbose = opt$verbose)
