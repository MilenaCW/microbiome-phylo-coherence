# 03_null_distribution.R — Step 3: For each shuffle, run hyperparameter search, then calculate loadings at best params; save null correlations and loadings (long, with seed).
# Usage: Rscript 03_null_distribution.R --config <path> --tax_level <level> --perc_identity <p> [--verbose]

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
  make_option(c("--n_cores"), type = "integer", default = NA, help = "Number of parallel workers (default: all available cores).", metavar = "N")
)
opt <- parse_args(OptionParser(option_list = option_list))
n_cores_arg <- if (is.na(opt$n_cores)) NULL else opt$n_cores
if (is.null(opt$config) || !nzchar(opt$config)) stop("Missing required --config", call. = FALSE)
if (!file.exists(opt$config)) stop("Config file not found: ", opt$config, call. = FALSE)
cfg <- source(opt$config, local = TRUE)$value
if (!is.list(cfg)) stop("Config must evaluate to an R list.", call. = FALSE)

root <- if (exists("REPO_ROOT")) REPO_ROOT else getwd()
results_path <- cfg$results_path
if (!grepl("^/", results_path)) results_path <- file.path(root, results_path)
step0_dir <- file.path(results_path, opt$perc_identity, opt$tax_level, "step0_data")
step3_dir <- file.path(results_path, opt$perc_identity, opt$tax_level, "step3_null")
dir.create(step3_dir, recursive = TRUE, showWarnings = FALSE)

X <- as.matrix(read_csv(file.path(step0_dir, "X_matrix.csv"), show_col_types = FALSE))
Y <- as.matrix(read_csv(file.path(step0_dir, "Y_matrix.csv"), show_col_types = FALSE))
hp <- cfg$hyperparam
k <- hp$k
lambda1_range <- hp$lambda1_range
lambda2_range <- hp$lambda2_range
null_cfg <- cfg$null
seeds <- null_cfg$seeds
if (is.null(seeds)) seeds <- seq(0, null_cfg$n_shuffles - 1)

null_corr_list <- list()
null_env_list <- list()
idx <- 1L
for (seed in seeds) {
  verbose_print(paste("Null seed", seed, "of", length(seeds)), verbose = opt$verbose)
  X_shuffled <- shuffle_sample(X, seed = seed)
  results <- k_fold_cv_rcca(X_shuffled, Y, ks = k, lambda1s = lambda1_range, lambda2s = lambda2_range, seed = seed, verbose = FALSE, n_cores = n_cores_arg)
  results_stats <- results %>%
    group_by(k, lambda1, lambda2, canonical_direction) %>%
    summarise(mean_test = mean(cor.test, na.rm = TRUE), .groups = "drop")
  best <- results_stats %>%
    filter(canonical_direction == 1) %>%
    filter(mean_test == max(mean_test)) %>%
    slice(1)
  best_lambda1 <- best$lambda1
  best_lambda2 <- best$lambda2
  loadings <- rcca_loadings(X_shuffled, Y, k = k, lambda1 = best_lambda1, lambda2 = best_lambda2, seed = seed)
  cor_df       <- loadings$corr
  aligned_null <- normalize_and_align_loadings(loadings$x_loadings, loadings$y_loadings)
  null_corr_list[[idx]] <- cor_df %>% mutate(seed = seed)
  null_env_list[[idx]]  <- aligned_null$y_loadings %>% mutate(seed = seed)
  idx <- idx + 1L
}
null_correlations_per_fold <- bind_rows(null_corr_list)
null_env_loadings <- bind_rows(null_env_list)
write.csv(null_correlations_per_fold, file.path(step3_dir, "null_correlations_per_fold.csv"), row.names = FALSE)
write.csv(null_env_loadings, file.path(step3_dir, "null_env_loadings.csv"), row.names = FALSE)
verbose_print(paste("Wrote", step3_dir), verbose = opt$verbose)
