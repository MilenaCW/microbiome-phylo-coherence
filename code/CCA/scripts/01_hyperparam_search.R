# 01_hyperparam_search.R — Step 1: k-fold CV over lambda1/lambda2, select best by mean test cor (canonical_direction 1), save CV results and plot.
# Usage: Rscript 01_hyperparam_search.R --config <path> --tax_level <level> --perc_identity <p> [--verbose]

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
source(file.path(dirname(script_path), "..", "..", "utility_functions.R"))
source(file.path(dirname(script_path), "..", "functions", "CCA_functions.R"))

option_list <- list(
  make_option(c("--config"), type = "character", default = NULL, help = "Path to CCA config R file.", metavar = "FILE"),
  make_option(c("--tax_level"), type = "character", default = "OTU", help = "Taxonomic level [default %default]"),
  make_option(c("--perc_identity"), type = "character", default = "0.90", help = "Perc identity [default %default]"),
  make_option(c("--verbose", "-v"), action = "store_true", default = FALSE, help = "Verbose output."),
  make_option(c("--n_cores"), type = "integer", default = NA, help = "Number of parallel workers (default: all available cores).", metavar = "N"),
  make_option(c("--run_label"), type = "character", default = "", 
            help = "Optional label to namespace outputs (e.g. benchmark scenario).")
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
step1_dir <- file.path(results_path, opt$perc_identity, opt$tax_level, "step1_hyperparam")
run_label <- opt$run_label
if (nzchar(run_label)) {
  step1_dir <- file.path(step1_dir, run_label)
}
dir.create(step1_dir, recursive = TRUE, showWarnings = FALSE)

X <- as.matrix(read_csv(file.path(step0_dir, "X_matrix.csv"), show_col_types = FALSE))
Y <- as.matrix(read_csv(file.path(step0_dir, "Y_matrix.csv"), show_col_types = FALSE))
hp <- cfg$hyperparam
k <- hp$k
lambda1_range <- hp$lambda1_range
lambda2_range <- hp$lambda2_range
seed <- hp$seed

verbose_print("Running k-fold CV hyperparameter search...", verbose = opt$verbose)
results <- k_fold_cv_rcca(X, Y, ks = k, lambda1s = lambda1_range, lambda2s = lambda2_range, seed = seed, verbose = opt$verbose, n_cores = n_cores_arg)
write.csv(results, file.path(step1_dir, "cv_results.csv"), row.names = FALSE)

results_stats <- results %>%
  group_by(k, lambda1, lambda2, canonical_direction) %>%
  summarise(
    mean_test = mean(cor.test, na.rm = TRUE),
    sd_test = sd(cor.test, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )
best <- results_stats %>%
  filter(canonical_direction == 1) %>%
  filter(mean_test == max(mean_test))
best <- best[1, ]
best_hyperparam <- data.frame(k = best$k, lambda1 = best$lambda1, lambda2 = best$lambda2)
write.csv(best_hyperparam, file.path(step1_dir, "best_hyperparam.csv"), row.names = FALSE)
verbose_print(paste("Best: k =", best$k, "lambda1 =", best$lambda1, "lambda2 =", best$lambda2), verbose = opt$verbose)

group_label <- paste0(opt$tax_level, " (", as.numeric(opt$perc_identity) * 100, "% ID)")
plot_cv_lambda(results, best$lambda1, step1_dir, group_label)
verbose_print(paste("Wrote", step1_dir), verbose = opt$verbose)