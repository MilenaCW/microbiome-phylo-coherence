# 03b_null_tax_shuffle.R — Step 3b: Taxonomy-shuffle null. For each shuffle, randomize
# OTU -> taxonomic lineage assignment (preserving group sizes at every rank), re-coarsen
# from OTU-level composition to the requested --tax_level, then run hyperparameter search
# + loadings as in steps 1-2. Saves null correlations and env loadings (long, with seed)
# for comparison against the real per-tax-level CCA results. Not defined for --tax_level OTU
# (every OTU is its own group there, so shuffling assignment is a no-op).
# Usage: Rscript 03b_null_tax_shuffle.R --config <path> --tax_level <level> --perc_identity <p> [--verbose]

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
  make_option(c("--tax_level"), type = "character", default = "Genus", help = "Taxonomic level to coarse-grain to: Domain, Phylum, Class, Order, Family, Genus, Species [default %default]"),
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

tax_level <- opt$tax_level
perc_identity <- opt$perc_identity
allowed_tax <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
if (!tax_level %in% allowed_tax) {
  stop("Invalid --tax_level. Allowed: ", paste(allowed_tax, collapse = ", "),
       ". (OTU is excluded: every OTU is its own group, so shuffling assignment is a no-op.)", call. = FALSE)
}

root <- if (exists("REPO_ROOT")) REPO_ROOT else getwd()
data_path <- cfg$data_path
if (!grepl("^/", data_path)) data_path <- file.path(root, data_path)
results_path <- cfg$results_path
if (!grepl("^/", results_path)) results_path <- file.path(root, results_path)

env_file <- cfg$env_file
if (is.null(env_file) || !nzchar(env_file)) env_file <- file.path(data_path, "environmental", "filtered", "envdata.csv")
if (!file.exists(env_file)) stop("envdata.csv not found; run read_data first.\n  Expected: ", env_file, call. = FALSE)

composition_final <- file.path(data_path, "16S", "GG2", perc_identity, "final")
seqtab_fn <- file.path(composition_final, "seqtab.csv")
if (!file.exists(seqtab_fn)) {
  stop("Cleaned up sequence table not found at: ", seqtab_fn, ". Have you run the full GG2 pipeline for this percent identity?", call. = FALSE)
}
taxonomy_fn <- file.path(composition_final, "taxonomy.csv")
if (!file.exists(taxonomy_fn)) stop("Taxonomy file not found for coarse-graining.", call. = FALSE)

# Raw OTU-level env + composition (pre-coarsening); each null replicate re-coarsens
# from here under a fresh taxonomy shuffle.
data <- load_datasets_cca(env_file, seqtab_fn, read_threshold = 10, verbose = opt$verbose)
env_df <- data$env_df
composition_df_otu <- data$composition_df
taxonomy <- read_csv(taxonomy_fn, show_col_types = FALSE)

otus <- setdiff(names(composition_df_otu), "sample_id")
missing <- setdiff(otus, taxonomy$Feature_ID)
if (length(missing) > 0) stop(length(missing), " features in table not in taxonomy. Run GG2/clean first.", call. = FALSE)

Y <- as.matrix(env_df %>% select(-sample_id))
Y_scaled <- scale(Y)

hp <- cfg$hyperparam
k <- hp$k
lambda1_range <- hp$lambda1_range
lambda2_range <- hp$lambda2_range
null_cfg <- cfg$null
seeds <- null_cfg$seeds
if (is.null(seeds)) seeds <- seq(0, null_cfg$n_shuffles - 1)

step3b_dir <- file.path(results_path, perc_identity, tax_level, "step3b_null_tax_shuffle")
dir.create(step3b_dir, recursive = TRUE, showWarnings = FALSE)

taxshuffle_corr_list <- list()
taxshuffle_env_list <- list()
idx <- 1L
for (seed in seeds) {
  verbose_print(paste("Taxonomy-shuffle seed", seed, "of", length(seeds)), verbose = opt$verbose)
  taxonomy_shuffled <- shuffle_taxonomy_labels(taxonomy, seed = seed)
  composition_df_shuffled <- coarsen_composition(composition_df_otu, taxonomy_shuffled, tax_level)
  composition_df_shuffled <- composition_df_shuffled %>% arrange(sample_id)
  if (!identical(composition_df_shuffled$sample_id, env_df$sample_id)) {
    stop("Sample order mismatch between shuffled composition and env data.", call. = FALSE)
  }
  X <- as.matrix(composition_df_shuffled %>% select(-sample_id))
  X_scaled <- scale(X)

  results <- k_fold_cv_rcca(X_scaled, Y_scaled, ks = k, lambda1s = lambda1_range, lambda2s = lambda2_range, seed = seed, verbose = FALSE, n_cores = n_cores_arg)
  results_stats <- results %>%
    group_by(k, lambda1, lambda2, canonical_direction) %>%
    summarise(mean_test = mean(cor.test, na.rm = TRUE), .groups = "drop")
  best <- results_stats %>%
    filter(canonical_direction == 1) %>%
    filter(mean_test == max(mean_test)) %>%
    slice(1)
  best_lambda1 <- best$lambda1
  best_lambda2 <- best$lambda2
  loadings <- rcca_loadings(X_scaled, Y_scaled, k = k, lambda1 = best_lambda1, lambda2 = best_lambda2, seed = seed)
  cor_df <- loadings$corr
  aligned_null <- normalize_and_align_loadings(loadings$x_loadings, loadings$y_loadings)
  taxshuffle_corr_list[[idx]] <- cor_df %>% mutate(seed = seed)
  taxshuffle_env_list[[idx]] <- aligned_null$y_loadings %>% mutate(seed = seed)
  idx <- idx + 1L
}
taxshuffle_correlations_per_fold <- bind_rows(taxshuffle_corr_list)
taxshuffle_env_loadings <- bind_rows(taxshuffle_env_list)
write.csv(taxshuffle_correlations_per_fold, file.path(step3b_dir, "taxshuffle_correlations_per_fold.csv"), row.names = FALSE)
write.csv(taxshuffle_env_loadings, file.path(step3b_dir, "taxshuffle_env_loadings.csv"), row.names = FALSE)
verbose_print(paste("Wrote", step3b_dir), verbose = opt$verbose)
