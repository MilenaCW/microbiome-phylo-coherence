#!/usr/bin/env Rscript
# compute_reference_dists.R
# Computes inter and intra-group cophenetic distances at each taxonomic level
# for one dataset, using the CCA-filtered feature set. Run on the cluster.
#
# Usage (from repo root):
#   Rscript code/reference_distances/compute_reference_dists.R --dataset soil
#   Rscript code/reference_distances/compute_reference_dists.R --dataset ocean --verbose
#   Rscript code/reference_distances/compute_reference_dists.R --dataset soil --n_sample 100000
#
# Output: {dataset}/results/reference_distances/{perc_identity}/dist_{level}.rds

suppressPackageStartupMessages({
  library(optparse)
  library(ape)
  library(dplyr)
  library(readr)
})

stopf <- function(...) stop(sprintf(...), call. = FALSE)

# ---------------------------------------------------------------------------
parse_cli_args <- function(argv = commandArgs(trailingOnly = TRUE)) {
  option_list <- list(
    optparse::make_option(
      "--dataset",
      type    = "character",
      default = NULL,
      help    = "Dataset: soil or ocean [required]"
    ),
    optparse::make_option(
      "--perc_identity",
      type    = "character",
      default = "0.90",
      help    = "Percent identity cutoff for tree/CCA files [default %default]"
    ),
    optparse::make_option(
      "--n_sample",
      type    = "integer",
      default = -1L,
      help    = "Pairs to sample per level per type. -1 = use all pairs [default %default]"
    ),
    optparse::make_option(
      "--tax_levels",
      type    = "character",
      default = "Phylum,Class,Order,Family,Genus",
      help    = "Comma-separated taxonomic levels to compute [default %default]"
    ),
    optparse::make_option(
      "--output_path",
      type    = "character",
      default = NULL,
      help    = "Output directory [default: {dataset}/results/reference_distances/{perc_identity} under repo root]"
    ),
    optparse::make_option(
      "--verbose",
      action  = "store_true",
      default = FALSE,
      help    = "Print progress messages"
    )
  )
  optparse::parse_args(optparse::OptionParser(option_list = option_list), args = argv)
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

vcat <- function(verbose, fmt, ...) {
  if (verbose) cat(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), sprintf(fmt, ...)))
}

# ---------------------------------------------------------------------------
args      <- parse_cli_args()
REPO_ROOT <- get_repo_root()

if (is.null(args$dataset) || !args$dataset %in% c("soil", "ocean")) {
  stopf("--dataset must be 'soil' or 'ocean'. Got: %s",
        if (is.null(args$dataset)) "(missing)" else args$dataset)
}
DATASET  <- args$dataset
PERC_ID  <- args$perc_identity
N_SAMPLE <- args$n_sample
VERBOSE  <- args$verbose
TAX_LEVELS <- trimws(strsplit(args$tax_levels, ",")[[1]])

tree_fn     <- file.path(REPO_ROOT, DATASET, "data/processed_data/16S/GG2",
                         PERC_ID, "final/tree.nwk")
taxonomy_fn <- file.path(REPO_ROOT, DATASET, "data/processed_data/16S/GG2",
                         PERC_ID, "final/taxonomy.csv")
xmat_fn     <- file.path(REPO_ROOT, DATASET, "results/CCA",
                         PERC_ID, "OTU/step0_data/X_matrix.csv")

for (fn in c(tree_fn, taxonomy_fn, xmat_fn)) {
  if (!file.exists(fn)) stopf("File not found: %s", fn)
}

out_dir <- if (!is.null(args$output_path)) {
  args$output_path
} else {
  file.path(REPO_ROOT, DATASET, "results", "reference_distances", PERC_ID)
}
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

vcat(VERBOSE, "Dataset: %s  perc_identity: %s  n_sample: %d\n",
     DATASET, PERC_ID, N_SAMPLE)

# ---------------------------------------------------------------------------
# Step 1: Features from X_matrix columns
vcat(VERBOSE, "Loading features from X_matrix...\n")
features <- colnames(read.csv(xmat_fn, nrows = 0, check.names = FALSE))
cat(sprintf("%d features loaded from X_matrix.\n", length(features)))

# Step 2: Prune tree to CCA features
vcat(VERBOSE, "Loading tree...\n")
tree <- read.tree(tree_fn)

missing <- setdiff(features, tree$tip.label)
if (length(missing) > 0) {
  stopf("%d features in X_matrix not found in tree. First few: %s",
        length(missing), paste(head(missing, 5), collapse = ", "))
}
vcat(VERBOSE, "Pruning tree to %d features...\n", length(features))
tree_sub <- keep.tip(tree, features)
rm(tree); gc()

# Step 3: Cophenetic distance matrix
vcat(VERBOSE, "Computing cophenetic distance matrix...\n")
phylo_dist <- cophenetic(tree_sub)
rm(tree_sub); gc()
vcat(VERBOSE, "Cophenetic matrix size: %d x %d (%.0f MB)\n",
     nrow(phylo_dist), ncol(phylo_dist),
     object.size(phylo_dist) / 1024^2)

# Step 4: Taxonomy
vcat(VERBOSE, "Loading taxonomy...\n")
taxonomy <- read_csv(taxonomy_fn, show_col_types = FALSE) %>%
  filter(Feature_ID %in% features)

missing_tax <- setdiff(features, taxonomy$Feature_ID)
if (length(missing_tax) > 0) {
  stopf("%d features in X_matrix not found in taxonomy.", length(missing_tax))
}

# ---------------------------------------------------------------------------
# Step 5: Compute inter/intra distances per taxonomic level

intra_for_level <- function(tax_df, level_col, dist_mat, n_sample, verbose) {
  groups <- unique(tax_df[[level_col]])
  # Drop groups where label is NA or empty
  groups <- groups[!is.na(groups) & nchar(trimws(as.character(groups))) > 0]
  res <- list()
  for (grp in groups) {
    feats <- tax_df$Feature_ID[tax_df[[level_col]] == grp]
    feats <- feats[!is.na(feats)]
    if (length(feats) < 2) next
    pairs <- expand.grid(otu1 = feats, otu2 = feats, stringsAsFactors = FALSE) %>%
      filter(as.character(otu1) > as.character(otu2))
    if (nrow(pairs) == 0) next
    pairs$dist <- dist_mat[cbind(pairs$otu1, pairs$otu2)]
    res[[length(res) + 1]] <- pairs
  }
  if (length(res) == 0) return(data.frame(dist = numeric(0)))
  out <- bind_rows(res) %>% select(dist)
  if (n_sample > 0 && nrow(out) > n_sample) {
    out <- out[sample(nrow(out), n_sample), , drop = FALSE]
  }
  out
}

inter_for_level <- function(tax_df, level_col, dist_mat, n_sample, verbose) {
  groups <- unique(tax_df[[level_col]])
  groups <- groups[!is.na(groups) & nchar(trimws(as.character(groups))) > 0]
  if (length(groups) < 2) return(data.frame(dist = numeric(0)))
  group_pairs <- combn(groups, 2, simplify = FALSE)
  res <- list()
  for (pair in group_pairs) {
    feats1 <- tax_df$Feature_ID[tax_df[[level_col]] == pair[1]]
    feats2 <- tax_df$Feature_ID[tax_df[[level_col]] == pair[2]]
    feats1 <- feats1[!is.na(feats1)]
    feats2 <- feats2[!is.na(feats2)]
    if (length(feats1) == 0 || length(feats2) == 0) next
    pairs  <- expand.grid(otu1 = feats1, otu2 = feats2, stringsAsFactors = FALSE)
    pairs$dist <- dist_mat[cbind(pairs$otu1, pairs$otu2)]
    res[[length(res) + 1]] <- pairs
  }
  if (length(res) == 0) return(data.frame(dist = numeric(0)))
  out <- bind_rows(res) %>% select(dist)
  if (n_sample > 0 && nrow(out) > n_sample) {
    out <- out[sample(nrow(out), n_sample), , drop = FALSE]
  }
  out
}

for (lvl in TAX_LEVELS) {
  vcat(VERBOSE, "Computing level: %s\n", lvl)
  if (!lvl %in% colnames(taxonomy)) {
    cat(sprintf("WARNING: Level '%s' not found in taxonomy columns, skipping.\n", lvl))
    next
  }

  intra <- intra_for_level(taxonomy, lvl, phylo_dist, N_SAMPLE, VERBOSE) %>%
    mutate(level = lvl, type = "intra")
  inter <- inter_for_level(taxonomy, lvl, phylo_dist, N_SAMPLE, VERBOSE) %>%
    mutate(level = lvl, type = "inter")

  cat(sprintf("  %s: intra n=%d (mean=%.4f), inter n=%d (mean=%.4f)\n",
              lvl,
              nrow(intra), mean(intra$dist, na.rm = TRUE),
              nrow(inter), mean(inter$dist, na.rm = TRUE)))

  level_fn <- file.path(out_dir, sprintf("dist_%s.rds", lvl))
  saveRDS(bind_rows(intra, inter) %>% mutate(dataset = DATASET), level_fn)
  cat(sprintf("  Saved: %s\n", level_fn))
  rm(intra, inter); gc()
}
cat(sprintf("All levels saved to: %s\n", out_dir))
