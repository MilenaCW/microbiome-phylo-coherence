# 00_read_data.R — Step 0: Read env (envdata.csv) and composition, coarse-grain by tax_level, align, scale, save matrices.
# Usage: Rscript 00_read_data.R --config <path> --tax_level <level> --perc_identity <p> [--verbose]

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
})

source(file.path(script_dir, "..", "functions", "CCA_functions.R"))
source(file.path(script_dir, "..", "..", "utility_functions.R"))

option_list <- list(
  make_option(c("--config"), type = "character", default = NULL, help = "Path to CCA config R file (evaluates to a list).", metavar = "FILE"),
  make_option(c("--tax_level"), type = "character", default = "OTU", help = "Taxonomic level: OTU, Domain, Phylum, Class, Order, Family, Genus, Species [default %default]"),
  make_option(c("--perc_identity"), type = "character", default = "0.90", help = "Perc identity (e.g. 0.99) [default %default]"),
  make_option(c("--verbose", "-v"), action = "store_true", default = FALSE, help = "Verbose output.")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$config) || !nzchar(opt$config)) stop("Missing required --config", call. = FALSE)
if (!file.exists(opt$config)) stop("Config file not found: ", opt$config, call. = FALSE)
cfg <- source(opt$config, local = TRUE)$value
if (!is.list(cfg)) stop("Config must evaluate to an R list.", call. = FALSE)

tax_level <- opt$tax_level
perc_identity <- opt$perc_identity
allowed_tax <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "OTU")
if (!tax_level %in% allowed_tax) stop("Invalid --tax_level. Allowed: ", paste(allowed_tax, collapse = ", "), call. = FALSE)

root <- if (exists("REPO_ROOT")) REPO_ROOT else getwd()
data_path <- cfg$data_path
if (!grepl("^/", data_path)) data_path <- file.path(root, data_path)
results_path <- cfg$results_path
if (!grepl("^/", results_path)) results_path <- file.path(root, results_path)

env_file <- cfg$env_file
if (is.null(env_file) || !nzchar(env_file)) env_file <- file.path(data_path, "environmental", "filtered", "envdata.csv")
if (!file.exists(env_file)) stop("envdata.csv not found; run read_data first.\n  Expected: ", env_file, call. = FALSE)
verbose_print(paste("Env file:", env_file), verbose = opt$verbose)

composition_final <- file.path(data_path, "16S", "GG2", perc_identity, "final")
seqtab_fn <- file.path(composition_final, "seqtab.csv")
if (!file.exists(seqtab_fn)) {
  stop("Cleaned up sequence table not found at: ", seqtab_fn, ". Have you run the full GG2 pipeline for this percent identity?", call. = FALSE)
}
verbose_print(paste("Sequence table:", seqtab_fn), verbose = opt$verbose)

data <- load_datasets_cca(env_file, seqtab_fn, read_threshold = 10, verbose = opt$verbose)
env_df <- data$env_df
composition_df <- data$composition_df

if (tax_level != "OTU") {
  taxonomy_fn <- file.path(composition_final, "taxonomy.csv")
  if (!file.exists(taxonomy_fn)) stop("Taxonomy file not found for coarse-graining.", call. = FALSE)
  taxonomy <- read_csv(taxonomy_fn, show_col_types = FALSE)
  needed <- c("Feature_ID", tax_level)
  if (!all(needed %in% names(taxonomy))) stop("Taxonomy must contain Feature_ID and ", tax_level, call. = FALSE)
  otus <- setdiff(names(composition_df), "sample_id")
  missing <- setdiff(otus, taxonomy$Feature_ID)
  if (length(missing) > 0) stop(length(missing), " features in table not in taxonomy. Run GG2/clean first.", call. = FALSE)
  composition_df <- composition_df %>%
    pivot_longer(cols = -sample_id, names_to = "Feature_ID", values_to = "rel_abundance") %>%
    left_join(taxonomy %>% select(Feature_ID, all_of(tax_level)), by = "Feature_ID") %>%
    group_by(!!sym(tax_level), sample_id) %>%
    summarise(rel_abundance = sum(rel_abundance), .groups = "drop") %>%
    pivot_wider(names_from = !!sym(tax_level), values_from = rel_abundance, values_fill = 0)
}
verbose_print(paste0("Compositional Data (",tax_level," level) - Samples: ", nrow(composition_df), " Features: ", ncol(composition_df) - 1), verbose = opt$verbose)

# When using OTU-level data, ensure every composition feature is in the GG2 tree (same final/ dir)
# so downstream step 4 (tree plots) can assume loadings match tree tips.
if (tax_level == "OTU") {
  tree_path <- file.path(composition_final, "tree.nwk")
  if (file.exists(tree_path)) {
    tree <- read.tree(tree_path)
    comp_features <- setdiff(names(composition_df), "sample_id")
    missing_in_tree <- setdiff(comp_features, tree$tip.label)
    if (length(missing_in_tree) > 0) {
      stop(
        length(missing_in_tree), " composition feature(s) are not in the tree. ",
        "Tree and seqtab must come from the same GG2 run (same ", perc_identity, "/final). ",
        "Missing (first 5): ", paste(head(missing_in_tree, 5), collapse = ", ")
      )
    }
  }
}

X <- as.matrix(composition_df %>% select(-sample_id))
Y <- as.matrix(env_df %>% select(-sample_id))
X_scaled <- scale(X)
Y_scaled <- scale(Y)

out_dir <- file.path(results_path, perc_identity, tax_level, "step0_data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(X_scaled, file.path(out_dir, "X_matrix.csv"), row.names = FALSE)
write.csv(Y_scaled, file.path(out_dir, "Y_matrix.csv"), row.names = FALSE)
write.csv(data.frame(sample_id = env_df$sample_id), file.path(out_dir, "sample_id.csv"), row.names = FALSE)
verbose_print(paste("Wrote", out_dir), verbose = opt$verbose)
