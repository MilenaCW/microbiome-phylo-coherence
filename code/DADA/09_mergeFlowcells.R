# Script to merge flowcell sequence tables and rank ASVs by mean relative abundance.

# Source repo setup (sets REPO_ROOT and working directory)
# Resolve this script's directory (works when run via Rscript)
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) != 1) stop("This script must be run with Rscript so --file= is available.")
script_path <- normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)
script_dir  <- dirname(script_path)
setup_path <- normalizePath(file.path(script_dir, "..", "setup.R"),
                            winslash = "/", mustWork = TRUE)
source(setup_path)

library(dada2)
library(optparse)
library(Biostrings)
library(tidyverse)
source('./code/utility_functions.R')

# Define command-line options
option_list <- list(
  make_option(c("--config"), action = "store", default = NULL,
              help = "Path to an R config file that evaluates to a list (required)."),
  make_option(c("--verbose"), action = "store_true", default = TRUE,
              help = "Print verbose output (Default: TRUE).")
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

load_config <- function(config_path) {
  if (is.null(config_path) || !nzchar(config_path)) {
    stop("Missing required --config <path>", call. = FALSE)
  }
  if (!file.exists(config_path)) {
    stop(paste0("Config file not found: ", config_path), call. = FALSE)
  }
  cfg <- source(config_path, local = TRUE)$value
  if (!is.list(cfg)) {
    stop("Config must evaluate to an R list", call. = FALSE)
  }
  cfg
}

cfg <- load_config(opt$config)

output_directory <- cfg$output$directory
verbose <- if (!is.null(cfg$verbose)) cfg$verbose else opt$verbose

if (is.null(output_directory)) {
  stop("Config missing required field: output$directory", call. = FALSE)
}

find_flowcells <- function(base_dir) {
  flowcells <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE)
  flowcells <- flowcells[flowcells != "merged"]
  flowcells
}

merge_flow_cells <- function(flow_cells_to_merge) {
  print("Merging flow cells...")
  seqtab_list <- list()

  for (i in seq_along(flow_cells_to_merge)) {
    flow_cell <- flow_cells_to_merge[i]
    seqtab_path <- file.path(output_directory, flow_cell, "07_removeChimera",
                             paste0(flow_cell, "_seqtab_noChimeras.rds"))

    if (file.exists(seqtab_path)) {
      cat("Reading sequence table for flow cell:", flow_cell, "\n")
      seqtab_list[[i]] <- readRDS(seqtab_path)
    } else {
      cat("Warning: Sequence table not found for flow cell:", flow_cell, "\n")
    }
  }

  seqtab_list <- seqtab_list[!sapply(seqtab_list, is.null)]
  if (length(seqtab_list) == 0) {
    stop("No valid sequence tables found to merge")
  }

  cat("Merging", length(seqtab_list), "sequence tables...\n")
  merged_seqtab <- mergeSequenceTables(tables = seqtab_list)
  cat("Merged sequence table has", nrow(merged_seqtab), "samples and",
      ncol(merged_seqtab), "ASVs\n")
  merged_seqtab
}

rank_asvs_by_mean_rel_abundance <- function(seqtab) {
  if (ncol(seqtab) == 0) {
    stop("Merged sequence table has no ASVs to rank", call. = FALSE)
  }
  totals <- rowSums(seqtab)
  rel_abund <- sweep(seqtab, 1, totals, "/")
  mean_rel <- colMeans(rel_abund)
  order_idx <- order(-mean_rel, colnames(seqtab)) # order by descending mean relative abundance and break ties by alpha order of sequences
  list(order_idx = order_idx, mean_rel = mean_rel)
}

verbose_print("Finding flowcell directories in output directory...", verbose)
flowcells <- find_flowcells(output_directory)
verbose_print("Flowcell directories found in output_directory:", verbose)
cat(paste0("  ", flowcells, collapse = "\n"), "\n")

flow_cells_to_merge <- c()
for (flow_cell in flowcells) {
  seqtab_path <- file.path(output_directory, flow_cell, "07_removeChimera",
                           paste0(flow_cell, "_seqtab_noChimeras.rds"))
  if (file.exists(seqtab_path)) {
    cat("Flow cell", flow_cell, "successfully completed chimera removal, merging...\n")
    flow_cells_to_merge <- c(flow_cells_to_merge, flow_cell)
  } else {
    cat("Flow cell", flow_cell, "missing chimera-removed seqtab, skipping...\n")
  }
}

merged_dir <- file.path(output_directory, "merged")
dir.create(merged_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(flow_cells_to_merge, file.path(merged_dir, "flow_cells_to_merge.csv"),
          row.names = FALSE)

merged_seqtab <- merge_flow_cells(flow_cells_to_merge)

ranking <- rank_asvs_by_mean_rel_abundance(merged_seqtab)
ordered_sequences <- colnames(merged_seqtab)[ranking$order_idx]
asv_names <- paste0("ASV_", seq_along(ordered_sequences))

# Convert the column names to ASV_
merged_seqtab <- merged_seqtab[, ranking$order_idx, drop = FALSE]
colnames(merged_seqtab) <- asv_names
# sort the samples in ascending order (numeric if possible, otherwise the characters)
merged_seqtab <- as.data.frame(merged_seqtab) %>%
  rownames_to_column(var = "sample_id") %>%
  mutate(
    .sample_id_num = suppressWarnings(as.numeric(sample_id))
  ) %>%
  arrange(
    if (all(!is.na(.sample_id_num))) .sample_id_num else sample_id
  ) %>%
  dplyr::select(-.sample_id_num)

write.csv(merged_seqtab, file.path(merged_dir, "merged_seqtab_asv.csv"),
          row.names = FALSE)

asv_set <- DNAStringSet(ordered_sequences)
names(asv_set) <- asv_names
writeXStringSet(asv_set, file.path(merged_dir, "merged_asv_sequences.fna"))

print("Merged sequence table CSV and ASV sequences FNA saved")
