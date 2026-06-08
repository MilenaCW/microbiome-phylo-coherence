# Script to do the sixth step of DADA2 pipeline: build the sequence table

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

library(tidyverse)
library(dada2)
library(optparse)
source('./code/utility_functions.R')

# Define command-line options
option_list <- list(
  make_option(c("--config"), action = "store", default = NULL,
              help = "Path to an R config file that evaluates to a list (required)."),
  make_option(c("--flowcell"), action = "store", default = NULL,
              help = "Build the sequence table for a given flowcell (Default: NULL, all flow cells)."),
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
seqtab_cfg <- cfg$sequence_table
verbose <- if (!is.null(cfg$verbose)) cfg$verbose else opt$verbose

if (is.null(output_directory)) {
  stop("Config missing required field: output$directory", call. = FALSE)
}
if (is.null(seqtab_cfg) || is.null(seqtab_cfg$length_range)) {
  stop("Config missing required field: sequence_table$length_range", call. = FALSE)
}

build_sequence_table <- function(flowcell, verbose = TRUE) {
  start_time <- Sys.time()

  merged_dir <- file.path(output_directory, flowcell, "05_merge")
  if (!dir.exists(merged_dir)) {
    stop(paste0("Merged directory not found: ", merged_dir), call. = FALSE)
  }

  # get the merged pairs
  mergers <- readRDS(file.path(merged_dir, paste0(flowcell, "_mergers.rds")))

  # create directory for the sequence table
  sequence_table_dir <- file.path(output_directory, flowcell, "06_buildTable")
  dir.create(sequence_table_dir, recursive = TRUE, showWarnings = FALSE)
  verbose_print(paste0("Created directory for sequence table: ", sequence_table_dir), verbose) # nolint

  # build the sequence table
  seqtab <- makeSequenceTable(mergers)

  # filter out sequences that have length outside of a reasonable range around the expected length
  seqtab.filtered <- seqtab[, nchar(colnames(seqtab)) %in% seqtab_cfg$length_range]

  # save the sequence table
  saveRDS(seqtab, file.path(sequence_table_dir, paste0(flowcell, "_seqtab.rds")))
  saveRDS(seqtab.filtered, file.path(sequence_table_dir, paste0(flowcell, "_seqtab_filtered.rds")))

  end_time <- Sys.time()
  verbose_print(paste0("Sequence table build complete. Time elapsed: ", hms_elapsed(start_time, end_time)), verbose) # nolint
}

if (!is.null(opt$flowcell)) {
  build_sequence_table(opt$flowcell, verbose = verbose)
} else {
  verbose_print("No flow cell specified, building sequence table for all flow cells...", verbose) # nolint
  flowcells <- list.dirs(output_directory, recursive = FALSE, full.names = FALSE)
  verbose_print("Flowcell directories found in output_directory:", verbose) # nolint
  cat(paste0("  ", flowcells, collapse = "\n"), "\n")
  for (flowcell in flowcells) {
    build_sequence_table(flowcell, verbose = verbose)
  }
}