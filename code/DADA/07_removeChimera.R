# Script to do the seventh step of DADA2 pipeline: remove chimeras

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
              help = "Remove chimeras for a given flowcell (Default: NULL, all flow cells)."),
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

remove_chimeras <- function(flowcell, verbose = TRUE) {
  start_time <- Sys.time()

  sequence_table_dir <- file.path(output_directory, flowcell, "06_buildTable")
  if (!dir.exists(sequence_table_dir)) {
    stop(paste0("Sequence table directory not found: ", sequence_table_dir), call. = FALSE)
  }

  # get the sequence table
  seqtab.filtered <- readRDS(file.path(sequence_table_dir, paste0(flowcell, "_seqtab_filtered.rds")))

  seqtab.nochim <- removeBimeraDenovo(seqtab.filtered, method = "consensus", multithread = TRUE, verbose = verbose)

  # save the sequence table
  chimera_dir <- file.path(output_directory, flowcell, "07_removeChimera")
  dir.create(chimera_dir, recursive = TRUE, showWarnings = FALSE)
  verbose_print(paste0("Created directory for chimera removal: ", chimera_dir), verbose) # nolint
  saveRDS(seqtab.nochim, file.path(chimera_dir, paste0(flowcell, "_seqtab_noChimeras.rds")))

  end_time <- Sys.time()
  verbose_print(paste0("Chimeras removed. Time elapsed: ", hms_elapsed(start_time, end_time)), verbose) # nolint
}

if (!is.null(opt$flowcell)) {
  remove_chimeras(opt$flowcell, verbose = verbose)
} else {
  verbose_print("No flow cell specified, removing chimeras for all flow cells...", verbose) # nolint
  flowcells <- list.dirs(output_directory, recursive = FALSE, full.names = FALSE)
  verbose_print("Flowcell directories found in output_directory:", verbose) # nolint
  cat(paste0("  ", flowcells, collapse = "\n"), "\n")
  for (flowcell in flowcells) {
    remove_chimeras(flowcell, verbose = verbose)
  }
}