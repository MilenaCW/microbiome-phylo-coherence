# Script to do the fifth step of DADA2 pipeline: merge paired reads

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
              help = "Merge forward and reverse reads for a given flowcell (Default: NULL, all flow cells)."),
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

merge_pairs <- function(flowcell, verbose = TRUE) {
    start_time <- Sys.time()

    filtered_dir <- file.path(output_directory, flowcell, "02_filterAndTrim")
    if (!dir.exists(filtered_dir)) {
      stop(paste0("Filtered directory not found: ", filtered_dir), call. = FALSE)
    }

    sample_inference_dir <- file.path(output_directory, flowcell, "04_sampleInference")
    if (!dir.exists(sample_inference_dir)) {
      stop(paste0("Sample inference directory not found: ", sample_inference_dir), call. = FALSE)
    }

    # get the filtered files
    forward_files <- list.files(filtered_dir, pattern = "_F_filt.fastq.gz", full.names = TRUE)
    reverse_files <- list.files(filtered_dir, pattern = "_R_filt.fastq.gz", full.names = TRUE)
    verbose_print(paste0("Found ", length(forward_files), " forward files and ", length(reverse_files),
                         " reverse files in ", filtered_dir), verbose)
    if (length(forward_files) == 0 || length(reverse_files) == 0) {
      stop(paste0("No filtered files found in: ", filtered_dir), call. = FALSE)
    }
    names(forward_files) <- sub("_F_filt.fastq.gz", "", basename(forward_files))
    names(reverse_files) <- sub("_R_filt.fastq.gz", "", basename(reverse_files))

    # get the sample inference results
    dadaFs <- readRDS(file.path(sample_inference_dir, paste0(flowcell, "_dadaFs.rds")))
    dadaRs <- readRDS(file.path(sample_inference_dir, paste0(flowcell, "_dadaRs.rds")))

    # create directory for the merged pairs
    merged_dir <- file.path(output_directory, flowcell, "05_merge")
    dir.create(merged_dir, recursive = TRUE, showWarnings = FALSE)
    verbose_print(paste0("Created directory for merged pairs: ", merged_dir), verbose)

    # merge the pairs
    mergers <- mergePairs(dadaFs, forward_files, dadaRs, reverse_files, verbose = verbose)
    # save the merged pairs
    saveRDS(mergers, file.path(merged_dir, paste0(flowcell, "_mergers.rds")))

    end_time <- Sys.time()
    verbose_print(paste0("Forward + reverse reads merged. Time elapsed: ", hms_elapsed(start_time, end_time)), verbose)
}

if (!is.null(opt$flowcell)) {
    merge_pairs(opt$flowcell, verbose = verbose)
} else {
    verbose_print("No flow cell specified, merging forward and reverse reads for all flow cells...", verbose)
    flowcells <- list.dirs(output_directory, recursive = FALSE, full.names = FALSE)
    verbose_print("Flowcell directories found in output_directory:", verbose)
    cat(paste0("  ", flowcells, collapse = "\n"), "\n")
    for (flowcell in flowcells) {
        merge_pairs(flowcell, verbose = verbose)
    }
}