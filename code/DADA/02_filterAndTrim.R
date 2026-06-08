# Script to do the second step of DADA2 pipeline: filter and trim the reads

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
              help = "Filter and trim the reads for a given flowcell (Default: NULL, all flow cells)."),
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

input_directory <- cfg$input$directory
output_directory <- cfg$output$directory
forward_pattern <- cfg$patterns$forward
reverse_pattern <- cfg$patterns$reverse
filter_cfg <- cfg$filtering
verbose <- if (!is.null(cfg$verbose)) cfg$verbose else opt$verbose

if (is.null(input_directory) || is.null(output_directory) ||
    is.null(forward_pattern) || is.null(reverse_pattern)) {
  stop("Config missing required fields: input$directory, output$directory, patterns$forward, patterns$reverse", call. = FALSE)
}
if (is.null(filter_cfg) ||
    is.null(filter_cfg$truncQ) ||
    is.null(filter_cfg$truncLen) ||
    is.null(filter_cfg$maxN) ||
    is.null(filter_cfg$maxEE)) {
  stop("Config missing required fields: filtering$truncQ, filtering$truncLen, filtering$maxN, filtering$maxEE", call. = FALSE)
}

filter_and_trim <- function(flowcell) {
    start_time <- Sys.time()
    
    flow_cell_dir <- file.path(input_directory, flowcell)
    # get the forward and reverse files
    forward_files <- list.files(flow_cell_dir, pattern = forward_pattern, full.names = TRUE)
    reverse_files <- list.files(flow_cell_dir, pattern = reverse_pattern, full.names = TRUE)
    verbose_print(paste0('Found ', length(forward_files), ' forward files and ', length(reverse_files), ' reverse files in ', flow_cell_dir), verbose)
    # get the sample names
    forward_sample_names <- sub(forward_pattern, "", basename(forward_files))
    reverse_sample_names <- sub(reverse_pattern, "", basename(reverse_files))
    # set the names of the forward and reverse files to the sample names
    names(forward_files) <- forward_sample_names
    names(reverse_files) <- reverse_sample_names

    # create directory for the filtered and trimmed files
    filtered_dir <- file.path(output_directory, flowcell, '02_filterAndTrim')
    dir.create(filtered_dir, recursive = TRUE, showWarnings = FALSE)
    verbose_print(paste0('Created directory for filtered and trimmed files: ', filtered_dir), verbose)

    # filter and trim the reads
    filtered_forward_files <- file.path(filtered_dir, paste0(forward_sample_names, '_F_filt.fastq.gz'))
    filtered_reverse_files <- file.path(filtered_dir, paste0(reverse_sample_names, '_R_filt.fastq.gz'))
    names(filtered_forward_files) <- forward_sample_names
    names(filtered_reverse_files) <- reverse_sample_names

    # trim to the length we said and then also filter with a max number of expected errors
    # (maxEE). maxEE is calculated based on the function of what quality score is (as
    # reported by illumina).
    filter_args <- list(
      forward_files, filtered_forward_files, rev = reverse_files,
      filt.rev = filtered_reverse_files,
      truncQ = filter_cfg$truncQ,
      truncLen = filter_cfg$truncLen,
      maxN = filter_cfg$maxN,
      maxEE = filter_cfg$maxEE,
      rm.phix = TRUE,
      compress = TRUE,
      multithread = TRUE
    )
    if (!is.null(filter_cfg$trimLeft)) {
      filter_args$trimLeft <- filter_cfg$trimLeft
    }

    out <- do.call(filterAndTrim, filter_args)
    write.csv(out, file.path(filtered_dir, paste0(flowcell, '_filter_trim.csv')), row.names = TRUE)
    verbose_print(paste("Filter and trim summary saved to:", filtered_dir), verbose)

    end_time <- Sys.time()
    verbose_print(paste0("Filter + trim complete. Time elapsed: ", hms_elapsed(start_time,end_time)), verbose)
}

if (!is.null(opt$flowcell)) {
    filter_and_trim(opt$flowcell)
} else {
    verbose_print('No flow cell specified, filtering all flow cells...', verbose = verbose)
    flowcells <- list.dirs(input_directory, recursive = FALSE, full.names = FALSE)
    verbose_print('Flowcell directories found in input_directory:', verbose = verbose)
    cat(paste0("  ", flowcells, collapse = "\n"), "\n")
    for (flowcell in flowcells) {
        filter_and_trim(flowcell)
    }
}