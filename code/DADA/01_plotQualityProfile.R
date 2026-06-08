# Script to do the first step of DADA2 pipeline: plot the q-score profiles

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
              help = "Plot the q-score profiles for a given flowcell (Default: NULL, all flow cells)."),
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
verbose <- if (!is.null(cfg$verbose)) cfg$verbose else opt$verbose

if (is.null(input_directory) || is.null(output_directory) ||
    is.null(forward_pattern) || is.null(reverse_pattern)) {
    stop("Config missing required fields: input$directory, output$directory, patterns$forward, patterns$reverse", call. = FALSE)
}

plot_q_score_profiles <- function(flowcell, input_directory, output_directory, forward_pattern, reverse_pattern, verbose = TRUE) {
    flow_cell_dir <- file.path(input_directory, flowcell)
    # get the forward and reverse files
    forward_files <- list.files(file.path(flow_cell_dir), pattern = forward_pattern, full.names = TRUE)
    reverse_files <- list.files(file.path(flow_cell_dir), pattern = reverse_pattern, full.names = TRUE)
    verbose_print(paste0('Found ', length(forward_files), ' forward files and ', length(reverse_files), ' reverse files in ', flow_cell_dir), verbose)
    # get the sample names by removing the pattern portion
    forward_sample_names <- sub(forward_pattern, "", basename(forward_files))
    reverse_sample_names <- sub(reverse_pattern, "", basename(reverse_files))

    # create directory for the forward and backward q-score profiles
    forward_output_dir <- file.path(output_directory, flowcell, '01_qualityProfile', 'forward')
    reverse_output_dir <- file.path(output_directory, flowcell, '01_qualityProfile', 'reverse')
    dir.create(forward_output_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(reverse_output_dir, recursive = TRUE, showWarnings = FALSE)

    # plot the q-score profiles
    for (i in seq_along(forward_files)) {
        p <- plotQualityProfile(forward_files[i])
        ggplot2::ggsave(file.path(forward_output_dir, paste0(forward_sample_names[i], '_forward.jpeg')),
                plot = p, width = 7.5, height = 7.5, dpi = 300)
    }
    for (i in seq_along(reverse_files)) {
        p <- plotQualityProfile(reverse_files[i])
        ggplot2::ggsave(file.path(reverse_output_dir, paste0(reverse_sample_names[i], '_reverse.jpeg')),
                plot = p, width = 7.5, height = 7.5, dpi = 300)
    }
    verbose_print(paste0('Q-score profiles for flowcell ', flowcell, ' saved to ', forward_output_dir, ' and ', reverse_output_dir, '\n'), verbose)
}

if (!is.null(opt$flowcell)) {
    plot_q_score_profiles(opt$flowcell, input_directory, output_directory, forward_pattern, reverse_pattern, verbose = verbose)
} else {
    verbose_print('No flow cell specified, filtering all flow cells...', verbose = verbose)
    flowcells <- list.dirs(input_directory, recursive = FALSE, full.names = FALSE)
    verbose_print('Flowcell directories found in input_directory:', verbose = verbose)
    cat(paste0("  ", flowcells, collapse = "\n"), "\n")
    for (flowcell in flowcells) {
        plot_q_score_profiles(flowcell, input_directory, output_directory, forward_pattern, reverse_pattern, verbose = verbose)
    }
}

# warnings()