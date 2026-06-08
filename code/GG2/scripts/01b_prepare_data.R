
# 01b_prepare_data.R
# ------------------
# Read a GG2 config, standardize the input feature table + sequences into
# a BIOM-ready TSV and filtered FASTA, and create a sample_id_mapping.csv.
#
# This script:
#   - Loads a config (list) from the given path.
#   - Calls `reformat_input()` with the appropriate format ("mitag" or "dada2").
#   - Writes:
#       <output$directory>/input/table_biom_ready.tsv
#       <output$directory>/input/sequences_filtered.fna
#       <output$directory>/input/sample_id_mapping.csv
#
# Environment activation is handled by the calling sbatch script; this
# script assumes the correct R environment is already active.

# Source repo setup (sets REPO_ROOT and working directory)
# Resolve this script's directory (works when run via Rscript)
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) != 1) stop("This script must be run with Rscript so --file= is available.")
script_path <- normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)
script_dir  <- dirname(script_path)
setup_path <- normalizePath(file.path(script_dir, "..", "..", "setup.R"),
                            winslash = "/", mustWork = TRUE)
source(setup_path)

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
})

source("./code/utility_functions.R")
source("./code/GG2/functions/reformat_input.R")

option_list <- list(
  make_option(c("--config"), type = "character", default = NULL,
              help = "Path to an R config file that evaluates to a list.",
              metavar = "FILE"),
  make_option(c("--verbose", "-v"), action = "store_true", default = FALSE,
              help = "Print detailed progress information.")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$config) || !nzchar(opt$config)) {
  stop("Missing required --config <path>", call. = FALSE)
}
if (!file.exists(opt$config)) {
  stop("Config file not found: ", opt$config, call. = FALSE)
}

cfg <- source(opt$config, local = TRUE)$value
if (!is.list(cfg)) {
  stop("Config must evaluate to an R list.", call. = FALSE)
}

dataset <- cfg$dataset
fmt     <- cfg$input$format
in_tab  <- cfg$input$table
in_seqs <- cfg$input$sequences
out_dir <- cfg$output$directory

if (is.null(fmt) || !fmt %in% c("mitag", "dada2")) {
  stop("Config$input$format must be 'mitag' or 'dada2'.")
}

verbose_print(paste0("[01b_prepare_data] Dataset: ", dataset,
                     " | format: ", fmt), opt$verbose)

input_dir <- file.path(out_dir, "input")
if (!dir.exists(input_dir)) {
  dir.create(input_dir, recursive = TRUE)
}

out_table  <- file.path(input_dir, "table_biom_ready.tsv")
out_seqs   <- file.path(input_dir, "sequences_filtered.fna")
map_path   <- file.path(input_dir, "sample_id_mapping.csv")

sample_mapping_cfg <- NULL
if (fmt == "mitag") {
  sample_mapping_cfg <- cfg$sample_mapping
}

res <- reformat_input(
  format            = fmt,
  input_table       = in_tab,
  input_sequences   = in_seqs,
  output_table      = out_table,
  output_sequences  = out_seqs,
  sample_mapping_cfg = sample_mapping_cfg,
  verbose           = opt$verbose
)

sample_mapping <- res$sample_mapping

if (!is.null(sample_mapping)) {
  write_csv(sample_mapping, map_path)
  verbose_print(paste0("[01b_prepare_data] Wrote sample_id_mapping.csv with ",
                       nrow(sample_mapping), " rows to ", map_path),
                opt$verbose)
} else {
  verbose_print("[01b_prepare_data] No sample mapping returned by reformat_input().",
                opt$verbose)
}

verbose_print("[01b_prepare_data] Completed successfully.", opt$verbose)

