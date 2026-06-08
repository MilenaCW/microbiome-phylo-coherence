# Script to do diagnostics on the DADA2 pipeline

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
              help = "Run diagnostics for a given flowcell (Default: NULL, all flow cells)."),
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
patterns_cfg <- cfg$patterns
verbose <- if (!is.null(cfg$verbose)) cfg$verbose else opt$verbose

if (is.null(output_directory)) {
  stop("Config missing required field: output$directory", call. = FALSE)
}
if (is.null(seqtab_cfg) || is.null(seqtab_cfg$expected_length)) {
  stop("Config missing required field: sequence_table$expected_length", call. = FALSE)
}

run_diagnostics <- function(flowcell, verbose = TRUE) {
    start_time <- Sys.time()

    flow_cell_dir <- file.path(output_directory, flowcell)
    if (!dir.exists(flow_cell_dir)) {
      stop(paste0("Flowcell output directory not found: ", flow_cell_dir), call. = FALSE)
    }
    # make a directory for the diagnostics
    diagnostics_dir <- file.path(flow_cell_dir, "08a_dada_diagnostics")
    dir.create(diagnostics_dir, recursive = TRUE, showWarnings = FALSE)

    # read in the filter + trim file...
    filtered_dir <- file.path(flow_cell_dir, "02_filterAndTrim")
    if (!dir.exists(filtered_dir)) {
      stop(paste0("Filtered directory not found: ", filtered_dir), call. = FALSE)
    }
    if (is.null(patterns_cfg) ||
        is.null(patterns_cfg$forward) ||
        is.null(patterns_cfg$reverse)) {
      stop("Config missing required fields: patterns$forward, patterns$reverse", call. = FALSE)
    }
    ft_out <- read.csv(file.path(filtered_dir, paste0(flowcell, "_filter_trim.csv")), header = TRUE)
    ft_out <- ft_out %>%
        rename(sample_id = X) %>%  # Rename the empty first column to sample_id
        mutate(sample_id = as.character(sample_id)) %>%
        mutate(sample_id = basename(sample_id)) %>%
        mutate(sample_id = sub(patterns_cfg$reverse, "", sub(patterns_cfg$forward, "", sample_id)))
    # print("filter + trim out file (head):")
    # print(head(ft_out))

    # read in the sample inference files...
    sample_inference_dir <- file.path(flow_cell_dir, "04_sampleInference")
    if (!dir.exists(sample_inference_dir)) {
      stop(paste0("Sample inference directory not found: ", sample_inference_dir), call. = FALSE)
    }
    dadaFs <- readRDS(file.path(sample_inference_dir, paste0(flowcell, "_dadaFs.rds")))
    dadaRs <- readRDS(file.path(sample_inference_dir, paste0(flowcell, "_dadaRs.rds")))
    # print("dadaFs (head):")
    # print(head(dadaFs))
    # print("dadaRs (head):")
    # print(head(dadaRs))

    # read in the merged pairs file...
    merged_dir <- file.path(flow_cell_dir, "05_merge")
    if (!dir.exists(merged_dir)) {
      stop(paste0("Merged directory not found: ", merged_dir), call. = FALSE)
    }
    mergers <- readRDS(file.path(merged_dir, paste0(flowcell, "_mergers.rds")))
    # print("mergers (head):")
    # print(head(mergers))

    # read in the sequence table file...
    sequence_table_dir <- file.path(flow_cell_dir, "06_buildTable")
    if (!dir.exists(sequence_table_dir)) {
      stop(paste0("Sequence table directory not found: ", sequence_table_dir), call. = FALSE)
    }
    chimera_dir <- file.path(flow_cell_dir, "07_removeChimera")
    if (!dir.exists(chimera_dir)) {
      stop(paste0("Chimera directory not found: ", chimera_dir), call. = FALSE)
    }
    seqtab <- readRDS(file.path(sequence_table_dir, paste0(flowcell, "_seqtab.rds")))
    seqtab.filtered <- readRDS(file.path(sequence_table_dir, paste0(flowcell, "_seqtab_filtered.rds")))
    seqtab.nochim <- readRDS(file.path(chimera_dir, paste0(flowcell, "_seqtab_noChimeras.rds")))
    # print("seqtab (head):")
    # print(head(seqtab))
    # print("seqtab.filtered (head):")
    # print(head(seqtab.filtered))
    # print("seqtab.nochim (head):")
    # print(head(seqtab.nochim))

    # count reads through the pipeline
    getN <- function(x) sum(getUniques(x))
    track_reads <- ft_out %>%
        rename(raw = reads.in, filtered = reads.out) %>%
        left_join(data.frame(sample_id = names(dadaFs), inference_forward = sapply(dadaFs, getN)), by = "sample_id") %>%
        left_join(data.frame(sample_id = names(dadaRs), inference_reverse = sapply(dadaRs, getN)), by = "sample_id") %>%
        left_join(data.frame(sample_id = names(mergers), merged = sapply(mergers, getN)), by = "sample_id") %>%
        left_join(data.frame(sample_id = rownames(seqtab), seqtab = rowSums(seqtab)), by = "sample_id") %>%
        left_join(data.frame(sample_id = rownames(seqtab.filtered), seqtab_filtered = rowSums(seqtab.filtered)), by = "sample_id") %>%
        left_join(data.frame(sample_id = rownames(seqtab.nochim), seqtab_nochim = rowSums(seqtab.nochim)), by = "sample_id") %>%
        select(sample_id, raw, filtered, inference_forward, inference_reverse, merged, seqtab, seqtab_filtered, seqtab_nochim)
    print("track_reads (head):")
    print(head(track_reads))
    write.csv(track_reads, file.path(diagnostics_dir, paste0(flowcell, "_track_reads.csv")), row.names = FALSE)
    # plot the reads through the pipeline
    track_reads %>%
        pivot_longer(cols = -sample_id, names_to = "step", values_to = "reads") %>%
        mutate(step = factor(step, levels = c("raw", "filtered", "inference_forward", "inference_reverse", "merged", "seqtab", "seqtab_filtered", "seqtab_nochim"))) %>%
        ggplot(aes(x = sample_id, y = reads, color = step)) +
        geom_point() +
        labs(x = "Sample", y = "Reads", color = "Step") +
        scale_y_log10()
    ggsave(file.path(diagnostics_dir, paste0(flowcell, "_track_reads.jpeg")),
           width = 10, height = 6, dpi = 300)

    # plot the number of reads through the pipeline (as a fraction of the raw reads)
    track_reads %>%
        mutate(across(-sample_id, ~ . / raw)) %>%
        mutate(inference = pmin(inference_forward, inference_reverse)) %>%
        select(-c(inference_forward, inference_reverse)) %>%
        mutate(raw = raw - filtered,
               filtered = filtered - inference,
               inference = inference - merged,
               merged = merged - seqtab,
               seqtab = seqtab - seqtab_filtered,
               seqtab_filtered = seqtab_filtered - seqtab_nochim) %>%
        pivot_longer(cols = -sample_id, names_to = "step", values_to = "frac") %>%
        mutate(step = factor(step, levels = c("raw", "filtered", "inference", "merged", "seqtab", "seqtab_filtered", "seqtab_nochim"))) %>%
        ggplot(aes(x = factor(sample_id), y = frac, fill = step)) +
        geom_bar(stat = "identity", position = "stack") +
        geom_hline(yintercept = c(0.25, 0.5, 0.75), linetype = "dashed", color = "grey") +
        ylim(0, 1) +
        labs(x = "Sample", y = "Fraction of Raw Reads", fill = "Step",
             title = paste0(flowcell, " (", nrow(track_reads), " samples)"))
    ggsave(file.path(diagnostics_dir, paste0(flowcell, "_track_reads_fraction.jpeg")),
           width = 10, height = 6, dpi = 300)

    # plot the distribution of the sequence lengths for each stage of filtering
    seqtab_lengths <- data.frame(table(nchar(getSequences(seqtab)))) %>%
        rename(length = Var1, raw = Freq) %>%
        left_join(data.frame(table(nchar(getSequences(seqtab.filtered)))) %>% 
                    rename(length = Var1, filtered = Freq), by = "length") %>%
        left_join(data.frame(table(nchar(getSequences(seqtab.nochim)))) %>% 
                    rename(length = Var1, nochim = Freq), by = "length") %>% 
        pivot_longer(cols = -length, names_to = "step", values_to = "frequency") %>%
        mutate(length = as.numeric(as.character(length)), frequency = as.numeric(frequency)) %>%
        mutate(step = factor(step, levels = c("raw", "filtered", "nochim")))
    ggplot(seqtab_lengths, aes(x = length, y = frequency, color = step)) +
        geom_point() +
        geom_line(aes(group = step)) +
        geom_vline(xintercept = seqtab_cfg$expected_length, linetype = "dashed", color = "grey") +
        labs(x = "Sequence length", y = "Frequency", color = "Step") +
        scale_color_manual(
            labels = c("raw" = "Sequences after sample inference", "filtered" = "Sequences after bandpass filter", "nochim" = "Sequences after chimera removal"),
            values = c("raw" = "black", "filtered" = "red", "nochim" = "blue")
        )
    ggsave(file.path(diagnostics_dir, paste0(flowcell, "_sequence_lengths.jpeg")),
           width = 10, height = 6, dpi = 300)

    end_time <- Sys.time()
    verbose_print(paste0("Diagnostics complete. Time elapsed: ", hms_elapsed(start_time, end_time)), verbose)
}

if (!is.null(opt$flowcell)) {
    run_diagnostics(opt$flowcell, verbose = verbose)
} else {
    verbose_print("No flow cell specified, running diagnostics for all flow cells...", verbose)
    flowcells <- list.dirs(output_directory, recursive = FALSE, full.names = FALSE)
    verbose_print("Flowcell directories found in output_directory:", verbose)
    cat(paste0("  ", flowcells, collapse = "\n"), "\n")
    for (flowcell in flowcells) {
        run_diagnostics(flowcell, verbose = verbose)
    }
}