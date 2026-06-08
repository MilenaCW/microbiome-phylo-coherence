# Script to do diagnostics on the DADA2 sequence table

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
verbose <- if (!is.null(cfg$verbose)) cfg$verbose else opt$verbose

if (is.null(output_directory)) {
  stop("Config missing required field: output$directory", call. = FALSE)
}

seqtab_diagnostics <- function(flowcell, verbose = TRUE) {
    start_time <- Sys.time()

    flow_cell_dir <- file.path(output_directory, flowcell)
    if (!dir.exists(flow_cell_dir)) {
      stop(paste0("Flowcell output directory not found: ", flow_cell_dir), call. = FALSE)
    }
    # make a directory for the diagnostics
    diagnostics_dir <- file.path(flow_cell_dir, "08b_seqtab_diagnostics")
    dir.create(diagnostics_dir, recursive = TRUE, showWarnings = FALSE)

    # read in the sequence table (ASV results from DADA output)
    chimera_dir <- file.path(flow_cell_dir, "07_removeChimera")
    if (!dir.exists(chimera_dir)) {
      stop(paste0("Chimera directory not found: ", chimera_dir), call. = FALSE)
    }
    seqtab <- readRDS(file.path(chimera_dir, paste0(flowcell, "_seqtab_noChimeras.rds")))
    seqtype <- "ASV"
    # print(rownames(seqtab))
    # print(colnames(seqtab))

    # move the sample id to a column
    seqtab <- seqtab %>%
        as.data.frame() %>%
        rownames_to_column(var = "sample_id") %>%
        mutate(sample_id = as.character(sample_id))

    # Define taxon numbering (with mapping to sequences) such that taxon_N is the Nth 
    # most prevalent taxon (number of samples with non-zero abundance)
    num_samples <- seqtab %>%
        pivot_longer(cols = -sample_id, names_to = "sequence", values_to = "reads") %>%
        group_by(sequence) %>%
        summarize(n_samples = sum(reads > 0)) %>%
        arrange(desc(n_samples)) %>%
        mutate(taxon_num = row_number())
    # Plot the number of samples with non-zero abundance for each taxon
    print("Plotting number of samples (with non-zero abundance) per taxon...")
    num_samples %>%
        mutate(rare = if_else(n_samples == 1, 'rare', 'shared')) %>%
        ggplot(aes(x = taxon_num, y = n_samples, color = rare)) +
        geom_point() +
        geom_hline(yintercept = c(1,nrow(seqtab)), linetype = "dashed", color = "grey") +
        scale_color_manual(values = c("shared" = "black", "rare" = "grey"), 
                        labels = c("shared" = "Shared (n_samples > 1)", "rare" = "Rare (n_samples == 1)")) +
        labs(x = seqtype, y = 'Number of samples with non-zero abundance',
            color = '', title = paste0("Flow cell ", flowcell, "\n", num_samples %>% filter(n_samples == 1) %>% nrow(), " rare ", seqtype, "s; ", num_samples %>% filter(n_samples > 1) %>% nrow(), " shared ", seqtype, "s")) +
        scale_x_log10()
    ggsave(file.path(diagnostics_dir, paste0(flowcell, "_num_samples_per_", seqtype, ".jpg")),width=6,height=5,dpi=300)

    # non-rare taxa
    no_rare_taxa <- num_samples %>%
        filter(n_samples > 1) %>%
        pull(taxon_num)
    seqtab <- seqtab %>%
        pivot_longer(cols = -sample_id, names_to = "sequence", values_to = "reads") %>%
        left_join(num_samples, by = "sequence") %>%
        select(sample_id, taxon_num, reads)
    seqtab_no_rare <- seqtab %>%
        filter(taxon_num %in% no_rare_taxa)
    
    # Plot the number of taxa per sample for the two conditions (full and no-rare)
    print("Plotting number of taxa per sample...")
    n_taxa_full <- seqtab %>% 
        pivot_longer(cols = !sample_id, names_to = 'taxon', values_to = 'reads') %>% 
        group_by(sample_id) %>% 
        summarize(n_taxa = sum(reads > 0)) %>%
        mutate(dataset = 'full')
    n_taxa_no_rare <- seqtab_no_rare %>% 
        pivot_longer(cols = !sample_id, names_to = 'taxon', values_to = 'reads') %>% 
        group_by(sample_id) %>% 
        summarize(n_taxa = sum(reads > 0)) %>%
        mutate(dataset = 'no_rare')
    n_taxa_all <- bind_rows(n_taxa_full, n_taxa_no_rare)
    # plot by number
    n_taxa_all %>% 
        ggplot(aes(x = sample_id, y = n_taxa, color = dataset)) +
        geom_point() +
        labs(x = 'Sample ID', y = paste0("Number of ", seqtype, "s per sample")) +
        ylim(0,NA) +
        scale_color_manual(values = c("full" = "black", "no_rare" = "grey"), labels = c("full" = "Full", "no_rare" = "No rare"))
    ggsave(file.path(diagnostics_dir, paste0(flowcell, "_n_", seqtype, "s_per_sample.jpg")),width=6,height=5,dpi=300)
    # plot by fraction
    n_taxa_all %>% 
        pivot_wider(names_from = dataset, values_from = n_taxa) %>%
        mutate(frac_full = full/full, frac_no_rare = no_rare/full) %>%
        select(-c(full,no_rare)) %>% 
        pivot_longer(cols = starts_with("frac_"), names_to = "dataset", values_to = "frac_taxa", names_prefix = "frac_") %>%
        ggplot(aes(x = sample_id, y = frac_taxa, color = dataset)) +
        geom_point() +
        labs(x = 'Sample ID', y = paste0("Fraction of ", seqtype, "s kept per sample")) +
        ylim(0,1) +
        scale_color_manual(values = c("full" = "black", "no_rare" = "grey"), labels = c("full" = "Full", "no_rare" = "No rare"))
    ggsave(file.path(diagnostics_dir, paste0(flowcell, "_frac_", seqtype, "s_per_sample.jpg")),width=6,height=5,dpi=300)

    # Plot the number of reads per sample for the two conditions (full and no-rare)
    print("Plotting number of reads per sample...")
    n_reads_full <- seqtab %>% 
        group_by(sample_id) %>% 
        summarize(n_reads = sum(reads)) %>%
        mutate(dataset = 'full')
    n_reads_no_rare <- seqtab_no_rare %>% 
        group_by(sample_id) %>% 
        summarize(n_reads = sum(reads)) %>%
        mutate(dataset = 'no_rare')
    n_reads_all <- bind_rows(n_reads_full, n_reads_no_rare)
    # plot by number
    n_reads_all %>% 
        ggplot(aes(x = sample_id, y = n_reads, color = dataset)) +
        geom_point() +
        labs(x = 'Sample ID', y = 'Number of reads per sample') +
        ylim(1,NA) +
        scale_y_log10() +
        scale_color_manual(values = c("full" = "black", "no_rare" = "grey"), labels = c("full" = "Full", "no_rare" = "No rare"))
    ggsave(file.path(diagnostics_dir, paste0(flowcell, "_", seqtype, "_n_reads_per_sample.jpg")),width=6,height=5,dpi=300)
    # plot by fraction
    n_reads_all %>% 
        pivot_wider(names_from = dataset, values_from = n_reads) %>%
        mutate(frac_full = full/full, frac_no_rare = no_rare/full) %>%
        select(-c(full,no_rare)) %>% 
        pivot_longer(cols = starts_with("frac_"), names_to = "dataset", values_to = "frac_reads", names_prefix = "frac_") %>%
        ggplot(aes(x = sample_id, y = frac_reads, color = dataset)) +
        geom_point() +
        labs(x = 'Sample ID', y = 'Fraction of reads kept per sample') +
        ylim(0,1) +
        scale_color_manual(values = c("full" = "black", "no_rare" = "grey"), labels = c("full" = "Full", "no_rare" = "No rare"))
    ggsave(file.path(diagnostics_dir, paste0(flowcell, "_", seqtype, "_frac_reads_per_sample.jpg")),width=6,height=5,dpi=300)

    # Plot the average relative abundance for the two conditions (full and no-rare)
    print("Plotting average relative abundance...")
    rel_abu_full <- seqtab %>% 
        group_by(sample_id) %>% 
        mutate(rel_abu = reads/sum(reads)) %>%
        ungroup() %>%
        group_by(taxon_num) %>%
        summarize(mean = mean(rel_abu), se = sd(rel_abu)/sqrt(n())) %>%
        mutate(dataset = 'full')
    rel_abu_no_rare <- seqtab_no_rare %>% 
        group_by(sample_id) %>% 
        mutate(rel_abu = reads/sum(reads)) %>%
        ungroup() %>%
        group_by(taxon_num) %>%
        summarize(mean = mean(rel_abu), se = sd(rel_abu)/sqrt(n())) %>%
        mutate(dataset = 'no_rare')
    rel_abu_all <- bind_rows(rel_abu_full, rel_abu_no_rare)
    # plot by number
    rel_abu_all %>% 
        ggplot(aes(x = taxon_num, y = mean, color = dataset)) +
        geom_point() +
        geom_errorbar(aes(ymin = mean - se, ymax = mean + se)) +
        labs(x = seqtype, y = 'Average relative abundance', color = '') +
        ylim(0,NA) +
        scale_x_log10() +
        scale_color_manual(values = c("full" = "black", "no_rare" = "grey"), labels = c("full" = "Full", "no_rare" = "No rare")) +
        facet_wrap(~dataset)
    ggsave(file.path(diagnostics_dir, paste0(flowcell, "_", seqtype, "_mean_relative_abundance.jpg")),width=6,height=5,dpi=300)

    end_time <- Sys.time()
    verbose_print(paste0("Seqtab diagnostics complete. Time elapsed: ", hms_elapsed(start_time, end_time)), verbose)
}

if (!is.null(opt$flowcell)) {
    seqtab_diagnostics(opt$flowcell, verbose = verbose)
} else {
    verbose_print("No flow cell specified, running diagnostics for all flow cells...", verbose)
    flowcells <- list.dirs(output_directory, recursive = FALSE, full.names = FALSE)
    verbose_print("Flowcell directories found in output_directory:", verbose)
    cat(paste0("  ", flowcells, collapse = "\n"), "\n")
    for (flowcell in flowcells) {
        seqtab_diagnostics(flowcell, verbose = verbose)
    }
}