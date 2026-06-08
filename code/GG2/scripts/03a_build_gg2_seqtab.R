# 03a_build_gg2_seqtab.R
# -----------------------
# Build GG2 sequence table + representative sequences from vsearch UC output
# and the standardized BIOM-ready table produced in step 01b.
#
# Inputs:
#   --uc-file      : UC file from vsearch (otu_to_gg2.uc)
#   --table        : BIOM-ready TSV (input features × sample_N; first col 'Feature ID')
#   --gg2-backbone : GG2 backbone FASTA
#   --output-dir   : Directory where intermediate/exports live
#
# Outputs (written into <output-dir>/intermediate/exports/):
#   feature-table.tsv      : GG2 feature table (rows = GG2 features, first col 'Feature ID')
#   dna-sequences.fasta    : GG2 representative sequences for mapped features
#   asv_to_gg2_map.tsv     : mapping between input Feature ID and GG2 Feature ID

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
  library(dplyr)
  library(readr)
  library(tidyr)
  library(Biostrings)
  library(optparse)
})

option_list <- list(
  make_option(c("--uc-file", "-u"), type = "character", help = "vsearch UC output file"),
  make_option(c("--table", "-t"), type = "character", help = "BIOM-ready input table (Feature ID × sample_N)"),
  make_option(c("--gg2-backbone", "-g"), type = "character", help = "GG2 backbone sequences FASTA"),
  make_option(c("--output-dir", "-o"), type = "character", help = "Base output directory (PI-specific directory)"),
  make_option(c("--verbose", "-v"), action = "store_true", default = FALSE,
              help = "Print detailed progress information")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$`uc-file`) || is.null(opt$table) || is.null(opt$`gg2-backbone`) || is.null(opt$`output-dir`)) {
  stop("Required arguments: --uc-file, --table, --gg2-backbone, --output-dir", call. = FALSE)
}

uc_file    <- opt$`uc-file`
in_table   <- opt$table
backbone   <- opt$`gg2-backbone`
pi_dir     <- opt$`output-dir`
export_dir <- file.path(pi_dir, "intermediate", "exports")

if (!file.exists(uc_file))   stop("UC file not found: ", uc_file)
if (!file.exists(in_table))  stop("Input table not found: ", in_table)
if (!file.exists(backbone))  stop("GG2 backbone FASTA not found: ", backbone)

if (!dir.exists(export_dir)) {
  dir.create(export_dir, recursive = TRUE)
}

#############################
## Helper functions       ##
#############################

parse_uc_file <- function(uc_path, verbose = FALSE) {
  if (verbose) {
    cat("[03a_build_gg2_seqtab] Parsing UC file...\n")
  }
  uc_lines <- readLines(uc_path)
  hit_lines <- uc_lines[grepl("^H", uc_lines)]

  mappings <- data.frame(
    InputFeatureID = character(),
    GG2FeatureID  = character(),
    stringsAsFactors = FALSE
  )

  for (line in hit_lines) {
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) >= 10) {
      input_id <- parts[9]   # query ID
      gg2_id   <- parts[10]  # target ID
      if (gg2_id != "*") { # only include successful hits
        mappings <- rbind(mappings,
                          data.frame(InputFeatureID = input_id,
                                     GG2FeatureID  = gg2_id,
                                     stringsAsFactors = FALSE))
      }
    }
  }

  # Remove duplicates
  mappings <- unique(mappings)
  if (verbose) {
    cat("[03a_build_gg2_seqtab] Found", nrow(mappings), "input → GG2 mappings\n")
  }
  mappings
}

aggregate_to_gg2 <- function(mapping_df, in_table, verbose = FALSE) {
  if (verbose) {
    cat("[03a_build_gg2_seqtab] Aggregating counts to GG2 features...\n")
  }

  if (!"Feature ID" %in% colnames(in_table)) {
    stop("Input table must have a 'Feature ID' column.")
  }

  sample_cols <- colnames(in_table)[startsWith(colnames(in_table), "sample_")]
  if (length(sample_cols) == 0) {
    stop("No sample_ columns found in input table.")
  }

  in_table2 <- in_table %>%
    dplyr::rename("InputFeatureID" = "Feature ID")

  # Join the input table with the GG2 feature mapping
  joined <- in_table2 %>%
    inner_join(mapping_df, by = "InputFeatureID")

  gg2_table <- joined %>%
    group_by(GG2FeatureID) %>%
    summarise(across(all_of(sample_cols), ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
    dplyr::rename("Feature ID" = GG2FeatureID)

  list(table = gg2_table, sample_cols = sample_cols)
}

create_representative_sequences <- function(gg2_table, backbone_path, verbose = FALSE) {
  if (verbose) {
    cat("[03a_build_gg2_seqtab] Creating representative sequences...\n")
  }
  gg2_seqs <- readDNAStringSet(backbone_path)
  reps <- gg2_seqs[names(gg2_seqs) %in% gg2_table$`Feature ID`]
  if (verbose) {
    cat("[03a_build_gg2_seqtab] Kept", length(reps), "sequences from backbone\n")
  }
  reps
}

#############################
## Main body              ##
#############################

if (opt$verbose) {
  cat("[03a_build_gg2_seqtab] Using:\n")
  cat("  UC file     :", uc_file, "\n")
  cat("  Input table :", in_table, "\n")
  cat("  GG2 backbone:", backbone, "\n")
  cat("  Output dir  :", export_dir, "\n")
}

mapping_df <- parse_uc_file(uc_file, verbose = opt$verbose)
mapping_path <- file.path(export_dir, "asv_to_gg2_map.tsv")
write_tsv(mapping_df, mapping_path)

table_in <- read_tsv(in_table, col_types = cols(.default = "c")) %>%
  mutate(across(starts_with("sample_"), as.numeric))

agg <- aggregate_to_gg2(mapping_df, table_in, verbose = opt$verbose)
gg2_table <- agg$table
sample_cols <- agg$sample_cols

feature_table_path <- file.path(export_dir, "feature-table.tsv")
write_tsv(gg2_table, feature_table_path)

reps <- create_representative_sequences(gg2_table, backbone, verbose = opt$verbose)
seqs_path <- file.path(export_dir, "dna-sequences.fasta")
writeXStringSet(reps, seqs_path)

cat("\n[03a_build_gg2_seqtab] Summary:\n")
cat("  GG2 features :", nrow(gg2_table), "\n")
cat("  Samples      :", length(sample_cols), "\n")
cat("  feature-table:", feature_table_path, "\n")
cat("  sequences    :", seqs_path, "\n")
cat("  mapping      :", mapping_path, "\n")

