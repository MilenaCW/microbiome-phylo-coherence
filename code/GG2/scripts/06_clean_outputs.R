# 06_clean_outputs.R
# -------------------
# Create cleaned, CCA-ready seqtab.csv and taxonomy.csv from GG2 intermediate
# outputs:
#   - PI/intermediate/exports/feature-table.tsv
#   - PI/intermediate/exports/taxonomy.tsv
#   - input/sample_id_mapping.csv
#
# This script:
#   1) Converts feature-table.tsv (Feature ID × sample_N) into seqtab.csv
#      (sample_id × GG2FeatureID columns, with actual feature IDs as headers).
#   2) Parses taxonomy.tsv (QIIME2 export) into taxonomy.csv with:
#        Feature_ID, Domain, Phylum, Class, Order, Family, Genus, Species
#      and ensures Feature_ID matches seqtab feature columns.

# Source repo setup (sets REPO_ROOT and working directory)
# Resolve this script's directory (works when run via Rscript)
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) != 1) stop("This script must be run with Rscript so --file= is available.")
script_path <- normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)
script_dir  <- dirname(script_path)
# Path to setup.R relative to THIS script (adjust ../.. as needed)
setup_path <- normalizePath(file.path(script_dir, "..", "..", "setup.R"),
                            winslash = "/", mustWork = TRUE)

source(setup_path)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(optparse)
})

source("./code/utility_functions.R")

option_list <- list(
  make_option(c("--config"), type = "character", default = NULL,
              help = "Path to GG2 config R file (evaluates to a list).",
              metavar = "FILE"),
  make_option(c("--perc-identity", "-p"), type = "character", default = "0.99",
              help = "Percentage identity threshold used in GG2 mapping (default: 0.99).",
              metavar = "CHAR"),
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

normalize_perc_identity <- function(input) {
  x <- trimws(input)
  if (grepl("^\\.", x)) x <- paste0("0", x)
  if (!grepl("^0\\.", x)) x <- paste0("0.", x)
  if (grepl("^0\\.[0-9]$", x)) x <- paste0(x, "0")
  x
}

pi <- normalize_perc_identity(opt$`perc-identity`)
out_base <- cfg$output$directory

input_dir <- file.path(out_base, "input")
pi_dir    <- file.path(out_base, pi)
inter_dir <- file.path(pi_dir, "intermediate")
export_dir <- file.path(inter_dir, "exports")
final_dir  <- file.path(pi_dir, "final")

if (!dir.exists(final_dir)) {
  dir.create(final_dir, recursive = TRUE)
}

feat_tsv   <- file.path(export_dir, "feature-table.tsv")
map_csv    <- file.path(input_dir, "sample_id_mapping.csv")

tax_src    <- file.path(export_dir, "taxonomy.tsv")
tax_dest   <- file.path(final_dir, "taxonomy.csv")

seqs_src   <- file.path(export_dir, "dna-sequences.fasta")
seqs_dest  <- file.path(final_dir, "dna-sequences.fasta")

tree_src   <- file.path(export_dir, "tree.nwk")
tree_dest  <- file.path(final_dir, "tree.nwk")

for (p in c(feat_tsv, tax_src, map_csv)) {
  if (!file.exists(p)) {
    stop("Required file not found: ", p)
  }
}

verbose_print(paste("[06_clean_outputs] Using perc-identity:", pi), opt$verbose)

############################
## 1. Build seqtab.csv   ##
############################

feat_tbl <- read_tsv(feat_tsv, col_types = cols(.default = "c"))

if (!"Feature ID" %in% colnames(feat_tbl)) {
  stop("feature-table.tsv must have 'Feature ID' as first column.")
}

sample_cols <- colnames(feat_tbl)[startsWith(colnames(feat_tbl), "sample_")]
if (length(sample_cols) == 0) {
  stop("feature-table.tsv must contain sample_ columns.")
}

feat_tbl_num <- feat_tbl %>%
  mutate(across(all_of(sample_cols), as.numeric))

sample_map <- read_csv(map_csv, show_col_types = FALSE)
if (!"sample_N" %in% colnames(sample_map)) {
  stop("sample_id_mapping.csv must contain 'sample_N' column.")
}

# Use env_sample_id when present; otherwise fall back to sample_N as character
# Mapping file from 01b has sample_N as "sample_1", "sample_2", ... (see reformat_input.R)
if ("env_sample_id" %in% colnames(sample_map)) {
  sample_map <- sample_map %>%
    mutate(sample_id = as.character(env_sample_id))
} else {
  sample_map <- sample_map %>%
    mutate(sample_id = as.character(sample_N))
}

# Keep sample_N as-is (01b writes "sample_1", "sample_2", ...)
sample_map <- sample_map %>%
  mutate(sample_N = as.character(sample_N))

gg2_long <- feat_tbl_num %>%
  pivot_longer(
    cols      = all_of(sample_cols),
    names_to  = "sample_N",
    values_to = "reads"
  )

gg2_long <- gg2_long %>%
  left_join(sample_map, by = "sample_N")

if (any(is.na(gg2_long$sample_id))) {
  n_unmapped <- sum(is.na(gg2_long$sample_id))
  unmapped_sample_N <- gg2_long %>%
    filter(is.na(sample_id)) %>%
    distinct(sample_N) %>%
    pull(sample_N)
  stop(
    "Some sample_N values in feature-table.tsv could not be mapped to sample_id. ",
    "Unmapped rows: ", n_unmapped, ". ",
    "Unmapped sample_N: ", paste(unmapped_sample_N, collapse = ", ")
  )
}

seqtab <- gg2_long %>%
  select(sample_id, `Feature ID`, reads) %>%
  pivot_wider(
    names_from   = `Feature ID`,
    values_from  = reads,
    values_fill  = 0
  ) %>%
  arrange(sample_id)

seqtab_path <- file.path(final_dir, "seqtab.csv")
write_csv(seqtab, seqtab_path)
verbose_print(paste("[06_clean_outputs] Wrote seqtab.csv to", seqtab_path), opt$verbose)

############################
## 2. Build taxonomy.csv ##
############################

# Read and parse QIIME2 taxonomy.tsv (Feature ID, Taxon, [Confidence]) into rank columns.
read_gg2_taxonomy <- function(gg2_taxonomy_path) {
  read_tsv(gg2_taxonomy_path, col_types = cols(.default = "c")) %>%
    separate(Taxon, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
             sep = ";\\s*", remove = FALSE) %>%
    mutate(
      Domain = gsub("^d__", "", Domain),
      Phylum = gsub("^p__", "", Phylum),
      Class  = gsub("^c__", "", Class),
      Order  = gsub("^o__", "", Order),
      Family = gsub("^f__", "", Family),
      Genus  = gsub("^g__", "", Genus),
      Species = gsub("^s__", "", Species)
    ) %>%
    mutate(across(c(Domain, Phylum, Class, Order, Family, Genus, Species), clean_tax_name)) %>%
    select(-Taxon)
}

raw_tax <- read_tsv(tax_src, col_types = cols(.default = "c"))

if (!"Feature ID" %in% colnames(raw_tax)) {
  stop("taxonomy.tsv must have 'Feature ID' as first column.")
}

parsed_tax <- read_gg2_taxonomy(tax_src)

clean_tax <- parsed_tax %>%
  rename(Feature_ID = `Feature ID`)

feature_cols <- setdiff(colnames(seqtab), "sample_id")

clean_tax <- clean_tax %>%
  filter(Feature_ID %in% feature_cols)

missing_in_tax <- setdiff(feature_cols, clean_tax$Feature_ID)
if (length(missing_in_tax) > 0) {
  warning("[06_clean_outputs] Some features in seqtab.csv are missing from taxonomy.tsv: ",
          paste(head(missing_in_tax, 10), collapse = ", "),
          if (length(missing_in_tax) > 10) " ...")
}

clean_tax <- clean_tax %>%
  arrange(Feature_ID)

write_csv(clean_tax, tax_dest)
verbose_print(paste("[06_clean_outputs] Wrote taxonomy.csv to", tax_dest), opt$verbose)

############################
## 3. Copy sequences      ##
############################

if (file.exists(seqs_src)) {
  file.copy(seqs_src, seqs_dest, overwrite = TRUE)
  verbose_print(paste("[06_clean_outputs] Copied dna-sequences.fasta to", seqs_dest), opt$verbose)
} else {
  warning("[06_clean_outputs] dna-sequences.fasta not found at ", seqs_src)
}

if (file.exists(tree_src)) {
  file.copy(tree_src, tree_dest, overwrite = TRUE)
  verbose_print(paste("[06_clean_outputs] Copied tree.nwk to", tree_dest), opt$verbose)
} else {
  warning("[06_clean_outputs] tree.nwk not found at ", tree_src)
}

verbose_print("[06_clean_outputs] Completed successfully.", opt$verbose)