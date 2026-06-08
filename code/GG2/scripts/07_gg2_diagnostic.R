# 07_gg2_diagnostic.R
# --------------------
# Generalized GG2 diagnostic script operating on the cleaned outputs:
#   - PI/final/seqtab.csv
#   - PI/final/taxonomy.csv
#   - PI/intermediate/exports/tree.nwk
#   - PI/intermediate/exports/asv_to_gg2_map.tsv (for occupancy = ASVs per GG2 feature)
#   - input/table_biom_ready.tsv (for "original" table)
#
# Produces:
#   - gg2_occupancy_distribution.jpg
#   - gg2_occupancy_data.tsv
#   - read_distribution_counts.jpg
#   - read_distribution_fractions.jpg
#   - read_distribution_scatter.jpg
#   - gg2_tree_phyla.jpg

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
  library(ape)
  library(ggplot2)
  library(ggtree)
  library(dplyr)
  library(readr)
  library(optparse)
  library(tidyr)
})

source("./code/utility_functions.R")
source("./code/GG2/functions/diagnostic_functions.R")

option_list <- list(
  make_option(c("--config"), type = "character", default = NULL,
              help = "Path to GG2 config R file (evaluates to a list).",
              metavar = "FILE"),
  make_option(c("--output_path", "-o"), type = "character", default = NULL,
              help = "Optional: base output directory (default: from config output$directory). Diagnostics go to <base>/<perc-identity>/diagnostics/.",
              metavar = "DIR"),
  make_option(c("--perc-identity", "-p"), type = "character", default = "0.99",
              help = "Percentage identity threshold used in GG2 mapping (default: 0.99).",
              metavar = "CHAR"),
  make_option(c("--verbose", "-v"), action = "store_true", default = FALSE,
              help = "Print detailed progress information.")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$config) || !nzchar(opt$config)) {
  print_help(opt_parser)
  stop("--config is required.", call. = FALSE)
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

data_base <- if (!is.null(opt$output_path) && nzchar(opt$output_path)) opt$output_path else cfg$output$directory
input_dir <- file.path(data_base, "input")
pi_dir    <- file.path(data_base, pi)
inter_dir <- file.path(pi_dir, "intermediate")
export_dir <- file.path(inter_dir, "exports")
final_dir  <- file.path(pi_dir, "final")
diag_dir   <- file.path(pi_dir, "diagnostics")

if (!dir.exists(diag_dir)) {
  dir.create(diag_dir, recursive = TRUE)
}

verbose_print(paste("Starting GG2 diagnostics for perc-identity:", pi), opt$verbose)

############################
## 1. Read inputs        ##
############################

orig_table_path <- file.path(input_dir, "table_biom_ready.tsv")
map_csv_path   <- file.path(input_dir, "sample_id_mapping.csv")
seqtab_path     <- file.path(final_dir, "seqtab.csv")
tax_csv_path    <- file.path(final_dir, "taxonomy.csv")
tree_path       <- file.path(export_dir, "tree.nwk")
asv_to_gg2_path <- file.path(export_dir, "asv_to_gg2_map.tsv")

required_files <- c(orig_table_path, map_csv_path, seqtab_path, tax_csv_path, tree_path, asv_to_gg2_path)
for (fp in required_files) {
  if (!file.exists(fp)) {
    stop("Required file not found for diagnostics: ", fp)
  }
}
verbose_print("All required diagnostic input files found.", opt$verbose)

orig_table <- read_tsv(orig_table_path, col_types = cols(.default = "c")) %>%
  mutate(across(starts_with("sample_"), as.numeric))

seqtab <- read_csv(seqtab_path, show_col_types = FALSE)

# Convert seqtab (sample_id × features) back to Feature × sample_ format
feature_cols <- setdiff(colnames(seqtab), "sample_id")

gg2_long <- seqtab %>%
  pivot_longer(cols = all_of(feature_cols),
               names_to = "Feature_ID",
               values_to = "reads") %>%
  group_by(Feature_ID, sample_id) %>%
  summarise(reads = sum(reads, na.rm = TRUE), .groups = "drop")

# Use sample_id -> sample_N mapping from input (matches 06_clean_outputs)
sample_map <- read_csv(map_csv_path, show_col_types = FALSE)
if (!"sample_N" %in% colnames(sample_map)) {
  stop("sample_id_mapping.csv must contain 'sample_N' column.")
}
if ("env_sample_id" %in% colnames(sample_map)) {
  sample_map <- sample_map %>%
    mutate(sample_id = as.character(env_sample_id))
} else {
  sample_map <- sample_map %>%
    mutate(sample_id = as.character(sample_N))
}
sample_map <- sample_map %>%
  mutate(sample_N = as.character(sample_N)) %>%
  select(sample_id, sample_N)

# Coerce sample_id to character so join matches mapping (seqtab may have numeric sample_id)
gg2_long <- gg2_long %>%
  mutate(sample_id = as.character(sample_id)) %>%
  left_join(sample_map, by = "sample_id")

if (any(is.na(gg2_long$sample_N))) {
  unmapped <- gg2_long %>% filter(is.na(sample_N)) %>% distinct(sample_id) %>% pull(sample_id)
  stop("Some sample_id in seqtab not found in sample_id_mapping.csv: ", paste(head(unmapped, 10), collapse = ", "),
       if (length(unmapped) > 10) " ..." else "")
}

gg2_table <- gg2_long %>%
  mutate(sample_col = sample_N) %>%
  select(`Feature ID` = Feature_ID, sample_col, reads) %>%
  pivot_wider(names_from = sample_col, values_from = reads, values_fill = 0) %>%
  arrange(`Feature ID`)

gg2_tax <- read_clean_taxonomy(tax_csv_path)
gg2_tree <- read.tree(tree_path)

############################
## 2. Occupancy          ##
############################
# Occupancy = number of original ASVs/OTUs grouped into each GG2 feature (from asv_to_gg2_map).

verbose_print("Step 2: Calculating occupancy distribution...", opt$verbose)

asv_to_gg2_map <- read_tsv(asv_to_gg2_path, col_types = cols(.default = "c"))
if (!all(c("InputFeatureID", "GG2FeatureID") %in% colnames(asv_to_gg2_map))) {
  stop("asv_to_gg2_map.tsv must contain columns InputFeatureID and GG2FeatureID.")
}

occupancy_data <- asv_to_gg2_map %>%
  group_by(Feature_ID = GG2FeatureID) %>%
  summarise(occupancy = n(), .groups = "drop")

write_tsv(occupancy_data, file.path(diag_dir, "gg2_occupancy_data.tsv"))

if (nrow(occupancy_data) == 0) {
  stop("No occupancy data found; check GG2 mapping.")
}

p_occupancy <- create_occupancy_distribution_plot(occupancy_data, pi)
ggsave(file.path(diag_dir, "gg2_occupancy_distribution.jpg"),
       plot = p_occupancy, width = 10, height = 6, dpi = 300)

############################
## 3. Read distribution  ##
############################

verbose_print("Step 3: Calculating read distributions...", opt$verbose)

read_distribution <- calculate_read_distribution(orig_table, gg2_table)
write_tsv(read_distribution, file.path(diag_dir, "read_distribution.tsv"))

p_reads_count    <- create_read_distribution_plot(read_distribution, use_fraction = FALSE)
p_reads_fraction <- create_read_distribution_plot(read_distribution, use_fraction = TRUE)
p_reads_scatter  <- create_read_distribution_scatter(orig_table, gg2_table)

ggsave(file.path(diag_dir, "read_distribution_counts.jpg"),
       plot = p_reads_count, width = 12, height = 8, dpi = 300)
ggsave(file.path(diag_dir, "read_distribution_fractions.jpg"),
       plot = p_reads_fraction, width = 12, height = 8, dpi = 300)
ggsave(file.path(diag_dir, "read_distribution_scatter.jpg"),
       plot = p_reads_scatter, width = 12, height = 8, dpi = 300)

############################
## 4. Tree plot          ##
############################

verbose_print("Step 4: Creating phylogenetic tree plot...", opt$verbose)

p_tree <- create_tree_plot(gg2_tree, gg2_tax)
ggsave(file.path(diag_dir, "gg2_tree_phyla.jpg"),
       plot = p_tree, width = 12, height = 12, dpi = 300)

cat("\n=== DIAGNOSTIC SUMMARY ===\n")
cat("✓ Occupancy distribution plot created\n")
cat("✓ Read distribution plots created (counts, fractions, scatter)\n")
cat("✓ Circular tree plot created\n")
cat("✓ All plots saved to:", diag_dir, "\n")
cat("\nDiagnostic completed successfully!\n")

