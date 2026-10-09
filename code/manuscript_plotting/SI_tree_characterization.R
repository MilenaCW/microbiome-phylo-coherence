# SI_tree_characterization.R
# ---------------------------
# Produces a three-panel SI figure characterizing the pruned GG2 tree and the
# OTU abundance distribution (CCA input: rare OTUs removed) for soil and ocean datasets:
#   Panel A: OTUs per taxon (rank-size distribution) at each taxonomic rank
#   Panel B: cumulative relative abundance (mean across samples) of OTUs ranked
#            from most to least abundant
#   Panel C: nominal vs abundance-weighted effective group size, one point per
#            taxon (effective size = inverse Simpson index of within-group
#            abundance shares, i.e. 1 / sum_i p_i^2); dashed line is y = x.
#            Medians per rank are in the summary CSV.
#
# Usage (from repo root):
#   Rscript code/manuscript_plotting/SI_tree_characterization.R [--perc_identity 0.90]
#
# Output: manuscript/SI/SI_tree_characterization.pdf  (7 x 8 in)
#         manuscript/SI/SI_tree_characterization_summary.csv

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(patchwork)
  library(optparse)
})

# =============================================================================
# USER-ADJUSTABLE PARAMETERS
# =============================================================================

ranks <- c("Phylum", "Class", "Order", "Family", "Genus")

# =============================================================================
# PATHS
# =============================================================================

opt <- parse_args(OptionParser(option_list = list(
  make_option("--perc_identity", type = "character", default = "0.90")
)))

get_repo_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) != 1) {
    stop("This script must be run with Rscript so --file= is available.", call. = FALSE)
  }
  script_path <- normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)
  script_dir <- dirname(script_path)
  setup_path <- normalizePath(file.path(script_dir, "..", "setup.R"), winslash = "/", mustWork = TRUE)
  source(setup_path, local = TRUE)
  get("REPO_ROOT", envir = .GlobalEnv)
}

root <- get_repo_root()
out_pdf <- file.path(root, "manuscript/SI/SI_tree_characterization.pdf")
out_csv <- file.path(root, "manuscript/SI/SI_tree_characterization_summary.csv")

# =============================================================================
# COLORS  (consistent with all other manuscript plots)
# =============================================================================

dataset_colors <- c(
  "soil"  = "#8f723d",
  "ocean" = "#90afa7"
)

# =============================================================================
# DATA
# =============================================================================

# Use the same loader as the CCA (samples with complete env data, rare OTUs removed,
# per-sample relative abundance) so the figure describes the actual analysis input.
source(file.path(root, "code/CCA/functions/CCA_functions.R"))

load_dataset <- function(dataset) {
  data_path <- file.path(root, dataset, "data/processed_data")
  final <- file.path(data_path, "16S/GG2", opt$perc_identity, "final")
  env_file <- file.path(data_path, "environmental/filtered/envdata.csv")
  tax <- read_csv(file.path(final, "taxonomy.csv"), show_col_types = FALSE)
  data <- load_datasets_cca(env_file, file.path(final, "seqtab.csv"), read_threshold = 10)
  rel <- as.matrix(data$composition_df[, setdiff(names(data$composition_df), "sample_id")])
  otu_abund <- tibble(Feature_ID = colnames(rel), abund = colMeans(rel)) %>%
    filter(abund > 0)
  tax %>%
    inner_join(otu_abund, by = "Feature_ID") %>%
    mutate(dataset = dataset)
}

otu_df <- bind_rows(load_dataset("soil"), load_dataset("ocean")) %>%
  mutate(dataset = factor(dataset, levels = c("soil", "ocean")))

# Long format: one row per OTU per rank; OTUs without an assignment at a rank are dropped
long_df <- otu_df %>%
  pivot_longer(all_of(ranks), names_to = "rank", values_to = "taxon") %>%
  filter(!is.na(taxon), nzchar(taxon)) %>%
  mutate(rank = factor(rank, levels = ranks))

group_df <- long_df %>%
  group_by(dataset, rank, taxon) %>%
  summarise(
    n_otus = n(),
    eff_size = sum(abund)^2 / sum(abund^2),
    top_share = max(abund) / sum(abund),
    .groups = "drop"
  )

summary_df <- group_df %>%
  group_by(dataset, rank) %>%
  summarise(
    n_taxa = n(),
    median_nominal = median(n_otus),
    median_effective = median(eff_size),
    median_top_share = median(top_share),
    frac_otus_in_top10pct_taxa = {
      s <- sort(n_otus, decreasing = TRUE)
      sum(s[seq_len(ceiling(0.1 * length(s)))]) / sum(s)
    },
    .groups = "drop"
  )
print(as.data.frame(summary_df), digits = 3)
write_csv(summary_df, out_csv)

# Number of OTUs needed to reach given fractions of total abundance
cum_summary <- otu_df %>%
  group_by(dataset) %>%
  arrange(desc(abund), .by_group = TRUE) %>%
  summarise(
    n_otus_total = n(),
    n_for_50pct = which(cumsum(abund) / sum(abund) >= 0.5)[1],
    n_for_90pct = which(cumsum(abund) / sum(abund) >= 0.9)[1],
    n_for_99pct = which(cumsum(abund) / sum(abund) >= 0.99)[1],
    .groups = "drop"
  )
print(as.data.frame(cum_summary))
write_csv(cum_summary, sub("\\.csv$", "_cumulative.csv", out_csv))

# =============================================================================
# THEME
# =============================================================================

theme_si <- theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor    = element_blank(),
    plot.title          = element_text(hjust = 0.5, size = 10),
    plot.title.position = "panel",
    axis.text           = element_text(size = 10),
    axis.title          = element_text(size = 10),
    strip.text          = element_text(size = 10),
    legend.text         = element_text(size = 10),
    plot.margin         = margin(2, 4, 12, 2)
  )

# Log10 axes labelled as powers of ten (10^0, 10^1, ...), breaks at integer exponents only
pow10_breaks <- function(lim) 10^seq(ceiling(log10(lim[1])), floor(log10(lim[2])))
scale_x_pow10 <- function(...) scale_x_log10(breaks = pow10_breaks, labels = scales::label_math(10^.x, format = log10), ...)
scale_y_pow10 <- function(...) scale_y_log10(breaks = pow10_breaks, labels = scales::label_math(10^.x, format = log10), ...)

# =============================================================================
# PANEL A: OTUs per taxon, rank-size
# =============================================================================

rank_size_df <- group_df %>%
  group_by(dataset, rank) %>%
  arrange(desc(n_otus), .by_group = TRUE) %>%
  mutate(taxon_rank = row_number()) %>%
  ungroup()

p_a <- ggplot(rank_size_df, aes(x = taxon_rank, y = n_otus, color = dataset)) +
  geom_point(size = 0.6, alpha = 0.7) +
  facet_wrap(~rank, nrow = 1, scales = "free_x") +
  scale_color_manual(values = dataset_colors) +
  scale_x_pow10() +
  scale_y_pow10() +
  labs(x = "taxon rank (by size)", y = "OTUs per taxon", color = NULL, title = "(a)") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme_si

# =============================================================================
# PANEL B: OTU rank-abundance
# =============================================================================

abund_rank_df <- otu_df %>%
  group_by(dataset) %>%
  arrange(desc(abund), .by_group = TRUE) %>%
  mutate(otu_rank = row_number(), cum_abund = cumsum(abund) / sum(abund)) %>%
  ungroup()

p_b <- ggplot(abund_rank_df, aes(x = otu_rank, y = cum_abund, color = dataset)) +
  geom_point(size = 0.6, alpha = 0.7) +
  scale_color_manual(values = dataset_colors) +
  scale_x_pow10() +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "OTU rank (by abundance)", y = "cumulative relative abundance", color = NULL, title = "(b)") +
  guides(color = "none") +
  theme_si

# =============================================================================
# PANEL C: nominal vs effective group size, one point per taxon
# =============================================================================

p_c <- ggplot(group_df, aes(x = n_otus, y = eff_size, color = dataset)) +
  geom_abline(slope = 1, intercept = 0, color = "grey50", linewidth = 0.4, linetype = "dashed") +
  geom_point(size = 0.6, alpha = 0.3) +
  facet_wrap(~rank, nrow = 1, scales = "free_x") +
  scale_color_manual(values = dataset_colors) +
  scale_x_pow10() +
  scale_y_pow10() +
  labs(x = "nominal group size (OTUs)", y = "effective group size", color = NULL,
       title = "(c)") +
  guides(color = "none") +
  theme_si +
  theme(plot.margin = margin(2, 4, 2, 2))

# =============================================================================
# COMBINE + SAVE
# =============================================================================

fig <- (p_a / p_b / p_c) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0))
ggsave(out_pdf, fig, width = 7, height = 8)
message("Wrote ", out_pdf)
