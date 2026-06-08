# 08_gg2_summary_diagnostic.R
# ---------------------------
# Generalized GG2 summary diagnostic: aggregates diagnostics across perc-identity
# values (e.g. after a sweep). Reads occupancy_data and read_distribution from
# each <data_base>/<perc-identity>/diagnostics/ (produced by 07_gg2_diagnostic.R)
# and produces:
#   1. Percent mapped/unmapped reads vs perc-identity
#   2. Mean occupancy (min–max error bars) vs perc-identity
#   3. Combined sweep summary table (gg2_sweep_summary.tsv)
#
# Usage: run after step 07 for multiple perc-identity values; uses --config
# (and optional --output_path) to locate the output directory.

# Source repo setup (sets REPO_ROOT and working directory)
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) != 1) stop("This script must be run with Rscript so --file= is available.")
script_path <- normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)
script_dir  <- dirname(script_path)
setup_path <- normalizePath(file.path(script_dir, "..", "..", "setup.R"),
                            winslash = "/", mustWork = TRUE)
source(setup_path)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(optparse)
  library(tidyr)
  library(scales)
})

option_list <- list(
  make_option(c("--config"), type = "character", default = NULL,
              help = "Path to GG2 config R file (evaluates to a list).",
              metavar = "FILE"),
  make_option(c("--output_path", "-o"), type = "character", default = NULL,
              help = "Optional: base output directory (default: from config output$directory).",
              metavar = "DIR"),
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

data_base <- if (!is.null(opt$output_path) && nzchar(opt$output_path)) opt$output_path else cfg$output$directory
summary_dir <- file.path(data_base, "summary_diagnostics")

verbose_print <- function(msg, verbose = FALSE) {
  if (verbose) cat(paste0("[INFO] ", msg, "\n"))
}

if (!dir.exists(data_base)) {
  stop("Output directory does not exist: ", data_base)
}
verbose_print(paste("Output base:", data_base), opt$verbose)

# -----------------------------------------------------------------------------
# 1. Find all perc-identity directories under data_base
# -----------------------------------------------------------------------------

verbose_print("Step 1: Finding perc-identity directories...", opt$verbose)

all_dirs <- list.dirs(data_base, recursive = FALSE, full.names = FALSE)
perc_identity_dirs <- all_dirs[grepl("^0\\.\\d{2}$", all_dirs)]

if (length(perc_identity_dirs) == 0) {
  stop("No perc-identity directories (0.XX) found in: ", data_base)
}

perc_identity_values <- sort(as.numeric(perc_identity_dirs))
perc_identity_dirs <- sprintf("%.2f", perc_identity_values)

verbose_print(paste("Found", length(perc_identity_dirs), "perc-identity dirs:",
                    paste(perc_identity_dirs, collapse = ", ")), opt$verbose)

# -----------------------------------------------------------------------------
# 2. Read occupancy data from each <data_base>/<pi>/diagnostics/
# -----------------------------------------------------------------------------

verbose_print("Step 2: Reading occupancy data from all directories...", opt$verbose)

occupancy_summary <- data.frame()

for (perc_id in perc_identity_dirs) {
  occupancy_file <- file.path(data_base, perc_id, "diagnostics", "gg2_occupancy_data.tsv")

  if (file.exists(occupancy_file)) {
    verbose_print(paste("Reading occupancy from:", perc_id), opt$verbose)
    occupancy_data <- read_tsv(occupancy_file, col_types = cols(.default = "c"))
    occupancy_data$occupancy <- as.numeric(occupancy_data$occupancy)

    summary_stats <- occupancy_data %>%
      summarise(
        perc_identity = as.numeric(perc_id),
        mean_occupancy = mean(occupancy, na.rm = TRUE),
        min_occupancy = min(occupancy, na.rm = TRUE),
        max_occupancy = max(occupancy, na.rm = TRUE),
        median_occupancy = median(occupancy, na.rm = TRUE),
        n_features = n(),
        .groups = "drop"
      )
    occupancy_summary <- rbind(occupancy_summary, summary_stats)
  } else {
    verbose_print(paste("Warning: occupancy file not found for", perc_id), opt$verbose)
  }
}

if (nrow(occupancy_summary) == 0) {
  stop("No occupancy data found in any perc-identity directory.")
}
verbose_print(paste("Read occupancy from", nrow(occupancy_summary), "directories"), opt$verbose)

# -----------------------------------------------------------------------------
# 3. Read read distribution data from each <data_base>/<pi>/diagnostics/
# -----------------------------------------------------------------------------

verbose_print("Step 3: Reading read distribution data from all directories...", opt$verbose)

read_distribution_summary <- data.frame()

for (perc_id in perc_identity_dirs) {
  read_dist_file <- file.path(data_base, perc_id, "diagnostics", "read_distribution.tsv")

  if (file.exists(read_dist_file)) {
    verbose_print(paste("Reading read distribution from:", perc_id), opt$verbose)
    read_dist_data <- read_tsv(read_dist_file, col_types = cols(.default = "c"))

    numeric_cols <- c("total_reads", "gg2_mapped_reads", "unmapped_reads",
                     "gg2_mapped_fraction", "unmapped_fraction")
    numeric_cols <- intersect(numeric_cols, colnames(read_dist_data))
    read_dist_data <- read_dist_data %>%
      mutate(across(all_of(numeric_cols), as.numeric))

    summary_stats <- read_dist_data %>%
      summarise(
        perc_identity = as.numeric(perc_id),
        total_reads = sum(total_reads, na.rm = TRUE),
        total_gg2_mapped = sum(gg2_mapped_reads, na.rm = TRUE),
        total_unmapped = sum(unmapped_reads, na.rm = TRUE),
        n_samples = n(),
        .groups = "drop"
      ) %>%
      mutate(
        total_gg2_mapped_fraction = total_gg2_mapped / total_reads,
        total_unmapped_fraction = total_unmapped / total_reads,
        gg2_mapped_percent = total_gg2_mapped_fraction * 100,
        unmapped_percent = total_unmapped_fraction * 100
      )

    read_distribution_summary <- rbind(read_distribution_summary, summary_stats)
  } else {
    verbose_print(paste("Warning: read distribution file not found for", perc_id), opt$verbose)
  }
}

if (nrow(read_distribution_summary) == 0) {
  stop("No read distribution data found in any perc-identity directory.")
}
verbose_print(paste("Read read distribution from", nrow(read_distribution_summary), "directories"), opt$verbose)

# -----------------------------------------------------------------------------
# 4. Create percent mapped/unmapped vs perc-identity plot
# -----------------------------------------------------------------------------

verbose_print("Step 4: Creating percent mapped/unmapped reads plot...", opt$verbose)

dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

plot_data <- read_distribution_summary %>%
  select(perc_identity, gg2_mapped_percent, unmapped_percent) %>%
  pivot_longer(c(gg2_mapped_percent, unmapped_percent),
               names_to = "category", values_to = "percentage") %>%
  mutate(category = ifelse(category == "gg2_mapped_percent", "GG2 Mapped", "Unmapped"))

p_reads <- ggplot(plot_data, aes(x = perc_identity * 100, y = percentage)) +
  geom_point(aes(shape = category), size = 3) +
  scale_shape_manual(values = c("GG2 Mapped" = 19, "Unmapped" = 21)) +
  ylim(0, 100) +
  xlim(80, 100) +
  labs(
    title = "Read Mapping Success",
    x = "Percentage Identity Match (%)",
    y = "Percentage of Reads (%)",
    shape = "Category"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(summary_dir, "gg2_mapping_success_vs_perc_identity.jpg"),
       plot = p_reads, width = 12, height = 8, dpi = 300)
verbose_print("Saved read mapping plot", opt$verbose)

# -----------------------------------------------------------------------------
# 5. Create mean occupancy with error bars vs perc-identity plot
# -----------------------------------------------------------------------------

verbose_print("Step 5: Creating mean occupancy with error bars plot...", opt$verbose)

p_occupancy <- ggplot(occupancy_summary, aes(x = 100 * perc_identity, y = mean_occupancy)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = min_occupancy, ymax = max_occupancy), width = 0, alpha = 0.7) +
  scale_y_log10() +
  xlim(80, 100) +
  labs(
    title = "ASV Occupancy per GG2 Feature",
    subtitle = "Point is mean; error bars show min to max occupancy range",
    x = "Percentage Identity Match (%)",
    y = "Occupancy (ASVs per GG2 Feature)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.title = element_text(size = 12),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(summary_dir, "gg2_occupancy_vs_perc_identity.jpg"),
       plot = p_occupancy, width = 12, height = 8, dpi = 300)
verbose_print("Saved occupancy plot", opt$verbose)

# -----------------------------------------------------------------------------
# 6. Create combined summary table
# -----------------------------------------------------------------------------

verbose_print("Step 6: Creating combined summary table...", opt$verbose)

combined_summary <- occupancy_summary %>%
  left_join(read_distribution_summary, by = "perc_identity") %>%
  select(perc_identity,
         mean_occupancy, min_occupancy, max_occupancy,
         gg2_mapped_percent, unmapped_percent,
         n_features, n_samples)

write_tsv(combined_summary, file.path(summary_dir, "gg2_sweep_summary.tsv"))
verbose_print("Saved summary table", opt$verbose)

# -----------------------------------------------------------------------------
# 7. Print summary
# -----------------------------------------------------------------------------

cat("\n=== GG2 SWEEP SUMMARY ===\n")
cat("Perc-identity range:", min(combined_summary$perc_identity), "to", max(combined_summary$perc_identity), "\n")
cat("Number of thresholds analyzed:", nrow(combined_summary), "\n")
cat("Mean occupancy range:", round(min(combined_summary$mean_occupancy), 2), "to", round(max(combined_summary$mean_occupancy), 2), "\n")
cat("GG2 mapping success range:", round(min(combined_summary$gg2_mapped_percent), 1), "% to", round(max(combined_summary$gg2_mapped_percent), 1), "%\n")
cat("\nOutputs written to:", summary_dir, "\n")
cat("  - gg2_mapping_success_vs_perc_identity.jpg\n")
cat("  - gg2_occupancy_vs_perc_identity.jpg\n")
cat("  - gg2_sweep_summary.tsv\n")
cat("\nSummary diagnostic completed successfully.\n")
