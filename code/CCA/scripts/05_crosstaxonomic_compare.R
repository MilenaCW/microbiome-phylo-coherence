#!/usr/bin/env Rscript
# 05_crosstaxonomic_compare.R — Compare CCA results across taxonomic levels; correlation and env loadings plots.
# Usage: Rscript 05_crosstaxonomic_compare.R --config <path> --perc_identity <p> [--verbose]

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) != 1) stop("Run with Rscript so --file= is available.")
script_path <- normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)
script_dir <- dirname(script_path)
setup_path <- normalizePath(file.path(script_dir, "..", "..", "setup.R"), winslash = "/", mustWork = TRUE)
source(setup_path)

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(readr)
  library(ggplot2)
})
source(file.path(script_dir, "..", "functions", "CCA_functions.R"))
source(file.path(script_dir, "..", "..", "utility_functions.R"))

# Hard-coded order of taxonomic levels (not dataset-specific)
tax_levels_order <- c("OTU", "Species", "Genus", "Family", "Order", "Class", "Phylum")

# Default tax palette if not in config
default_tax_palette <- c(
  "#332288", "#117733", "#44AA99", "#88CCEE",
  "#DDCC77", "#CC6677", "#AA4499", "#882255"
)

option_list <- list(
  make_option(c("--config"), type = "character", default = NULL, help = "Path to CCA config R file.", metavar = "FILE"),
  make_option(c("--perc_identity"), type = "character", default = "0.90", help = "Perc identity [default %default]."),
  make_option(c("--verbose", "-v"), action = "store_true", default = FALSE, help = "Verbose output.")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$config) || !nzchar(opt$config)) stop("Missing required --config", call. = FALSE)
if (is.null(opt$perc_identity) || !nzchar(opt$perc_identity)) stop("Missing required --perc_identity", call. = FALSE)
if (!file.exists(opt$config)) stop("Config file not found: ", opt$config, call. = FALSE)
cfg <- source(opt$config, local = TRUE)$value
if (!is.list(cfg)) stop("Config must evaluate to an R list.", call. = FALSE)

root <- if (exists("REPO_ROOT")) REPO_ROOT else getwd()
results_path <- cfg$results_path
if (!grepl("^/", results_path)) results_path <- file.path(root, results_path)
compare_base <- file.path(results_path, opt$perc_identity)
compare_out <- file.path(compare_base, "crosstax_compare")
dir.create(compare_out, recursive = TRUE, showWarnings = FALSE)

plot_cfg <- cfg$plotting_params
tax_palette <- if (!is.null(plot_cfg$tax_palette) && length(plot_cfg$tax_palette) > 0) {
  plot_cfg$tax_palette
} else {
  default_tax_palette
}

group_label <- paste0(as.numeric(opt$perc_identity) * 100, "% ID Matching")
verbose_print("Starting crosstaxonomic CCA compare...", verbose = opt$verbose)
verbose_print(paste("Parameters: perc_identity =", opt$perc_identity, ", output =", compare_out), verbose = opt$verbose)

# Discover tax levels: dirs under compare_base that have step2_loadings/correlations_per_fold.csv
all_dirs <- list.dirs(compare_base, full.names = FALSE, recursive = FALSE)
tax_dirs <- character(0)
for (d in all_dirs) {
  if (d == "" || d == "crosstax_compare" || d == "Domain") next
  corr_file <- file.path(compare_base, d, "step2_loadings", "correlations_per_fold.csv")
  if (file.exists(corr_file)) tax_dirs <- c(tax_dirs, d)
}

if (opt$verbose) {
  verbose_print(paste("Found", length(tax_dirs), "taxonomic level directories:", paste(tax_dirs, collapse = ", ")), verbose = TRUE)
}

if (length(tax_dirs) == 0) {
  stop("No tax-level directories with step2_loadings/correlations_per_fold.csv found under ", compare_base, call. = FALSE)
}

# Order factor levels: use tax_levels_order, then any discovered level not in the list
plot_tax_levels <- unique(c(tax_levels_order, tax_dirs))
plot_tax_levels <- plot_tax_levels[plot_tax_levels %in% tax_dirs]

corr_list <- list()
env_loadings_list <- list()
null_expectations_list <- list()

for (tax_level in tax_dirs) {
  verbose_print(paste("Processing", tax_level, "..."), verbose = opt$verbose)
  step2_dir <- file.path(compare_base, tax_level, "step2_loadings")
  step3_dir <- file.path(compare_base, tax_level, "step3_null")

  # Correlations: from correlations_per_fold (canonical_direction, fold, train, test)
  corr_file <- file.path(step2_dir, "correlations_per_fold.csv")
  corr_raw <- read_csv(corr_file, show_col_types = FALSE)
  corr_summary <- corr_raw %>%
    pivot_longer(cols = c(train, test), names_to = "type", values_to = "value") %>%
    group_by(canonical_direction, type) %>%
    summarise(mean = mean(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE), n = n(), .groups = "drop") %>%
    mutate(tax_level = tax_level) %>%
    select(tax_level, canonical_direction, type, mean, sd)
  corr_list[[tax_level]] <- corr_summary

  # Env loadings: from env_loadings.csv (canonical_direction, fold, var, value); already aligned in step 2
  env_file <- file.path(step2_dir, "env_loadings.csv")
  if (!file.exists(env_file)) {
    warning(paste("env_loadings.csv not found in", tax_level, "- skipping env loadings"))
  } else {
    env_raw <- read_csv(env_file, show_col_types = FALSE)
    env_summary <- env_raw %>%
      group_by(canonical_direction, var) %>%
      summarise(mean = mean(value, na.rm = TRUE), se = sd(value, na.rm = TRUE) / sqrt(n()), .groups = "drop") %>%
      mutate(tax_level = tax_level) %>%
      select(tax_level, canonical_direction, var, mean, se)
    env_loadings_list[[tax_level]] <- env_summary
  }

  # Null: from step3_null/null_correlations_per_fold.csv (canonical_direction, fold, train, test, seed)
  null_file <- file.path(step3_dir, "null_correlations_per_fold.csv")
  if (file.exists(null_file)) {
    null_raw <- read_csv(null_file, show_col_types = FALSE)
    null_cd1 <- null_raw %>%
      filter(canonical_direction == 1) %>%
      group_by(seed) %>%
      summarise(mean_test = mean(test, na.rm = TRUE), .groups = "drop") %>%
      summarise(mean = mean(mean_test), sd = sd(mean_test), n = n(), .groups = "drop") %>%
      mutate(tax_level = tax_level)
    null_expectations_list[[tax_level]] <- null_cd1
  }
}

all_correlations <- bind_rows(corr_list)
all_env_loadings <- bind_rows(env_loadings_list)
all_null_expectations <- bind_rows(null_expectations_list)

if (nrow(all_env_loadings) == 0) {
  stop("No env loadings found in any tax-level directory.", call. = FALSE)
}

# Align env loadings across tax levels (sign flip per canonical_direction)
aligned_env_loadings <- all_env_loadings
tax_levels <- unique(all_env_loadings$tax_level)
directions <- unique(all_env_loadings$canonical_direction)

verbose_print("Aligning environmental loadings across taxonomic levels...", verbose = opt$verbose)

for (CD in directions) {
  cd_data <- all_env_loadings %>% filter(canonical_direction == CD)
  ref_tax_level <- tax_levels[1]
  ref_data <- cd_data %>%
    filter(tax_level == ref_tax_level) %>%
    select(var, mean) %>%
    arrange(var)
  ref_values <- ref_data$mean

  for (tax in tax_levels) {
    if (tax == ref_tax_level) next
    tax_data <- cd_data %>%
      filter(tax_level == tax) %>%
      arrange(var)
    if (!setequal(ref_data$var, tax_data$var)) {
      stop(
        paste(
          "Variables do not match between", ref_tax_level, "and", tax, "for canonical_direction", CD
        ),
        call. = FALSE
      )
    }
    correlation_og <- cor(ref_data$mean, tax_data$mean)
    correlation_flipped <- cor(ref_data$mean, -tax_data$mean)
    if (correlation_og < correlation_flipped) {
      aligned_env_loadings <- aligned_env_loadings %>%
        mutate(mean = ifelse(canonical_direction == CD & tax_level == tax, -mean, mean))
      verbose_print(paste("  Canonical direction", CD, "- flipped", tax), verbose = opt$verbose)
    }
  }
}

verbose_print("Creating summary plots...", verbose = opt$verbose)

# ---- Plot 1: Correlation comparative ----
cor_plot_data <- all_correlations %>%
  filter(type == "test") %>%
  mutate(tax_level = factor(tax_level, levels = plot_tax_levels))

n_cd <- max(cor_plot_data$canonical_direction, 0)
ymin_val <- min(cor_plot_data$mean - cor_plot_data$sd, na.rm = TRUE)
ymin <- min(0, if (is.finite(ymin_val)) ymin_val else 0)
cor_plot <- ggplot(cor_plot_data, aes(x = factor(canonical_direction), y = mean, color = factor(canonical_direction))) +
  scale_color_manual(values = cov_palette, name = "Canonical direction") +
  labs(
    x = "Canonical direction",
    y = expression(rho[k] ~ "(mean ± sd)"),
    title = group_label
  ) +
  ylim(ymin, NA) +
  guides(color = "none") +
  facet_wrap(vars(tax_level)) +
  theme_minimal(base_size = 8)

if (nrow(all_null_expectations) > 0) {
  null_plot_data <- all_null_expectations %>%
    mutate(tax_level = factor(tax_level, levels = plot_tax_levels))
  cor_plot <- cor_plot +
    geom_rect(
      data = null_plot_data,
      aes(xmin = 0.5, xmax = n_cd + 0.5, ymin = mean - sd / sqrt(n), ymax = mean + sd / sqrt(n)),
      fill = "grey40", alpha = 0.15, inherit.aes = FALSE
    ) +
    geom_hline(
      data = null_plot_data,
      aes(yintercept = mean),
      linetype = "dashed", color = "grey40", linewidth = 0.7
    )
}
cor_plot <- cor_plot +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0, linewidth = 0.8)

ggsave(file.path(compare_out, "correlation_comparative.jpg"), plot = cor_plot, height = 4, width = 5, dpi = 150)
verbose_print(paste("Saved", file.path(compare_out, "correlation_comparative.jpg")), verbose = opt$verbose)

# ---- Plot 2: Env loadings comparative ----
# env_var_labels: named vector (names = var names in data; values = display labels), e.g. c(ph = "pH", soil_c = "C", ...)
env_var_labels <- plot_cfg$env_var_labels

env_plot_data <- aligned_env_loadings %>%
  mutate(tax_level = factor(tax_level, levels = plot_tax_levels))

if (!is.null(env_var_labels) && length(env_var_labels) > 0 && !is.null(names(env_var_labels))) {
  env_plot_data <- env_plot_data %>%
    mutate(var = factor(var, levels = names(env_var_labels), labels = unname(env_var_labels)))
} else {
  env_plot_data <- env_plot_data %>%
    mutate(var = factor(var, levels = sort(unique(var))))
}

# One color per tax level in plot order (recycle if config palette is shorter)
if (length(tax_palette) < length(plot_tax_levels)) {
  warning(
    "tax_palette has fewer colors (", length(tax_palette), ") than taxonomic levels (",
    length(plot_tax_levels), "); colors will be recycled."
  )
}
tax_palette_use <- rep(tax_palette, length.out = length(plot_tax_levels))
tax_palette_named <- setNames(tax_palette_use, plot_tax_levels)

# Set facet layout explicitly so we can match figure dimensions (height/width per facet)
n_facets <- dplyr::n_distinct(env_plot_data$canonical_direction)
n_facet_cols <- min(3L, max(2L, ceiling(sqrt(n_facets))))
n_facet_rows <- ceiling(n_facets / n_facet_cols)
# print(paste(n_facet_rows,"x", n_facet_cols))
inch_per_facet_height <- 1.5
inch_per_facet_width <- 1.5
env_plot_height <- 1 + n_facet_rows * inch_per_facet_height
env_plot_width <- 2 + n_facet_cols * inch_per_facet_width

env_plot <- env_plot_data %>%
  ggplot(aes(x = var, y = mean, color = tax_level)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey60") +
  geom_point(size = 1.5) +
  geom_line(aes(group = tax_level), linetype = "dotted", linewidth = 0.5) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0) +
  scale_color_manual(values = tax_palette_named, name = "taxonomic level") +
  facet_wrap(
    vars(canonical_direction),
    ncol = n_facet_cols,
    labeller = labeller(canonical_direction = function(x) paste0("canonical direction ", x))
  ) +
  labs(
    x = "feature",
    y = "environmental loading (mean ± se)",
    color = "taxonomic level"
  ) +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(file.path(compare_out, "env_loadings_comparative.jpg"), plot = env_plot, height = env_plot_height, width = env_plot_width, dpi = 150)
verbose_print(paste("Saved", file.path(compare_out, "env_loadings_comparative.jpg")), verbose = opt$verbose)

verbose_print("Crosstaxonomic compare done.", verbose = opt$verbose)
