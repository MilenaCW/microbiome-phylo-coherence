#!/usr/bin/env Rscript
# 04_model_performance_plots.R — Step 4: Test correlation vs null, env loadings scatter, abundance loadings on tree (one per CD, shared colormap, phylum ring).
# Usage: Rscript 04_model_performance_plots.R --config <path> --tax_level <level> --perc_identity <p> [--verbose] [--loadings_style continuous|binary]

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
  library(RColorBrewer)
  library(ggplot2)
  library(ggtree)
  library(ggnewscale)
  library(scales)
  library(ape)
})
source(file.path(script_dir, "..", "functions", "CCA_functions.R"))
source(file.path(script_dir, "..", "..", "utility_functions.R"))

option_list <- list(
  make_option(c("--config"), type = "character", default = NULL, help = "Path to CCA config R file.", metavar = "FILE"),
  make_option(c("--tax_level"), type = "character", default = "OTU", help = "Taxonomic level [default %default]"),
  make_option(c("--perc_identity"), type = "character", default = "0.90", help = "Perc identity [default %default]"),
  make_option(c("--verbose", "-v"), action = "store_true", default = FALSE, help = "Verbose output."),
  make_option(c("--loadings_style"), type = "character", default = "binary", help = "Tree loadings: continuous or binary [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$config) || !nzchar(opt$config)) stop("Missing required --config", call. = FALSE)
if (!file.exists(opt$config)) stop("Config file not found: ", opt$config, call. = FALSE)
cfg <- source(opt$config, local = TRUE)$value
if (!is.list(cfg)) stop("Config must evaluate to an R list.", call. = FALSE)

root <- if (exists("REPO_ROOT")) REPO_ROOT else getwd()
data_path <- cfg$data_path
if (!grepl("^/", data_path)) data_path <- file.path(root, data_path)
results_path <- cfg$results_path
if (!grepl("^/", results_path)) results_path <- file.path(root, results_path)
step2_dir <- file.path(results_path, opt$perc_identity, opt$tax_level, "step2_loadings")
step3_dir <- file.path(results_path, opt$perc_identity, opt$tax_level, "step3_null")
step4_dir <- file.path(results_path, opt$perc_identity, opt$tax_level, "step4_plots")
dir.create(step4_dir, recursive = TRUE, showWarnings = FALSE)

# ---- (1) Test correlation vs null ----
corr_per_fold <- read.csv(file.path(step2_dir, "correlations_per_fold.csv"))
true_stats <- corr_per_fold %>%
  group_by(canonical_direction) %>%
  summarise(mean = mean(test, na.rm = TRUE), sd = sd(test, na.rm = TRUE), .groups = "drop")
null_corr <- read.csv(file.path(step3_dir, "null_correlations_per_fold.csv"))
null_cd1 <- null_corr %>% filter(canonical_direction == 1)
null_per_seed <- null_cd1 %>% group_by(seed) %>% summarise(mean_test = mean(test, na.rm = TRUE), .groups = "drop")
null_mean <- mean(null_per_seed$mean_test)
null_se <- sd(null_per_seed$mean_test) / sqrt(nrow(null_per_seed))

y_raw_min <- min(true_stats$mean - true_stats$sd, na.rm = TRUE)
y_raw_max <- max(true_stats$mean + true_stats$sd, na.rm = TRUE)

p_corr <- ggplot(true_stats, aes(x = factor(canonical_direction), y = mean, color = factor(canonical_direction))) +
  geom_hline(yintercept = null_mean, linetype = "dashed", color = "grey40", linewidth = 0.7) +
  geom_rect(aes(xmin = 0.5, xmax = max(true_stats$canonical_direction) + 0.5, ymin = null_mean - null_se, ymax = null_mean + null_se),
            fill = "grey40", alpha = 0.15, inherit.aes = FALSE) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0, linewidth = 0.8) +
  scale_color_manual(values = cov_palette, name = "Canonical direction") +
  labs(x = "Canonical direction", y = "Test correlation (mean ± SD)", title = "True vs null (CD1)") +
  ylim(min(0,y_raw_min), max(1.0,y_raw_max)) +
  guides(color = "none")
ggsave(file.path(step4_dir, "test_correlation_vs_null.jpg"), plot = p_corr, width = 5, height = 4, dpi = 150)
verbose_print("Wrote test_correlation_vs_null.jpg", verbose = opt$verbose)

# ---- (2) Env loadings scatter ----
env_loadings <- read.csv(file.path(step2_dir, "env_loadings.csv"))
env_var_labels <- cfg$plotting_params$env_var_labels
if (!is.null(cfg$plotting_params) && !is.null(env_var_labels) && length(env_var_labels) > 0 && !is.null(names(env_var_labels))) {
  env_loadings <- env_loadings %>%
    mutate(var = factor(var, levels = names(env_var_labels), labels = unname(env_var_labels)))
}
env_summary <- env_loadings %>%
  group_by(canonical_direction, var) %>%
  summarise(mean = mean(value, na.rm = TRUE), se = sd(value, na.rm = TRUE) / sqrt(n()), .groups = "drop")
p_env <- ggplot(env_loadings, aes(x = var, y = value, color = factor(canonical_direction))) +
  geom_point(alpha = 0.3, size = 1) +
  geom_point(data = env_summary, aes(x = var, y = mean), color = 'black', size = 0.9, shape = 16, inherit.aes = FALSE) +
  geom_errorbar(data = env_summary, aes(x = var, y = mean, ymin = mean - se, ymax = mean + se), width = 0, linewidth = 0.6, color = 'black', inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey60") +
  scale_color_manual(values = cov_palette, name = "Canonical direction") +
  facet_wrap(~ canonical_direction, labeller = labeller(canonical_direction = function(x) paste0("CD ", x))) +
  labs(x = "Variable", y = "Loading (mean ± se)") +
  # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  guides(color = "none")
ggsave(file.path(step4_dir, "env_loadings.jpg"), plot = p_env, width = 8, height = 5, dpi = 150)
verbose_print("Wrote env_loadings.jpg", verbose = opt$verbose)

# ---- (3) Abundance loadings on tree (one per CD, shared colormap, phylum ring, mean loading) ----
if (opt$tax_level == "OTU") {
  abu_loadings <- read.csv(file.path(step2_dir, "abundance_loadings.csv"))
  abu_mean <- abu_loadings %>%
    group_by(canonical_direction, var) %>%
    summarise(mean_loading = mean(value, na.rm = TRUE), .groups = "drop")
  composition_final <- file.path(data_path, "16S", "GG2", opt$perc_identity, "final")
  tree_path <- file.path(composition_final, "tree.nwk")
  taxonomy_path <- file.path(composition_final, "taxonomy.csv")
  if (!file.exists(tree_path)) stop(paste("Tree file does not exist:", tree_path))
  if (!file.exists(taxonomy_path)) stop(paste("Taxonomy file does not exist:", taxonomy_path))

  tree <- read.tree(tree_path)
  taxonomy <- read_csv(taxonomy_path, show_col_types = FALSE)
  features <- unique(abu_mean$var)
  if (!all(features %in% tree$tip.label)) {
    missing <- setdiff(features, tree$tip.label)
    n_missing <- length(missing)
    if (n_missing > 0) {
      msg <- paste(n_missing, "out of", length(features), 
                   "feature(s) in loadings are missing from the tree.",
                   "Including:", paste(missing[1:10], collapse = ", ", "..."))
      if (n_missing < 10) {
        msg <- paste0(msg, " Missing feature(s): ", paste(missing, collapse = ", "))
      }
      stop(msg)
    }
  }
  tree_subset <- keep.tip(tree, features)
  if (!"Phylum" %in% names(taxonomy)) {
    stop("Taxonomy file does not contain a 'Phylum' column. Check to make sure GG2 clean-up was run.")
  }
  top_phyla <- NULL
  top_phyla <- taxonomy %>%
    filter(Feature_ID %in% tree_subset$tip.label) %>%
    select(Feature_ID, Phylum) %>%
    count(Phylum, sort = TRUE) %>%
    head(20) %>% pull(Phylum)
  taxonomy <- taxonomy %>%
    mutate(Phylum20 = ifelse(Phylum %in% top_phyla, Phylum, "Other")) %>%
    mutate(Phylum20 = factor(Phylum20, levels = c(top_phyla, "Other")))
  phylum_df <- taxonomy %>% filter(Feature_ID %in% tree_subset$tip.label) %>%
    select(Feature_ID, Phylum20) %>% tibble::column_to_rownames("Feature_ID")
  n_top <- length(top_phyla)
  phylum_cols <- c(hue_pal(l = 70, c = 100)(n_top), "grey80")

  load_mat <- abu_mean %>%
    pivot_wider(names_from = canonical_direction, values_from = mean_loading, names_prefix = "CD") %>%
    tibble::column_to_rownames("var")
  load_mat <- load_mat[tree_subset$tip.label, , drop = FALSE]
  maxmin <- max(abs(load_mat), na.rm = TRUE)
  if (maxmin <= 0) maxmin <- 0.01
  pal <- RColorBrewer::brewer.pal(11, "RdBu")
  cds <- colnames(load_mat)
  for (cd in cds) {
    m_cd <- load_mat[, cd, drop = TRUE]
    m_cd <- matrix(m_cd, ncol = 1)
    rownames(m_cd) <- rownames(load_mat)
    colnames(m_cd) <- cd
    circ <- ggtree(tree_subset, layout = "circular", size = 0.1)
    circ <- gheatmap(circ, phylum_df, offset = 0.02, width = 0.05, colnames = FALSE, color = NA) +
      scale_fill_manual(values = phylum_cols, breaks = c(top_phyla, "Other"), name = "Phyla")
    p_cd <- circ + new_scale_fill()
    df_cd <- as.data.frame(m_cd)
    if (opt$loadings_style == "binary") {
      df_cd[[1]] <- ifelse(df_cd[[1]] > 0, "positive", ifelse(df_cd[[1]] < 0, "negative", NA))
      p_cd <- gheatmap(p_cd, df_cd, offset = 22.5, width = 0.1, colnames = FALSE, color = NA) +
        scale_fill_manual(values = c(negative = "blue", positive = "red"), na.value = "grey80", name = "Loading")
    } else {
      p_cd <- gheatmap(p_cd, df_cd, offset = 22.5, width = 0.1, colnames = FALSE, color = NA) +
        scale_fill_gradientn(colours = pal, limits = c(-maxmin, maxmin), name = "Mean loading")
    }
    p_cd <- p_cd + ggtitle(cd)
    ggsave(file.path(step4_dir, paste0("abundance_loadings_tree_", tolower(cd), ".jpg")), plot = p_cd, width = 10, height = 8, dpi = 150)
  }
  verbose_print(paste("Wrote tree plots to", step4_dir), verbose = opt$verbose)
}

verbose_print("Step 4 done.", verbose = opt$verbose)
