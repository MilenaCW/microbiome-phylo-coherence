#!/usr/bin/env Rscript
# 06_consenTRAIT.R — consenTRAIT analysis on CCA OTU-level abundance loadings using consentrait_signed().
# Usage: Rscript 06_consenTRAIT.R --config <path> --perc_identity <p> [--verbose]

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) != 1) stop("Run with Rscript so --file= is available.")
script_path <- normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)
script_dir  <- dirname(script_path)
setup_path  <- normalizePath(file.path(script_dir, "..", "..", "setup.R"), winslash = "/", mustWork = TRUE)
source(setup_path)

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(readr)
  library(RColorBrewer)
  library(ape)
})
source(file.path(script_dir, "..", "functions", "CCA_functions.R"))
source(file.path(script_dir, "..", "functions", "consentrait_signed.R"))
source(file.path(script_dir, "..", "..", "utility_functions.R"))

option_list <- list(
  make_option(c("--config"),        type = "character", default = NULL,   help = "Path to CCA config R file.", metavar = "FILE"),
  make_option(c("--perc_identity"), type = "character", default = "0.90", help = "Perc identity [default %default]."),
  make_option(c("--verbose", "-v"), action = "store_true", default = FALSE, help = "Verbose output.")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$config) || !nzchar(opt$config)) stop("Missing required --config", call. = FALSE)
if (!file.exists(opt$config)) stop("Config file not found: ", opt$config, call. = FALSE)
cfg <- source(opt$config, local = TRUE)$value
if (!is.list(cfg)) stop("Config must evaluate to an R list.", call. = FALSE)

root         <- if (exists("REPO_ROOT")) REPO_ROOT else getwd()
results_path <- cfg$results_path
if (!grepl("^/", results_path)) results_path <- file.path(root, results_path)
data_path    <- cfg$data_path
if (!grepl("^/", data_path)) data_path <- file.path(root, data_path)

step2_dir <- file.path(results_path, opt$perc_identity, "OTU", "step2_loadings")
out_dir   <- file.path(results_path, opt$perc_identity, "OTU", "step6_consentrait")
tree_file <- file.path(data_path, "16S", "GG2", opt$perc_identity, "final", "tree.nwk")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load inputs ----
verbose_print("Loading abundance loadings...", verbose = opt$verbose)
abundance_loadings_file <- file.path(step2_dir, "abundance_loadings.csv")
if (!file.exists(abundance_loadings_file)) stop("Abundance loadings not found: ", abundance_loadings_file)
loadings <- read_csv(abundance_loadings_file, show_col_types = FALSE)
verbose_print(paste("Loaded", nrow(loadings), "loading rows;",
                    length(unique(loadings$var)), "unique OTUs,",
                    max(loadings$canonical_direction), "canonical directions,",
                    max(loadings$fold), "folds"), verbose = opt$verbose)

verbose_print("Loading phylogenetic tree...", verbose = opt$verbose)
if (!file.exists(tree_file)) stop("Tree file not found: ", tree_file)
tree <- read.tree(tree_file)
verbose_print(paste("Tree has", length(tree$tip.label), "tips"), verbose = opt$verbose)

# Require all loading OTUs to be in tree
loading_otus    <- unique(loadings$var)
missing_in_tree <- setdiff(loading_otus, tree$tip.label)
if (length(missing_in_tree) > 0) {
  stop("Error: ", length(missing_in_tree), " OTU(s) in loadings but not in tree. ",
       "Missing (first 20): ", paste(head(missing_in_tree, 20), collapse = ", "),
       if (length(missing_in_tree) > 20) paste0(" ... and ", length(missing_in_tree) - 20, " more"),
       call. = FALSE)
}
tree <- keep.tip(tree, loading_otus)
verbose_print(paste("Tree pruned to", length(tree$tip.label), "tips"), verbose = opt$verbose)

# ---- Run consentrait_signed for each (canonical_direction x fold) ----
cds       <- sort(unique(loadings$canonical_direction))
folds_all <- sort(unique(loadings$fold))
n_runs    <- length(cds) * length(folds_all)
verbose_print(paste("Running consentrait_signed for", length(cds), "canonical directions x",
                    length(folds_all), "folds =", n_runs, "runs"), verbose = opt$verbose)

tau_D_list      <- vector("list", n_runs)
null_tau_D_list <- vector("list", n_runs)
run_idx <- 0L

for (cd in cds) {
  for (f in folds_all) {
    run_idx <- run_idx + 1L
    verbose_print(paste0("  cd=", cd, " fold=", f, " (", run_idx, "/", n_runs, ")"),
                  verbose = opt$verbose)

    sub <- loadings %>% filter(canonical_direction == cd, fold == f)

    # Convert continuous loadings to +/-1 signed trait; drop exact zeros
    trait_vals_raw <- setNames(as.integer(sign(sub$value)), sub$var)
    n_zeros <- sum(trait_vals_raw == 0L)
    if (n_zeros > 0)
      warning("cd=", cd, " fold=", f, ": dropping ", n_zeros, " OTUs with loading == 0")
    trait_vals <- trait_vals_raw[trait_vals_raw != 0L]

    tree_sub <- keep.tip(tree, names(trait_vals))

    res <- consentrait_signed(
      tree                 = tree_sub,
      trait_values         = trait_vals,
      frac_consensus       = 0.9,
      n_shuffles           = 100L,
      singleton_depth_frac = 0.5,
      weight_clades        = TRUE,
      seed                 = 1L
    )

    p_val <- mean(res$null_tau_D >= res$tau_D)

    tau_D_list[[run_idx]] <- data.frame(
      canonical_direction = cd,
      fold    = f,
      tau_D   = res$tau_D,
      p_value = p_val
    )

    null_tau_D_list[[run_idx]] <- data.frame(
      canonical_direction = cd,
      fold       = f,
      shuffle_id = seq_along(res$null_tau_D),
      tau_D      = res$null_tau_D
    )

    # Save per-cell clades and null_clades
    cell_dir <- file.path(out_dir, paste0("cd", cd), paste0("fold", f))
    dir.create(cell_dir, recursive = TRUE, showWarnings = FALSE)
    write_csv(res$clades,      file.path(cell_dir, "clades.csv"))
    write_csv(res$null_clades, file.path(cell_dir, "null_clades.csv"))
  }
}

# ---- Bind and save aggregated CSVs ----
verbose_print("Saving aggregated CSVs...", verbose = opt$verbose)
tau_D_all      <- bind_rows(tau_D_list)
null_tau_D_all <- bind_rows(null_tau_D_list)

write_csv(tau_D_all,      file.path(out_dir, "tau_D.csv"))
write_csv(null_tau_D_all, file.path(out_dir, "null_tau_D.csv"))
verbose_print(paste("tau_D.csv:", nrow(tau_D_all), "rows"), verbose = opt$verbose)
verbose_print(paste("null_tau_D.csv:", nrow(null_tau_D_all), "rows"), verbose = opt$verbose)

# ---- Diagnostic plots ----
verbose_print("Generating plots...", verbose = opt$verbose)

# Read back per-cell clades to build aggregated frame for plotting
clades_all <- map_dfr(cds, function(cd) {
  map_dfr(folds_all, function(f) {
    p <- file.path(out_dir, paste0("cd", cd), paste0("fold", f), "clades.csv")
    read_csv(p, show_col_types = FALSE) %>%
      mutate(canonical_direction = cd, fold = f)
  })
})

# cluster_sizes: observed coherent clades (type == "clade")
clade_rows <- clades_all %>% filter(type == "clade")
if (nrow(clade_rows) > 0) {
  cd_breaks <- sort(unique(clades_all$canonical_direction))

  # Subtitle: % singletons per canonical direction (across all folds)
  singleton_pct_labels <- clades_all %>%
    group_by(canonical_direction) %>%
    summarise(pct = round(100 * mean(type == "singleton"), 1), .groups = "drop") %>%
    mutate(label = paste0(canonical_direction, "=", pct, "%")) %>%
    pull(label)
  subtitle_text <- paste(
    "% singletons by CD:",
    paste(
      sapply(split(singleton_pct_labels, ceiling(seq_along(singleton_pct_labels) / 3)),
             paste, collapse = ", "),
      collapse = "\n"
    )
  )

  p_sizes <- clade_rows %>%
    ggplot(aes(x = canonical_direction, y = size, color = factor(canonical_direction))) +
    geom_point(alpha = 0.3, position = position_jitter(width = 0.2, seed = 1)) +
    scale_x_continuous(breaks = cd_breaks) +
    scale_color_manual(values = cov_palette) +
    theme_minimal(base_size = 10) +
    labs(x = "canonical direction", y = "cluster size", subtitle = subtitle_text) +
    theme(legend.position = "none")
  ggsave(file.path(out_dir, "cluster_sizes.jpeg"), plot = p_sizes, width = 5, height = 5)
  verbose_print("Saved cluster_sizes.jpeg", verbose = opt$verbose)
}

# tauD: observed tau_D overlaid on null distribution
if (nrow(tau_D_all) > 0) {
  tau_breaks <- sort(unique(tau_D_all$canonical_direction))
  y_max <- max(c(tau_D_all$tau_D, null_tau_D_all$tau_D), na.rm = TRUE)
  p_tau <- ggplot() +
    geom_point(data = null_tau_D_all,
               aes(x = canonical_direction, y = tau_D),
               color = "grey", alpha = 0.1, position = position_jitter(width = 0.2, seed = 2)) +
    geom_point(data = tau_D_all,
               aes(x = canonical_direction, y = tau_D, color = factor(canonical_direction)),
               alpha = 0.25, position = position_jitter(width = 0.2, seed = 1)) +
    scale_x_continuous(breaks = tau_breaks) +
    scale_color_manual(values = cov_palette) +
    ylim(0, y_max) +
    theme_minimal(base_size = 10) +
    labs(x = "canonical direction", y = "mean genetic depth (tau_D)") +
    theme(legend.position = "none")
  ggsave(file.path(out_dir, "tauD.jpeg"), plot = p_tau, width = 5, height = 5)
  verbose_print("Saved tauD.jpeg", verbose = opt$verbose)
}

verbose_print("CCA consenTRAIT analysis complete.", verbose = opt$verbose)
