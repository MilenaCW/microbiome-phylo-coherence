# Plot abundance loadings on tree for one dataset, perc_identity, and canonical direction (manuscript).
# Based on CCA/scripts/04_model_performance_plots.R (OTU tree + phylum ring + loadings).
#
# Usage (from repo root):
#   Rscript code/manuscript_plotting/plot_abu_loadings.R --dataset soil --perc_identity 0.99 --direction 1
#
# Output: manuscript/<dataset_short>_abu_loadings_tree_<perc_short>_cd<direction>.pdf

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggtree)
  library(ggnewscale)
  library(ape)
})

stopf <- function(...) stop(sprintf(...), call. = FALSE)

# Manuscript palettes
phyla_colors_top10 <- c(
  "#7e9652", "#b8c98a", "#4f807b", "#5f97b8", "#4f709b",
  "#7067a5", "#9b6c8f", "#a65c71", "#b04b53", "#b93a34"
)
phyla_other <- "#cdcccc"
# loading_positive <- "#dd8286"
loading_positive <- "white"
# loading_negative <- "#2a254a"
loading_negative <- "black"

parse_cli_args <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (!requireNamespace("optparse", quietly = TRUE)) {
    stopf("R package 'optparse' is required. Install it with: install.packages('optparse')")
  }

  option_list <- list(
    optparse::make_option(
      "--dataset",
      type = "character",
      default = NA_character_,
      help = "Dataset: soil, ocean, or wwtp"
    ),
    optparse::make_option(
      "--perc_identity",
      type = "character",
      default = "0.90",
      help = "Perc identity (e.g. 0.99) [default %default]"
    ),
    optparse::make_option(
      "--direction",
      type = "integer",
      default = NA_integer_,
      help = "Canonical direction number (e.g. 1)"
    )
  )

  p <- optparse::OptionParser(option_list = option_list)
  a <- optparse::parse_args(p, args = argv, positional_arguments = FALSE)

  list(
    dataset = a$dataset,
    perc_identity = a$perc_identity,
    direction = a$direction
  )
}

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

main <- function() {
  args <- parse_cli_args()
  if (is.na(args$dataset) || !nzchar(trimws(args$dataset))) {
    stop("Required flag --dataset must be one of: soil, ocean, wwtp", call. = FALSE)
  }
  if (is.na(args$direction) || args$direction < 1L) {
    stop("Required flag --direction must be a positive integer (canonical direction number).", call. = FALSE)
  }

  ds <- trimws(tolower(args$dataset))
  if (!ds %in% c("soil", "ocean", "wwtp")) {
    stopf("Unknown --dataset '%s'. Must be one of: soil, ocean, wwtp", args$dataset)
  }
  perc <- args$perc_identity
  cd_num <- as.integer(args$direction)
  root <- get_repo_root()

  # Paths: results from dataset folder; data_path from CCA config when available
  results_path <- file.path(root, ds, "results", "CCA")
  step2_dir <- file.path(results_path, perc, "OTU", "step2_loadings")

  data_path <- NULL
  config_path <- file.path(root, "code", "CCA", "config", paste0(ds, ".R"))
  if (file.exists(config_path)) {
    cfg <- source(config_path, local = TRUE)$value
    if (is.list(cfg) && !is.null(cfg$data_path)) {
      data_path <- cfg$data_path
      if (!grepl("^/", data_path)) data_path <- file.path(root, data_path)
    }
  }
  if (is.null(data_path) || !nzchar(data_path)) {
    data_path <- file.path(root, ds, "data", "processed_data")
  }

  composition_final <- file.path(data_path, "16S", "GG2", perc, "final")
  tree_path <- file.path(composition_final, "tree.nwk")
  taxonomy_path <- file.path(composition_final, "taxonomy.csv")

  abu_file <- file.path(step2_dir, "abundance_loadings.csv")
  if (!file.exists(abu_file)) {
    stopf("Abundance loadings not found: %s\nRun CCA step 2 for this dataset/perc_identity/OTU first.", abu_file)
  }
  if (!file.exists(tree_path)) stopf("Tree file does not exist: %s", tree_path)
  if (!file.exists(taxonomy_path)) stopf("Taxonomy file does not exist: %s", taxonomy_path)

  abu_loadings <- read.csv(abu_file)
  abu_mean <- abu_loadings %>%
    group_by(canonical_direction, var) %>%
    summarise(mean_loading = mean(value, na.rm = TRUE), .groups = "drop")

  if (!cd_num %in% abu_mean[["canonical_direction"]]) {
    stopf("Canonical direction %d not found in loadings. Available: %s",
          cd_num, paste(sort(unique(abu_mean$canonical_direction)), collapse = ", "))
  }

  tree <- read.tree(tree_path)
  taxonomy <- read_csv(taxonomy_path, show_col_types = FALSE)
  features <- unique(abu_mean$var)

  if (!all(features %in% tree$tip.label)) {
    missing <- setdiff(features, tree$tip.label)
    n_missing <- length(missing)
    stopf("%d feature(s) in loadings are missing from the tree. First few: %s",
          n_missing, paste(head(missing, 10), collapse = ", "))
  }

  tree_subset <- keep.tip(tree, features)
  if (!"Phylum" %in% names(taxonomy)) {
    stop("Taxonomy file does not contain a 'Phylum' column. Check that GG2 clean-up was run.", call. = FALSE)
  }

  # Top 10 phyla (manuscript)
  top_phyla <- taxonomy %>%
    filter(Feature_ID %in% tree_subset$tip.label) %>%
    count(Phylum, sort = TRUE) %>%
    head(10) %>%
    pull(Phylum)

  taxonomy <- taxonomy %>%
    mutate(Phylum10 = if_else(Phylum %in% top_phyla, Phylum, "Other")) %>%
    mutate(Phylum10 = factor(Phylum10, levels = c(top_phyla, "Other")))

  phylum_df <- taxonomy %>%
    filter(Feature_ID %in% tree_subset$tip.label) %>%
    select(Feature_ID, Phylum10) %>%
    tibble::column_to_rownames("Feature_ID")

  n_top <- length(top_phyla)
  phylum_cols <- c(phyla_colors_top10[seq_len(n_top)], phyla_other)

  # Loadings matrix for the requested canonical direction only
  load_mat <- abu_mean %>%
    filter(canonical_direction == cd_num) %>%
    pivot_wider(names_from = canonical_direction, values_from = mean_loading, names_prefix = "CD") %>%
    tibble::column_to_rownames("var")
  load_mat <- load_mat[tree_subset$tip.label, , drop = FALSE]

  cd_col <- colnames(load_mat)[1L]
  m_cd <- load_mat[, cd_col, drop = TRUE]
  m_cd <- matrix(m_cd, ncol = 1)
  rownames(m_cd) <- rownames(load_mat)
  colnames(m_cd) <- cd_col

  df_cd <- as.data.frame(m_cd)
  df_cd[[1]] <- ifelse(df_cd[[1]] > 0, "positive", ifelse(df_cd[[1]] < 0, "negative", NA))

  circ <- ggtree(tree_subset, layout = "circular", size = 0.05)
  circ <- gheatmap(circ, phylum_df, offset = 0.02, width = 0.1, colnames = FALSE, color = NA) +
    scale_fill_manual(values = phylum_cols, breaks = c(top_phyla, "Other"), name = "Phyla", guide = guide_legend(ncol = 2))

  p_cd <- circ + new_scale_fill()
  p_cd <- gheatmap(p_cd, df_cd, offset = 22.5, width = 0.25, colnames = FALSE, color = NA) +
    scale_fill_manual(
      values = c(negative = loading_negative, positive = loading_positive),
      breaks = c("negative", "positive"),
      labels = c(expression(a[i] < 0), expression(a[i] > 0)),
      na.value = phyla_other,
      name = "Loading"
    ) +
    theme(
      text = element_text(size = 10),
      legend.key.size = unit(0.25, "cm"),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 10)
    )

  out_dir <- file.path(root, "manuscript", ds)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  perc_short <- paste0("p", sub("^[^.]*\\.?", "", perc))
  base_name <- sprintf("%s_abu_loadings_tree_%s_cd%d", ds, perc_short, cd_num)

  # With legends
  out_path <- file.path(out_dir, paste0(base_name, ".pdf"))
  ggsave(out_path, plot = p_cd, width = 6, height = 3)
  message("Saved: ", out_path)

  # Without legends: minimal padding so circle fills the device
  p_cd_nolegend <- p_cd +
    theme(
      legend.position = "none",
      plot.margin = margin(-3, -3, -4.5, -4.5, "mm"),
      panel.spacing = unit(0, "cm")
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0))

  out_path_nolegend <- file.path(out_dir, paste0(base_name, "_nolegend.pdf"))
  ggsave(out_path_nolegend, plot = p_cd_nolegend, width = 2, height = 2)
  message("Saved: ", out_path_nolegend)
}

main()
