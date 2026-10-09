# Rough exploratory comparison of CCA canonical directions across soil hemisphere splits.
# For each pair of fits (North, South, Full), aligns mean env loadings per canonical
# direction and builds a CD x CD sign-corrected correlation matrix — not just the
# diagonal (CD_i vs CD_i), since canonical dimension rank order isn't guaranteed to
# line up across independently-fit models. Sign-correction logic follows
# CCA/scripts/05_crosstaxonomic_compare.R:143-171.
#
# One-off exploratory script: not wired into the numbered CCA pipeline.
#
# Usage (from repo root):
#   Rscript code/manuscript_plotting/hemisphere_alignment_explore.R --perc_identity 0.90

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(tidyr)
})

stopf <- function(...) stop(sprintf(...), call. = FALSE)

parse_cli_args <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (!requireNamespace("optparse", quietly = TRUE)) {
    stopf("R package 'optparse' is required. Install it with: install.packages('optparse')")
  }
  option_list <- list(
    optparse::make_option(
      "--perc_identity",
      type = "character",
      default = "0.90",
      help = "Perc identity (e.g. 0.90) [default %default]"
    )
  )
  p <- optparse::OptionParser(option_list = option_list)
  a <- optparse::parse_args(p, args = argv, positional_arguments = FALSE)
  list(perc_identity = a$perc_identity)
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

# Reads env_loadings.csv (canonical_direction, fold, var, value) and collapses
# to one mean loading vector per canonical direction, matched by var name.
read_env_loadings <- function(env_file, group_label) {
  if (!file.exists(env_file)) {
    stopf("env_loadings.csv not found for '%s': %s", group_label, env_file)
  }
  read_csv(env_file, show_col_types = FALSE) %>%
    group_by(canonical_direction, var) %>%
    summarise(mean = mean(value, na.rm = TRUE), .groups = "drop") %>%
    mutate(group = group_label)
}

# Builds a CD_a x CD_b sign-corrected Pearson correlation matrix between two groups'
# mean env loading vectors. Matrix entries are |cor| (sign is arbitrary for CCA
# weight vectors); `signed` keeps the raw correlation for reference.
build_alignment_matrix <- function(loadings_a, loadings_b, label_a, label_b) {
  cds_a <- sort(unique(loadings_a$canonical_direction))
  cds_b <- sort(unique(loadings_b$canonical_direction))

  vars_a <- loadings_a %>% filter(canonical_direction == cds_a[1]) %>% pull(var) %>% sort()
  vars_b <- loadings_b %>% filter(canonical_direction == cds_b[1]) %>% pull(var) %>% sort()
  if (!setequal(vars_a, vars_b)) {
    stopf("Env variables do not match between '%s' and '%s'.", label_a, label_b)
  }

  results <- expand.grid(cd_a = cds_a, cd_b = cds_b) %>%
    rowwise() %>%
    mutate(
      vec_a = list((loadings_a %>% filter(canonical_direction == cd_a) %>% arrange(var))$mean),
      vec_b = list((loadings_b %>% filter(canonical_direction == cd_b) %>% arrange(var))$mean),
      cor_raw = cor(unlist(vec_a), unlist(vec_b)),
      cor_abs = abs(cor_raw)
    ) %>%
    ungroup() %>%
    select(cd_a, cd_b, cor_raw, cor_abs) %>%
    mutate(group_a = label_a, group_b = label_b)

  results
}

plot_alignment_heatmap <- function(mat_df, label_a, label_b, out_path) {
  p <- mat_df %>%
    ggplot(aes(x = factor(cd_b), y = factor(cd_a), fill = cor_abs)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", cor_raw)), size = 3, color = "black") +
    scale_fill_gradient(low = "white", high = "#332288", limits = c(0, 1), name = "|r|") +
    labs(
      title = sprintf("%s vs %s: env loading alignment by canonical direction", label_a, label_b),
      x = paste(label_b, "CD"),
      y = paste(label_a, "CD")
    ) +
    theme_minimal(base_size = 10) +
    coord_fixed()
  ggsave(out_path, plot = p, width = 5, height = 4.5, dpi = 150)
  message("Saved: ", out_path)
}

summarize_best_matches <- function(mat_df, label_a, label_b) {
  mat_df %>%
    group_by(cd_a) %>%
    slice_max(cor_abs, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(
      comparison = sprintf("%s vs %s", label_a, label_b),
      from_cd = cd_a,
      best_match_cd = cd_b,
      r = round(cor_raw, 3)
    )
}

main <- function() {
  args <- parse_cli_args()
  root <- get_repo_root()
  perc <- args$perc_identity

  paths <- list(
    North = file.path(root, "soil", "results", "CCA_hemisphere", "north", perc, "OTU", "step2_loadings", "env_loadings.csv"),
    South = file.path(root, "soil", "results", "CCA_hemisphere", "south", perc, "OTU", "step2_loadings", "env_loadings.csv"),
    Full  = file.path(root, "soil", "results", "CCA", perc, "OTU", "step2_loadings", "env_loadings.csv")
  )

  loadings <- Map(read_env_loadings, paths, names(paths))

  # Cross-reference: test-set canonical correlation strength by CD, per group
  # (flags CDs being compared that had weak/near-zero out-of-sample correlation).
  corr_paths <- list(
    North = file.path(root, "soil", "results", "CCA_hemisphere", "north", perc, "OTU", "step2_loadings", "correlations_per_fold.csv"),
    South = file.path(root, "soil", "results", "CCA_hemisphere", "south", perc, "OTU", "step2_loadings", "correlations_per_fold.csv"),
    Full  = file.path(root, "soil", "results", "CCA", perc, "OTU", "step2_loadings", "correlations_per_fold.csv")
  )
  test_cor_summary <- bind_rows(Map(function(path, label) {
    read_csv(path, show_col_types = FALSE) %>%
      group_by(canonical_direction) %>%
      summarise(mean_test_cor = round(mean(test, na.rm = TRUE), 3), .groups = "drop") %>%
      mutate(group = label)
  }, corr_paths, names(corr_paths)))

  out_dir <- file.path(root, "soil", "results", "CCA_hemisphere", "alignment_explore")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  pairs <- list(
    c("North", "South"),
    c("North", "Full"),
    c("South", "Full")
  )

  all_best_matches <- list()
  for (pair in pairs) {
    label_a <- pair[1]
    label_b <- pair[2]
    mat_df <- build_alignment_matrix(loadings[[label_a]], loadings[[label_b]], label_a, label_b)
    out_path <- file.path(out_dir, sprintf("alignment_%s_vs_%s.jpg", tolower(label_a), tolower(label_b)))
    plot_alignment_heatmap(mat_df, label_a, label_b, out_path)
    all_best_matches[[paste(label_a, label_b)]] <- summarize_best_matches(mat_df, label_a, label_b)
  }

  best_matches <- bind_rows(all_best_matches)
  message("\nBest-matching canonical direction per comparison (env loadings):")
  print(as.data.frame(best_matches))

  message("\nMean out-of-sample (test) canonical correlation by CD and group:")
  print(as.data.frame(test_cor_summary %>% pivot_wider(names_from = group, values_from = mean_test_cor)))

  write_csv(best_matches, file.path(out_dir, "best_matches_summary.csv"))
  write_csv(test_cor_summary, file.path(out_dir, "test_correlation_by_cd_and_group.csv"))
  message("\nSaved summary tables to: ", out_dir)
}

main()
