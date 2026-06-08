# SI_example_alignment.R — Demonstrate sign-flip alignment of loading vectors across
# cross-validation folds using a self-contained synthetic example.
# 5 environmental variables, 10 folds, 1 canonical direction.
# Folds that were sign-corrected are shown in black; others in grey.
#
# Usage (from repo root):
#   Rscript code/manuscript_plotting/SI_example_alignment.R
#   Rscript code/manuscript_plotting/SI_example_alignment.R --seed 42

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

stopf <- function(...) stop(sprintf(...), call. = FALSE)

parse_cli_args <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (!requireNamespace("optparse", quietly = TRUE))
    stopf("R package 'optparse' is required. Install it with: install.packages('optparse')")
  option_list <- list(
    optparse::make_option("--seed", type = "integer", default = 42L,
                          help = "RNG seed for synthetic data [default %default]"),
    optparse::make_option("--root", type = "character", default = NA_character_,
                          help = "Repo root path (optional; if omitted, uses code/setup.R)")
  )
  a <- optparse::parse_args(optparse::OptionParser(option_list = option_list), args = argv)
  list(seed = a$seed, root = a$root)
}

get_repo_root <- function() {
  args     <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) != 1)
    stop("This script must be run with Rscript so --file= is available.", call. = FALSE)
  script_path <- normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)
  script_dir  <- dirname(script_path)
  setup_path  <- normalizePath(file.path(script_dir, "..", "setup.R"), winslash = "/", mustWork = TRUE)
  source(setup_path, local = TRUE)
  get("REPO_ROOT", envir = .GlobalEnv)
}

add_flip_flag <- function(df, flips_df) {
  df %>%
    left_join(flips_df %>% mutate(flipped = TRUE),
              by = c("canonical_direction", "fold")) %>%
    mutate(flipped = tidyr::replace_na(flipped, FALSE))
}

make_panel <- function(df, title) {
  shared_theme <- theme_minimal() + theme(
    text            = element_text(size = 10),
    plot.title      = element_text(size = 10),
    axis.title      = element_text(size = 10),
    axis.text       = element_text(size = 10),
    axis.text.x     = element_text(size = 10),
    legend.text     = element_text(size = 10),
    strip.text      = element_text(face = "plain", size = 10),
    legend.position = "bottom"
  )
  ggplot(df, aes(x = var, y = value, group = fold, colour = flipped)) +
    geom_line(alpha = 0.8) +
    geom_point(size = 1.5) +
    scale_colour_manual(
      values = c("FALSE" = "grey60", "TRUE" = "black"),
      labels = c("FALSE" = "same orientation", "TRUE" = "sign-flipped"),
      name   = NULL
    ) +
    ylim(-1, 1) +
    labs(x = "environmental variable", y = "loading", title = title) +
    shared_theme
}

# Generate synthetic loading vectors: 5 variables, 10 folds, 1 canonical direction.
# Folds 3, 5, 8 have their signs deliberately inverted before alignment.
generate_synthetic_loadings <- function(seed = 42L) {
  set.seed(seed)
  n_vars     <- 5L
  n_folds    <- 10L
  true_v     <- c(0.4, 0.8, 0.3, -0.5, 0.6)
  true_v     <- true_v / sqrt(sum(true_v^2))
  flip_folds <- c(3L, 5L, 8L)
  fold_rows <- lapply(seq_len(n_folds), function(f) {
    noisy <- true_v + rnorm(n_vars, sd = 0.12)
    noisy <- noisy / sqrt(sum(noisy^2))
    if (f %in% flip_folds) noisy <- -noisy
    data.frame(
      canonical_direction = 1L,
      fold  = f,
      var   = paste0("V", seq_len(n_vars)),
      value = noisy,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, fold_rows)
}

main <- function() {
  args <- parse_cli_args()
  root <- args$root
  if (is.na(root) || !nzchar(trimws(root))) root <- get_repo_root()
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)

  cca_fn <- file.path(root, "code", "CCA", "functions", "CCA_functions.R")
  if (!file.exists(cca_fn)) stopf("CCA_functions.R not found: %s", cca_fn)
  source(cca_fn, local = TRUE)

  loadings_raw <- generate_synthetic_loadings(seed = args$seed)

  alignment    <- align_loading_vectors(loadings_raw, reference_method = "first_fold")
  loadings_aln <- alignment$aligned_loadings

  flips_df <- tibble::enframe(
    alignment$flipped_folds,
    name  = "canonical_direction",
    value = "fold"
  ) %>%
    tidyr::unnest_longer(fold) %>%
    dplyr::mutate(canonical_direction = as.integer(canonical_direction))

  orig    <- add_flip_flag(loadings_raw, flips_df)
  aligned <- add_flip_flag(loadings_aln, flips_df)

  p_orig    <- make_panel(orig,    "before alignment")
  p_aligned <- make_panel(aligned, "after alignment")

  p <- (p_orig | p_aligned) +
    plot_layout(guides = "collect") &
    theme(
      legend.position    = "bottom",
      legend.box.spacing = unit(2, "pt"),
      legend.margin      = margin(0, 0, 0, 0)
    )

  out_dir <- file.path(root, "manuscript", "SI")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(out_dir, "SI_example_alignment.pdf")
  ggsave(out_file, plot = p, width = 4, height = 2.5, units = "in")
  message("Saved: ", out_file)
}

main()
