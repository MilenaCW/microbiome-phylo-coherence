suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(patchwork)
  library(optparse)
})

stopf <- function(...) stop(sprintf(...), call. = FALSE)

TAX_LEVELS_ORDER <- c("OTU", "Species", "Genus", "Family", "Order", "Class", "Phylum")
TAX_ABBR         <- c(OTU = "OTU", Species = "Sp", Genus = "Ge",
                      Family = "Fa", Order = "Or", Class = "Cl", Phylum = "Ph")
DATASET_COLORS   <- c(soil = "#8f723d", ocean = "#90afa7")
DATASET_LABELS   <- c(soil = "soil", ocean = "ocean")

# ---------------------------------------------------------------------------
get_repo_root <- function() {
  args     <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) != 1) stop("Run with Rscript so --file= is available.", call. = FALSE)
  script_path <- normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)
  # two levels up: scratch/ -> manuscript_plotting/ -> code/
  setup_path <- normalizePath(file.path(dirname(script_path), "..", "setup.R"),
                              winslash = "/", mustWork = TRUE)
  source(setup_path, local = TRUE)
  get("REPO_ROOT", envir = .GlobalEnv)
}

# ---------------------------------------------------------------------------
parse_cli_args <- function(argv = commandArgs(trailingOnly = TRUE)) {
  option_list <- list(
    make_option("--perc_identity", type = "character", default = "0.90",
                help = "Percent identity [default %default]"),
    make_option("--root", type = "character", default = NA_character_,
                help = "Repo root path (inferred from script location if omitted)"),
    make_option("--n_cd_soil",  type = "integer", default = 4L,
                help = "Number of canonical directions to keep for soil [default %default]"),
    make_option("--n_cd_ocean", type = "integer", default = 3L,
                help = "Number of canonical directions to keep for ocean [default %default]")
  )
  a <- parse_args(OptionParser(option_list = option_list), args = argv)
  list(
    perc = a$perc_identity,
    root = a$root,
    n_cd_soil = a$n_cd_soil,
    n_cd_ocean = a$n_cd_ocean
  )
}

# ---------------------------------------------------------------------------
discover_tax_dirs <- function(base) {
  all_dirs <- list.dirs(base, full.names = FALSE, recursive = FALSE)
  found <- character(0)
  for (d in all_dirs) {
    if (d == "" || d == "crosstax_compare" || d == "Domain") next
    if (file.exists(file.path(base, d, "step2_loadings", "env_loadings.csv")))
      found <- c(found, d)
  }
  intersect(TAX_LEVELS_ORDER, found)
}

# ---------------------------------------------------------------------------
clamp01 <- function(x) pmin(1, pmax(0, x))

# ---------------------------------------------------------------------------
make_pair_key <- function(tax_dirs) {
  tax_dirs <- intersect(TAX_LEVELS_ORDER, tax_dirs)
  if (length(tax_dirs) < 2) stop("Need at least 2 taxonomic levels.", call. = FALSE)

  cm <- combn(tax_dirs, 2)
  tibble(
    tax_i = factor(cm[1, ], levels = TAX_LEVELS_ORDER),
    tax_j = factor(cm[2, ], levels = TAX_LEVELS_ORDER)
  ) %>%
    mutate(
      pair_id    = paste(as.character(tax_i), as.character(tax_j), sep = "__"),
      pair_label = paste(TAX_ABBR[as.character(tax_i)],
                         TAX_ABBR[as.character(tax_j)],
                         sep = "\u2013")
    )
}

# ---------------------------------------------------------------------------
angle_scale <- function() {
  scale_y_continuous(
    limits = c(0, pi / 2),
    expand = expansion(mult = c(0, 0.03)),
    breaks = c(0, pi / 4, pi / 2),
    labels = expression(0, pi/4, pi/2)
  )
}

# ---------------------------------------------------------------------------
theme_si <- function() {
  theme_minimal(base_size = 10) +
    theme(
      strip.text         = element_text(size = 10),
      axis.text.x        = element_text(size = 10),
      axis.text.y        = element_text(size = 10),
      axis.title         = element_text(size = 10),
      legend.text        = element_text(size = 10),
      legend.position    = "bottom",
      plot.margin        = margin(6, 6, 6, 6)
    )
}

# ---------------------------------------------------------------------------
# Adds an invisible phantom geom + scale_fill_manual to a panel so that
# patchwork's guides = "collect" can build a shared soil/ocean/null legend.
# Points are NA so nothing is drawn; override.aes controls legend appearance.
legend_layer <- function() {
  df <- data.frame(
    x   = NA_real_,
    y   = NA_real_,
    cat = factor(c("soil", "ocean", "null"), levels = c("soil", "ocean", "null"))
  )

  shared_guide <- guide_legend(
    override.aes = list(
      size   = 3,
      stroke = 0.75,
      shape  = 21,
      fill   = c("#8f723d", "#90afa7", "grey88"),
      color  = c("#8f723d", "#90afa7", "grey50")
    )
  )

  list(
    geom_point(
      data        = df,
      aes(x = x, y = y, color = cat, fill = cat),
      shape       = 21, size = 0, stroke = 0.75,
      inherit.aes = FALSE, show.legend = TRUE, na.rm = TRUE
    ),
    scale_color_manual(
      name   = NULL,                          # <-- same name on both scales
      values = c(soil = "#8f723d", ocean = "#90afa7", null = "grey50"),
      drop   = FALSE,
      guide  = shared_guide
    ),
    scale_fill_manual(
      name   = NULL,                          # <-- same name → merged legend
      values = c(soil = "#8f723d", ocean = "#90afa7", null = "grey88"),
      drop   = FALSE,
      guide  = shared_guide
    )
  )
}

# ---------------------------------------------------------------------------
load_dataset_inputs <- function(root, ds, perc) {
  base <- file.path(root, ds, "results", "CCA", perc)
  if (!dir.exists(base)) stopf("Base directory does not exist: %s", base)

  tax_dirs <- discover_tax_dirs(base)
  if (length(tax_dirs) < 2) stopf("Need >= 2 tax levels for %s/%s", ds, perc)

  message("  ", ds, ": ", paste(tax_dirs, collapse = ", "))

  real_raw <- bind_rows(lapply(tax_dirs, function(tax) {
    f <- file.path(base, tax, "step2_loadings", "env_loadings.csv")
    read_csv(f, show_col_types = FALSE) %>%
      mutate(tax_level = tax)
  }))

  null_parts <- lapply(tax_dirs, function(tax) {
    f <- file.path(base, tax, "step3_null", "null_env_loadings.csv")
    if (!file.exists(f)) {
      message("  No null file for ", ds, "/", tax, ", skipping")
      return(NULL)
    }
    read_csv(f, show_col_types = FALSE) %>%
      mutate(tax_level = tax)
  })

  null_parts <- Filter(Negate(is.null), null_parts)
  null_raw   <- if (length(null_parts) > 0) bind_rows(null_parts) else NULL

  real_raw <- real_raw %>%
    filter(tax_level %in% tax_dirs) %>%
    mutate(
      tax_level = factor(tax_level, levels = TAX_LEVELS_ORDER),
      canonical_direction = as.integer(canonical_direction),
      fold = as.integer(fold)
    )

  if (!is.null(null_raw)) {
    null_raw <- null_raw %>%
      filter(tax_level %in% tax_dirs) %>%
      mutate(
        tax_level = factor(tax_level, levels = TAX_LEVELS_ORDER),
        canonical_direction = as.integer(canonical_direction),
        fold = as.integer(fold),
        seed = as.integer(seed)
      )
  }

  list(
    ds = ds,
    perc = perc,
    base = base,
    tax_dirs = tax_dirs,
    pair_key = make_pair_key(tax_dirs),
    real_raw = real_raw,
    null_raw = null_raw
  )
}

# ---------------------------------------------------------------------------
compute_angles <- function(df, has_seed = FALSE) {
  # df has columns: canonical_direction, var, fold, value, tax_level (+ seed)
  shared <- setdiff(names(df), c("tax_level", "value"))

  df %>%
    inner_join(df, by = shared, suffix = c("_i", "_j"), relationship = "many-to-many") %>%
    filter(as.integer(tax_level_i) < as.integer(tax_level_j)) %>%
    rename(tax_i = tax_level_i, tax_j = tax_level_j) %>%
    group_by(canonical_direction, fold, tax_i, tax_j, across(any_of("seed"))) %>%
    summarise(
      dot    = sum(value_i * value_j, na.rm = TRUE),
      norm_i = sqrt(sum(value_i^2, na.rm = TRUE)),
      norm_j = sqrt(sum(value_j^2, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    filter(norm_i > 0, norm_j > 0) %>%
    mutate(
      pair_id = paste(as.character(tax_i), as.character(tax_j), sep = "__"),
      angle   = acos(clamp01(abs(dot / (norm_i * norm_j))))
    )
}

# ---------------------------------------------------------------------------
summarize_angles <- function(real_angles, null_angles = NULL, pair_key) {
  # ---------------- TRUE ----------------
  real_pair <- real_angles %>%
    group_by(canonical_direction, tax_i, tax_j, pair_id) %>%
    summarise(
      n_folds = n(),
      m_jp = mean(angle, na.rm = TRUE),
      sd_fold_within_pair = sd(angle, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(pair_key, by = c("tax_i", "tax_j", "pair_id"))

  real_overall <- real_angles %>%
    group_by(canonical_direction, fold) %>%
    summarise(M_jf = mean(angle, na.rm = TRUE), .groups = "drop") %>%
    group_by(canonical_direction) %>%
    summarise(
      n_folds = n(),
      M_j = mean(M_jf, na.rm = TRUE),
      sd_M_over_folds = sd(M_jf, na.rm = TRUE),
      .groups = "drop"
    )

  # ---------------- NULL ----------------
  null_pair_seed    <- NULL
  null_pair_summary <- NULL
  null_overall_seed <- NULL
  null_summary      <- NULL
  true_vs_null      <- NULL

  if (!is.null(null_angles) && nrow(null_angles) > 0) {
    null_pair_seed <- null_angles %>%
      group_by(canonical_direction, tax_i, tax_j, pair_id, seed) %>%
      summarise(
        m_jps = mean(angle, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      left_join(pair_key, by = c("tax_i", "tax_j", "pair_id"))

    null_pair_summary <- null_pair_seed %>%
      group_by(canonical_direction, tax_i, tax_j, pair_id, pair_label) %>%
      summarise(
        n_seeds = n(),
        null_pair_mean   = mean(m_jps, na.rm = TRUE),
        null_pair_median = median(m_jps, na.rm = TRUE),
        null_pair_q025   = quantile(m_jps, 0.025, na.rm = TRUE),
        null_pair_q975   = quantile(m_jps, 0.975, na.rm = TRUE),
        .groups = "drop"
      )

    null_overall_seed <- null_pair_seed %>%
      group_by(canonical_direction, seed) %>%
      summarise(
        M_js = mean(m_jps, na.rm = TRUE),
        .groups = "drop"
      )

    null_summary <- null_overall_seed %>%
      group_by(canonical_direction) %>%
      summarise(
        n_seeds = n(),
        null_mean   = mean(M_js, na.rm = TRUE),
        null_median = median(M_js, na.rm = TRUE),
        null_q025   = quantile(M_js, 0.025, na.rm = TRUE),
        null_q975   = quantile(M_js, 0.975, na.rm = TRUE),
        null_sd     = sd(M_js, na.rm = TRUE),
        .groups = "drop"
      )

    true_vs_null <- real_overall %>%
      left_join(null_summary, by = "canonical_direction") %>%
      rowwise() %>%
      mutate(
        z_like = (null_mean - M_j) / null_sd,
        p_empirical = {
          cd <- canonical_direction
          x_null <- null_overall_seed %>%
            filter(canonical_direction == cd) %>%
            pull(M_js)
          (1 + sum(x_null <= M_j)) / (length(x_null) + 1)
        }
      ) %>%
      ungroup()
  }

  list(
    real_pair = real_pair,
    real_overall = real_overall,
    null_pair_seed = null_pair_seed,
    null_pair_summary = null_pair_summary,
    null_overall_seed = null_overall_seed,
    null_summary = null_summary,
    true_vs_null = true_vs_null
  )
}

# ---------------------------------------------------------------------------
make_overall_panel <- function(res, ds) {
  ds_col <- DATASET_COLORS[[ds]]
  ds_lab <- DATASET_LABELS[[ds]]

  if (is.null(res$null_overall_seed) || nrow(res$null_overall_seed) == 0) {
    return(
      ggplot() +
        theme_void() +
        annotate("text", x = 0, y = 0, label = paste0(ds_lab, "\n(no null data)"), size = 4)
    )
  }

  ggplot(res$null_overall_seed, aes(x = factor(canonical_direction), y = M_js)) +
    geom_violin(
      draw_quantiles = c(0.5),  # add median line
      fill = "grey88", 
      color = "grey50", 
      linewidth = 0.25, width = 0.9, 
      scale = "width") +
    # geom_point(
    #   data = res$null_summary,
    #   aes(x = factor(canonical_direction), y = null_mean),
    #   inherit.aes = FALSE,
    #   shape = 21, fill = NA, color = "grey50", size = 3.0, stroke = 1.0
    # ) +
    geom_point(
      data = res$real_overall,
      aes(x = factor(canonical_direction), y = M_j),
      inherit.aes = FALSE,
      shape = 21, fill = ds_col, color = ds_col, size = 3.2, stroke = 0.3
    ) +
    angle_scale() +
    legend_layer() +
    labs(
      x = "canonical direction",
      y = "angle (rad)"
    ) +
    theme_si()
}

# ---------------------------------------------------------------------------
make_pair_panel <- function(res, ds, n_cd) {
  ds_col <- DATASET_COLORS[[ds]]
  ds_lab <- DATASET_LABELS[[ds]]

  if (is.null(res$null_pair_summary) || nrow(res$null_pair_summary) == 0) {
    return(
      ggplot() +
        theme_void() +
        annotate("text", x = 0, y = 0, label = paste0(ds_lab, "\n(no null pair data)"), size = 4)
    )
  }

  plot_null <- res$null_pair_summary %>%
    mutate(
      canonical_direction = factor(canonical_direction, levels = seq_len(n_cd)),
      pair_label = factor(pair_label, levels = unique(pair_label))
    )

  plot_true <- res$real_pair %>%
    mutate(
      canonical_direction = factor(canonical_direction, levels = seq_len(n_cd)),
      pair_label = factor(pair_label, levels = levels(plot_null$pair_label))
    )

  ggplot() +
    geom_linerange(
      data = plot_null,
      aes(x = pair_label, ymin = null_pair_q025, ymax = null_pair_q975),
      color = "grey50", linewidth = 0.5
    ) +
    geom_point(
      data = plot_null,
      aes(x = pair_label, y = null_pair_median),
      shape = 21, fill = "grey88", color = "grey50", size = 1.6, stroke = 0.5
    ) +
    geom_point(
      data = plot_true,
      aes(x = pair_label, y = m_jp),
      shape = 21, fill = ds_col, color = ds_col, size = 1.9, stroke = 0.2
    ) +
    facet_wrap(~ canonical_direction, nrow = 2, scales = "fixed",
               labeller = labeller(canonical_direction = function(x) paste("CD", x))) +
    angle_scale() +
    legend_layer() +
    labs(
      x = "taxonomic-level pair",
      y = "angle (rad)"
    ) +
    theme_si() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
      axis.text.y = element_text(size = 10),
      panel.grid.major.x = element_blank()
    )
}

# ---------------------------------------------------------------------------
run_dataset <- function(root, ds, perc, out_dir, n_cd = Inf) {
  dat <- load_dataset_inputs(root, ds, perc)

  real_angles <- compute_angles(dat$real_raw, has_seed = FALSE)
  null_angles <- NULL
  if (!is.null(dat$null_raw) && nrow(dat$null_raw) > 0) {
    null_angles <- compute_angles(dat$null_raw, has_seed = TRUE)
  }

  if (is.finite(n_cd)) {
    message("  Keeping canonical directions 1-", n_cd, " for ", ds)
    real_angles <- real_angles %>% filter(canonical_direction <= n_cd)
    if (!is.null(null_angles)) {
      null_angles <- null_angles %>% filter(canonical_direction <= n_cd)
    }
  }

  res <- summarize_angles(real_angles, null_angles, dat$pair_key)

  message("\n=== ", ds, " overall summary ===")
  print(res$real_overall)
  if (!is.null(res$null_summary)) print(res$null_summary)
  if (!is.null(res$true_vs_null)) print(res$true_vs_null)

  list(
    data = dat,
    real_angles = real_angles,
    null_angles = null_angles,
    res = res,
    overall_plot = make_overall_panel(res, ds),
    pair_plot = make_pair_panel(res, ds, n_cd = ifelse(is.finite(n_cd), n_cd, max(res$real_overall$canonical_direction)))
  )
}

# ---------------------------------------------------------------------------
main <- function() {
  args <- parse_cli_args()

  root <- if (is.na(args$root) || !nzchar(trimws(args$root))) {
    get_repo_root()
  } else {
    normalizePath(args$root, winslash = "/", mustWork = TRUE)
  }

  perc_short <- paste0("p", sub("^[^.]*\\.?", "", args$perc))
  message("perc: ", args$perc)

  out_dir <- file.path(root, "manuscript", "SI")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  soil_run  <- run_dataset(root, "soil",  args$perc, out_dir, n_cd = args$n_cd_soil)
  ocean_run <- run_dataset(root, "ocean", args$perc, out_dir, n_cd = args$n_cd_ocean)

  fig_violin <- (soil_run$overall_plot | ocean_run$overall_plot) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  fig_pair <- (soil_run$pair_plot / ocean_run$pair_plot) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  pdf_violin <- file.path(out_dir,
    sprintf("SI_crosstax_stability_mean_%s.pdf", perc_short))
  pdf_pairs  <- file.path(out_dir,
    sprintf("SI_crosstax_stability_pairs_%s.pdf",  perc_short))

  ggsave(pdf_violin, plot = fig_violin, width = 6.5, height = 3.0, units = "in", device = "pdf")
  ggsave(pdf_pairs,  plot = fig_pair,   width = 6.5, height = 7.0, units = "in", device = "pdf")
  message("Saved:\n  ", pdf_violin, "\n  ", pdf_pairs)
}

main()