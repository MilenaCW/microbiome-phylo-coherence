# SI_crosstax_corr.R — SI figure: cross-taxonomic CCA correlations for soil (top) and ocean (bottom).
# Colored by taxonomic level; first facet left empty as icon placeholder.
#
# Usage (from repo root):
#   Rscript code/manuscript_plotting/SI_crosstax_corr.R
#   Rscript code/manuscript_plotting/SI_crosstax_corr.R --perc_identity 0.90
#
# Output: manuscript/SI_crosstax.pdf

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(patchwork)
  library(optparse)
  library(grid)
})

stopf <- function(...) stop(sprintf(...), call. = FALSE)

# Canonical order of taxonomic levels (coarsest last = lightest colour)
TAX_LEVELS_ORDER <- c("OTU", "Species", "Genus", "Family", "Order", "Class", "Phylum")

# Manuscript colour palette (matches plot_crosstax_env_loadings.R)
TAX_PALETTE <- c(
  OTU     = "#812727",
  Species = "#a05532",
  Genus   = "#bd863e",
  Family  = "#dab449",
  Order   = "#dcbd6c",
  Class   = "#ddc68e",
  Phylum  = "#decfb2"
)

# Label used for the empty icon-placeholder facet
ICON_LABEL <- " "

# ---------------------------------------------------------------------------
get_repo_root <- function() {
  args     <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) != 1) {
    stop("This script must be run with Rscript so --file= is available.", call. = FALSE)
  }
  script_path <- normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)
  script_dir  <- dirname(script_path)
  setup_path  <- normalizePath(file.path(script_dir, "..", "setup.R"), winslash = "/", mustWork = TRUE)
  source(setup_path, local = TRUE)
  get("REPO_ROOT", envir = .GlobalEnv)
}

# ---------------------------------------------------------------------------
parse_cli_args <- function(argv = commandArgs(trailingOnly = TRUE)) {
  option_list <- list(
    optparse::make_option(
      "--perc_identity",
      type    = "character",
      default = "0.90",
      help    = "Percent identity used for CCA results (e.g. 0.90) [default %default]"
    ),
    optparse::make_option(
      "--root",
      type    = "character",
      default = NA_character_,
      help    = "Repo root path (optional; inferred from script location if omitted)"
    )
  )
  p <- optparse::OptionParser(option_list = option_list)
  a <- optparse::parse_args(p, args = argv, positional_arguments = FALSE)
  list(perc_identity = a$perc_identity, root = a$root)
}

# ---------------------------------------------------------------------------
load_dataset <- function(root, ds, perc) {
  results_path <- file.path(root, ds, "results", "CCA")
  compare_base <- file.path(results_path, perc)

  if (!dir.exists(compare_base)) {
    stopf("Results directory not found: %s", compare_base)
  }

  # Discover tax-level directories that have correlations_per_fold.csv
  all_dirs <- list.dirs(compare_base, full.names = FALSE, recursive = FALSE)
  tax_dirs <- character(0)
  for (d in all_dirs) {
    if (d == "" || d == "crosstax_compare" || d == "Domain") next
    corr_file <- file.path(compare_base, d, "step2_loadings", "correlations_per_fold.csv")
    if (file.exists(corr_file)) tax_dirs <- c(tax_dirs, d)
  }
  if (length(tax_dirs) == 0) {
    stopf("No tax-level dirs with correlations_per_fold.csv found under: %s", compare_base)
  }

  # Order factor levels canonically
  plot_tax_levels <- unique(c(TAX_LEVELS_ORDER, tax_dirs))
  plot_tax_levels <- plot_tax_levels[plot_tax_levels %in% tax_dirs]

  # Collect test-set correlation summaries per tax_level x canonical_direction
  corr_list <- list()
  for (tax_level in tax_dirs) {
    corr_file <- file.path(compare_base, tax_level, "step2_loadings", "correlations_per_fold.csv")
    corr_raw  <- read_csv(corr_file, show_col_types = FALSE)
    corr_sum  <- corr_raw %>%
      filter(!is.na(test)) %>%
      group_by(canonical_direction) %>%
      summarise(mean = mean(test, na.rm = TRUE),
                sd   = sd(test,   na.rm = TRUE),
                .groups = "drop") %>%
      mutate(tax_level = tax_level)
    corr_list[[tax_level]] <- corr_sum
  }
  corr_data <- bind_rows(corr_list)

  # Collect null expectations (CD=1 per tax_level) if available
  null_list <- list()
  for (tax_level in tax_dirs) {
    null_file <- file.path(compare_base, tax_level, "step3_null", "null_correlations_per_fold.csv")
    if (!file.exists(null_file)) next
    null_raw  <- read_csv(null_file, show_col_types = FALSE)
    null_cd1  <- null_raw %>%
      filter(canonical_direction == 1, !is.na(test)) %>%
      group_by(seed) %>%
      summarise(mean_test = mean(test, na.rm = TRUE), .groups = "drop") %>%
      summarise(mean = mean(mean_test), sd = sd(mean_test), n = dplyr::n(), .groups = "drop") %>%
      mutate(tax_level = tax_level)
    null_list[[tax_level]] <- null_cd1
  }
  null_data <- if (length(null_list) > 0) bind_rows(null_list) else NULL

  n_cd <- max(corr_data$canonical_direction, na.rm = TRUE)

  list(
    corr_data      = corr_data,
    null_data      = null_data,
    n_cd           = n_cd,
    plot_tax_levels = plot_tax_levels
  )
}

# ---------------------------------------------------------------------------
# Post-process the icon panel in a ggplotGrob:
#   1. Blank the icon panel background / grid.
#   2. Move the row-1 y-axis from the icon column to the OTU column.
#
# The gtable structure (verified by inspection for ncol = 4, 8 facets):
#   col 6  (icon_axis_col)  — y-axis for col-1 panels (both rows)
#   col 7  (icon_panel_col) — panel col 1 (icon row 1, Family row 2)
#   col 9  (5pt spacer)
#   col 10 (otu_axis_col)   — y-axis slot for col-2 panels, 0-width by default
#   col 11 (otu_panel_col)  — panel col 2 (OTU row 1, Order row 2)
#
# Steps:
#   a. Save the row-1 y-axis grob from icon_axis_col.
#   b. Replace it with zeroGrob (row-2 Family y-axis in the same column is untouched).
#   c. Expand otu_axis_col to the same width and plant the saved grob there.
#   d. Blank the icon panel grob (removes background + grid lines).
fix_icon_panel <- function(g) {
  layout <- g$layout

  # Locate the icon panel: minimum t among data panels, minimum l in that row
  panel_layout <- layout[grepl("^panel-", layout$name), ]
  min_t        <- min(panel_layout$t)
  top_panels   <- panel_layout[panel_layout$t == min_t, ]
  icon_l       <- min(top_panels$l)   # e.g. col 7
  otu_l        <- sort(top_panels$l)[2]  # e.g. col 11

  icon_axis_col <- icon_l - 1   # e.g. col 6  — holds axis-l for col-1 panels
  otu_axis_col  <- otu_l  - 1   # e.g. col 10 — holds axis-l slot for col-2 panels

  # (a+b) Grab and clear the row-1 icon y-axis
  axis_icon_idx <- which(layout$l == icon_axis_col & layout$r == icon_axis_col &
                           layout$t == min_t & layout$b == min_t)
  if (length(axis_icon_idx) > 0) {
    axis_grob <- g$grobs[[axis_icon_idx[1]]]
    g$grobs[[axis_icon_idx[1]]] <- zeroGrob()

    # (c) Widen the OTU axis slot and plant the grob
    g$widths[otu_axis_col] <- g$widths[icon_axis_col]
    otu_axis_idx <- which(layout$l == otu_axis_col & layout$r == otu_axis_col &
                            layout$t == min_t & layout$b == min_t)
    if (length(otu_axis_idx) > 0) g$grobs[[otu_axis_idx[1]]] <- axis_grob
  }

  # (d) Blank the icon panel (panel grob at icon_l, min_t)
  icon_panel_idx <- which(layout$l == icon_l & layout$r == icon_l &
                             layout$t == min_t & layout$b == min_t &
                             grepl("^panel-", layout$name))
  if (length(icon_panel_idx) > 0) g$grobs[[icon_panel_idx[1]]] <- zeroGrob()

  g
}

# ---------------------------------------------------------------------------
make_plot <- function(corr_data, null_data, n_cd, plot_tax_levels, title) {
  # Factor levels: icon placeholder first, then tax levels in canonical order
  fac_levels <- c(ICON_LABEL, plot_tax_levels)

  # Build complete colour palette (icon = transparent, tax levels = manuscript colours)
  missing_cols <- plot_tax_levels[!plot_tax_levels %in% names(TAX_PALETTE)]
  if (length(missing_cols) > 0) {
    stopf("No colour specified for taxonomic level(s): %s. Add to TAX_PALETTE.",
          paste(missing_cols, collapse = ", "))
  }
  pal <- c(setNames("transparent", ICON_LABEL),
           setNames(TAX_PALETTE[plot_tax_levels], plot_tax_levels))

  # Prepare correlation data with facet factor
  plot_data <- corr_data %>%
    mutate(tax_level = factor(tax_level, levels = fac_levels))

  ymin_val <- min(plot_data$mean - plot_data$sd, na.rm = TRUE)
  ymin     <- min(0, if (is.finite(ymin_val)) ymin_val else 0)

  p <- ggplot(plot_data,
              aes(x = factor(canonical_direction), y = mean, color = tax_level)) +
    scale_color_manual(values = pal, drop = FALSE) +
    facet_wrap(vars(tax_level), ncol = 4, drop = FALSE) +
    ylim(ymin, NA) +
    ggtitle(title) +
    # Preserve axis-label spacing without showing text
    labs(x = " ", y = " ") +
    theme_minimal(base_size = 10) +
    theme(
      plot.title      = element_text(hjust = 0.5, size = 10),
      axis.title      = element_text(colour = "transparent", size = 10),
      legend.position = "none"
    )

  # Null distribution band + line (skip for icon facet — icon level has no null data)
  if (!is.null(null_data) && nrow(null_data) > 0) {
    null_plot <- null_data %>%
      mutate(tax_level = factor(tax_level, levels = fac_levels))
    p <- p +
      geom_rect(
        data        = null_plot,
        aes(xmin = 0.5, xmax = n_cd + 0.5,
            ymin = mean - sd / sqrt(n), ymax = mean + sd / sqrt(n)),
        fill        = "grey",
        alpha       = 0.15,
        inherit.aes = FALSE
      ) +
      geom_hline(
        data      = null_plot,
        aes(yintercept = mean),
        linetype  = "dashed",
        color     = "grey",
        linewidth = 0.7
      )
  }

  p <- p +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0, linewidth = 0.8)

  wrap_elements(full = fix_icon_panel(ggplotGrob(p)))
}

# ---------------------------------------------------------------------------
main <- function() {
  args <- parse_cli_args()
  perc <- args$perc_identity

  root <- args$root
  if (is.na(root) || !nzchar(trimws(root))) {
    root <- get_repo_root()
  } else {
    root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  }

  soil  <- load_dataset(root, "soil",  perc)
  ocean <- load_dataset(root, "ocean", perc)

  p_soil  <- make_plot(soil$corr_data,  soil$null_data,  soil$n_cd,  soil$plot_tax_levels,  "Soil")
  p_ocean <- make_plot(ocean$corr_data, ocean$null_data, ocean$n_cd, ocean$plot_tax_levels, "Ocean")

  combined <- p_soil / p_ocean

  out_dir <- file.path(root, "manuscript", "SI")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_path <- file.path(out_dir, "SI_crosstax.pdf")

  ggsave(out_path, plot = combined, width = 6.5, height = 7, units = "in")
  message("Saved: ", out_path)
}

main()
