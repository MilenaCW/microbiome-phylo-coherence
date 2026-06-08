# SI_PCA_comparison.R — Plot angle and EVR SI figures using pre-computed results from 08_PCA_compare.R.
#
# Reads from: <dataset>/results/CCA/<perc>/<tax_level>/08_PCA_compare/
#   angle_results.csv, null_results.csv, evr_results.csv
#
# Output PDFs (6.5 x 3.5 in):
#   manuscript/<dataset>_SI_angle_PCA_<perc_short>.pdf
#   manuscript/<dataset>_SI_EVR_PCA_<perc_short>.pdf
#
# Usage (from repo root):
#   Rscript code/manuscript_plotting/SI_PCA_comparison.R \
#     --dataset soil --n_cds 4 --perc_identity 0.90

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(optparse)
  library(patchwork)
  library(grid)
})

stopf <- function(...) stop(sprintf(...), call. = FALSE)

TAX_LEVELS  <- c("OTU", "Species", "Genus", "Family", "Order", "Class", "Phylum")
ICON_LABEL  <- " "   # placeholder facet label — renders as blank strip

TAX_PALETTE <- c(
  OTU     = "#812727",
  Species = "#a05532",
  Genus   = "#bd863e",
  Family  = "#dab449",
  Order   = "#dcbd6c",
  Class   = "#ddc68e",
  Phylum  = "#decfb2"
)

# ---------------------------------------------------------------------------
get_repo_root <- function() {
  args     <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) != 1) stop("Run with Rscript so --file= is available.", call. = FALSE)
  script_path <- normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)
  setup_path  <- normalizePath(file.path(dirname(script_path), "..", "setup.R"),
                               winslash = "/", mustWork = TRUE)
  source(setup_path, local = TRUE)
  get("REPO_ROOT", envir = .GlobalEnv)
}

parse_cli_args <- function(argv = commandArgs(trailingOnly = TRUE)) {
  option_list <- list(
    make_option("--dataset",       type = "character", default = NA,
                help = "Dataset: soil or ocean"),
    make_option("--perc_identity", type = "character", default = "0.90",
                help = "Percent identity [default %default]"),
    make_option("--root",          type = "character", default = NA,
                help = "Repo root path (inferred from script location if omitted)")
  )
  a <- parse_args(OptionParser(option_list = option_list), args = argv)
  if (is.na(a$dataset) || !nzchar(trimws(a$dataset)))
    stop("--dataset is required", call. = FALSE)
  list(dataset = trimws(tolower(a$dataset)),
       perc    = a$perc_identity,
       root    = a$root)
}

# ---------------------------------------------------------------------------
read_results <- function(root, dataset, perc) {
  base_path  <- file.path(root, dataset, "results", "CCA", perc)
  angle_list <- list(); null_list <- list(); evr_list <- list()

  for (tax in TAX_LEVELS) {
    out_dir <- file.path(base_path, tax, "step8_PCA_compare")
    af <- file.path(out_dir, "angle_results.csv")
    nf <- file.path(out_dir, "null_results.csv")
    ef <- file.path(out_dir, "evr_results.csv")
    if (!all(file.exists(af, nf, ef))) {
      message("  Skipping ", tax, " (missing 08_PCA_compare results)")
      next
    }
    angle_list[[tax]] <- read_csv(af, show_col_types = FALSE) %>% mutate(tax_level = tax)
    null_list[[tax]]  <- read_csv(nf, show_col_types = FALSE) %>% mutate(tax_level = tax)
    evr_list[[tax]]   <- read_csv(ef, show_col_types = FALSE) %>% mutate(tax_level = tax)
  }

  if (length(angle_list) == 0) stopf("No 08_PCA_compare results found under %s", base_path)

  list(
    df_angle = bind_rows(angle_list),
    df_null  = bind_rows(null_list),
    df_evr   = bind_rows(evr_list)
  )
}

# ---------------------------------------------------------------------------
# Deterministic per-(cd, fold) jitter: true angles and null pointranges overlay
# at the same x-position for each fold (same convention as plot_coherence.R).
make_jitter_lookup <- function(cd_vals, fold_vals, seed = 1L, width = 0.2) {
  keys <- expand.grid(cd = sort(unique(cd_vals)), fold = sort(unique(fold_vals)),
                      stringsAsFactors = FALSE)
  set.seed(seed)
  keys$jitter <- runif(nrow(keys), -width, width)
  keys
}

# ---------------------------------------------------------------------------
# Post-process the ggplotGrob: blank the first panel (top-left) and move the
# y-axis from the icon column to the next data column.
# Copied from SI_crosstax.R — works for ncol=4 with icon as first facet level.
fix_icon_panel <- function(g) {
  layout <- g$layout

  panel_layout <- layout[grepl("^panel-", layout$name), ]
  min_t        <- min(panel_layout$t)
  top_panels   <- panel_layout[panel_layout$t == min_t, ]
  icon_l       <- min(top_panels$l)
  otu_l        <- sort(top_panels$l)[2]

  icon_axis_col <- icon_l - 1
  otu_axis_col  <- otu_l  - 1

  axis_icon_idx <- which(layout$l == icon_axis_col & layout$r == icon_axis_col &
                           layout$t == min_t & layout$b == min_t)
  if (length(axis_icon_idx) > 0) {
    axis_grob <- g$grobs[[axis_icon_idx[1]]]
    g$grobs[[axis_icon_idx[1]]] <- zeroGrob()
    g$widths[otu_axis_col] <- g$widths[icon_axis_col]
    otu_axis_idx <- which(layout$l == otu_axis_col & layout$r == otu_axis_col &
                            layout$t == min_t & layout$b == min_t)
    if (length(otu_axis_idx) > 0) g$grobs[[otu_axis_idx[1]]] <- axis_grob
  }

  icon_panel_idx <- which(layout$l == icon_l & layout$r == icon_l &
                             layout$t == min_t & layout$b == min_t &
                             grepl("^panel-", layout$name))
  if (length(icon_panel_idx) > 0) g$grobs[[icon_panel_idx[1]]] <- zeroGrob()

  g
}

# Inject a legend grob into the first panel (top-left) and move the y-axis from
# the icon column to the next data column — same y-axis logic as fix_icon_panel,
# but replaces the panel with legend_grob instead of zeroGrob.
inject_first_panel_legend <- function(g, legend_grob) {
  layout <- g$layout

  panel_layout <- layout[grepl("^panel-", layout$name), ]
  min_t        <- min(panel_layout$t)
  top_panels   <- panel_layout[panel_layout$t == min_t, ]
  icon_l       <- min(top_panels$l)
  next_l       <- sort(top_panels$l)[2]

  icon_axis_col <- icon_l - 1
  next_axis_col <- next_l - 1

  axis_icon_idx <- which(layout$l == icon_axis_col & layout$r == icon_axis_col &
                           layout$t == min_t & layout$b == min_t)
  if (length(axis_icon_idx) > 0) {
    axis_grob <- g$grobs[[axis_icon_idx[1]]]
    g$grobs[[axis_icon_idx[1]]] <- zeroGrob()
    g$widths[next_axis_col] <- g$widths[icon_axis_col]
    next_axis_idx <- which(layout$l == next_axis_col & layout$r == next_axis_col &
                             layout$t == min_t & layout$b == min_t)
    if (length(next_axis_idx) > 0) g$grobs[[next_axis_idx[1]]] <- axis_grob
  }

  icon_panel_idx <- which(layout$l == icon_l & layout$r == icon_l &
                             layout$t == min_t & layout$b == min_t &
                             grepl("^panel-", layout$name))
  if (length(icon_panel_idx) > 0) g$grobs[[icon_panel_idx[1]]] <- legend_grob

  g
}

# ---------------------------------------------------------------------------
# Build a manual legend grob for the EVR plot.
# Rows: shape section (circle=train, square=test) then color section (grey=PCA, pie=CCA).
# All radial distances use "snpc" (smaller of panel width/height) so symbols
# stay circular regardless of panel aspect ratio.
build_evr_legend_grob <- function(tax_palette) {
  sym_x <- 0.40
  txt_x <- 0.53
  r     <- 0.04   # snpc units
  fsz   <- 10
  lwd   <- 0.7

  ys <- c(train = 0.82, test = 0.66, pca = 0.34, cca = 0.17)

  hdr_shape <- textGrob("data subset", x = 0.25, y = 1.00,
                         just = c("left", "top"),
                         gp = gpar(fontsize = fsz, col = "black", fontface = "italic"))
  hdr_color <- textGrob("technique", x = 0.25, y = 0.52,
                         just = c("left", "top"),
                         gp = gpar(fontsize = fsz, col = "black", fontface = "italic"))

  sym_train <- circleGrob(x = sym_x, y = ys["train"], r = unit(r, "snpc"),
                           gp = gpar(fill = "grey80", col = "black", lwd = lwd))
  sym_test  <- rectGrob(x = sym_x, y = ys["test"],
                         width = unit(r * 2, "snpc"), height = unit(r * 2, "snpc"),
                         gp = gpar(fill = "grey80", col = "black", lwd = lwd))
  sym_pca   <- circleGrob(x = sym_x, y = ys["pca"], r = unit(r, "snpc"),
                           gp = gpar(fill = "grey60", col = "black", lwd = lwd))

  # CCA: pie circle. Radial offsets use "snpc" so wedges match the circleGrob outline.
  cols     <- unname(tax_palette)
  n_pie    <- length(cols)
  thetas   <- seq(pi / 2, pi / 2 + 2 * pi, length.out = n_pie + 1)
  cy       <- unname(ys["cca"])
  cx_u     <- unit(sym_x, "npc")
  cy_u     <- unit(cy,    "npc")
  pie_segs <- lapply(seq_len(n_pie), function(i) {
    th <- seq(thetas[i], thetas[i + 1], length.out = 15)
    polygonGrob(
      x = unit.c(cx_u, cx_u + unit(r * cos(th), "snpc")),
      y = unit.c(cy_u, cy_u + unit(r * sin(th), "snpc")),
      gp = gpar(fill = cols[i], col = NA)
    )
  })
  sym_cca <- grobTree(
    do.call(grobTree, pie_segs),
    circleGrob(x = sym_x, y = cy, r = unit(r, "snpc"),
               gp = gpar(fill = NA, col = "black", lwd = lwd))
  )

  lbl_train <- textGrob("train", x = txt_x, y = ys["train"],
                         just = c("left", "center"), gp = gpar(fontsize = fsz, col = "black"))
  lbl_test  <- textGrob("test",  x = txt_x, y = ys["test"],
                         just = c("left", "center"), gp = gpar(fontsize = fsz, col = "black"))
  lbl_pca   <- textGrob("PCA",   x = txt_x, y = ys["pca"],
                         just = c("left", "center"), gp = gpar(fontsize = fsz, col = "black"))
  lbl_cca   <- textGrob("CCA",   x = txt_x, y = cy,
                         just = c("left", "center"), gp = gpar(fontsize = fsz, col = "black"))

  grobTree(hdr_shape, hdr_color,
           sym_train, sym_test, sym_pca, sym_cca,
           lbl_train, lbl_test, lbl_pca, lbl_cca)
}

# ---------------------------------------------------------------------------
make_angle_plot <- function(df_angle, df_null) {
  present   <- intersect(TAX_LEVELS, unique(df_angle$tax_level))
  fac_levels <- c(ICON_LABEL, present)   # icon first -> empty top-left panel

  df_angle <- df_angle %>% mutate(tax_level = factor(tax_level, levels = fac_levels))
  df_null  <- df_null  %>% mutate(tax_level = factor(tax_level, levels = fac_levels))

  jitter_df <- make_jitter_lookup(df_angle$cd, df_angle$fold)
  df_angle  <- df_angle %>% left_join(jitter_df, by = c("cd", "fold")) %>%
    mutate(x_pos = cd + jitter)
  df_null   <- df_null  %>% left_join(jitter_df, by = c("cd", "fold")) %>%
    mutate(x_pos = cd + jitter)

  # palette includes transparent for the icon facet so scale_color_manual is happy
  pal <- c(setNames("transparent", ICON_LABEL), TAX_PALETTE[present])

  cd_breaks <- sort(unique(df_angle$cd))
  y_breaks  <- c(0, pi / 4, pi / 2)
  y_labels  <- expression(0, pi/4, pi/2)

  xmin <- min(cd_breaks)
  xmax <- max(cd_breaks)
  range_x <- xmax - xmin

  ymin <- min(pi/4, min(df_angle$angle))
  ymax <- pi/2
  range_y <- ymax - ymin

  p <- ggplot() +
    geom_pointrange(
      data = df_null,
      aes(x = x_pos, y = null_mean,
          ymin = null_mean - null_sd, ymax = null_mean + null_sd),
      color = "grey60", alpha = 0.5, size = 0.2, fatten = 1.5,
      show.legend = FALSE
    ) +
    geom_point(
      data = df_angle,
      aes(x = x_pos, y = angle, color = tax_level),
      alpha = 0.7, size = 1.2,
      show.legend = FALSE
    ) +
    geom_hline(yintercept = pi / 2, linetype = "dashed", color = "grey40",
               linewidth = 0.5) +
    scale_color_manual(values = pal, drop = FALSE) +
    scale_x_continuous(breaks = cd_breaks, 
                       limits = c(xmin - range_x * 0.1, xmax + range_x * 0.1)) +
    scale_y_continuous(breaks = y_breaks, labels = y_labels,
                       limits = c(ymin - range_y * 0.05, ymax + range_y * 0.05)) +
    facet_wrap(vars(tax_level), ncol = 4, drop = FALSE) +
    labs(x = "canonical direction", y = "angle with PC1 (rad)") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "none", 
          strip.text = element_text(size = 10),
          panel.spacing.x = unit(0.25, "in"))

  grDevices::pdf(NULL)
  g <- fix_icon_panel(ggplotGrob(p))
  grDevices::dev.off()
  wrap_elements(full = g)
}

# ---------------------------------------------------------------------------
make_evr_plot <- function(df_evr) {
  present    <- intersect(TAX_LEVELS, unique(df_evr$tax_level))
  fac_levels <- c(ICON_LABEL, present)   # icon first -> empty top-left panel for legend

  df_evr <- df_evr %>%
    mutate(
      tax_level   = factor(tax_level, levels = fac_levels),
      color_group = factor(
        if_else(method == "CCA", as.character(tax_level), "PCA"),
        levels = c(present, "PCA")
      ),
      split = factor(split, levels = c("train", "test"))
    )

  pal_evr <- c(setNames("transparent", ICON_LABEL), TAX_PALETTE[present], PCA = "grey60")

  df_summ <- df_evr %>%
    group_by(tax_level, color_group, split, j) %>%
    summarise(mean_evr = mean(evr, na.rm = TRUE),
              sd_evr   = sd(evr,   na.rm = TRUE),
              .groups  = "drop")

  cd_breaks <- sort(unique(df_evr$j))
  xmin <- min(cd_breaks)
  xmax <- max(cd_breaks)
  range_x <- xmax - xmin

  p <- ggplot() +
    geom_point(
      data = df_evr,
      aes(x = j, y = evr, fill = color_group, shape = split),
      color = "black", stroke = 0.2, alpha = 0.25, size = 2.0,
      position = position_jitter(width = 0.08, seed = 2L),
      show.legend = FALSE
    ) +
    # geom_errorbar(
    #   data = df_summ,
    #   aes(x = j, ymin = mean_evr - sd_evr, ymax = mean_evr + sd_evr,
    #       color = color_group),
    #   alpha = 0.75, linewidth = 0.25, width = 0.5,
    #   show.legend = FALSE
    # ) +
    geom_point(
      data = df_summ,
      aes(x = j, y = mean_evr, fill = color_group, shape = split),
      color = "black", stroke = 0.2, alpha = 1.0, size = 2.5,
      show.legend = FALSE
    ) +
    scale_fill_manual(values  = pal_evr, drop = FALSE) +
    scale_color_manual(values = pal_evr, drop = FALSE) +
    scale_shape_manual(values = c(train = 21L, test = 22L)) +
    scale_x_continuous(breaks = cd_breaks, 
                       limits = c(xmin - range_x * 0.1, xmax + range_x * 0.1)) +
    scale_y_continuous(limits = c(0, NA)) +
    facet_wrap(vars(tax_level), ncol = 4, drop = FALSE) +
    labs(x = "dimension of latent space", y = "cumulative explained variance ratio") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "none", 
          strip.text = element_text(size = 10),
          panel.spacing.x = unit(0.25, "in"))

  grDevices::pdf(NULL)
  g <- ggplotGrob(p)
  grDevices::dev.off()
  g <- inject_first_panel_legend(g, build_evr_legend_grob(TAX_PALETTE[present]))
  wrap_elements(full = g)
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
  message("Dataset: ", args$dataset, " | perc: ", args$perc)

  res <- read_results(root, args$dataset, args$perc)

  p_angle <- make_angle_plot(res$df_angle, res$df_null)
  p_evr   <- make_evr_plot(res$df_evr)

  out_dir <- file.path(root, "manuscript", "SI")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  out_angle <- file.path(out_dir,
    sprintf("SI_%s_angle_PCA_%s.pdf", args$dataset, perc_short))
  out_evr <- file.path(out_dir,
    sprintf("SI_%s_EVR_PCA_%s.pdf",   args$dataset, perc_short))

  ggsave(out_angle, plot = p_angle, width = 6.5, height = 3.5, units = "in")
  ggsave(out_evr,   plot = p_evr,   width = 6.5, height = 3.5, units = "in")
  message("Saved:\n  ", out_angle, "\n  ", out_evr)
}

main()
