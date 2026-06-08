#!/usr/bin/env Rscript
# SI_interintra_dist.R — SI figure: inter vs intra-group cophenetic distances by taxonomic level.
# One panel per taxonomic level; each panel overlays soil and ocean density curves
# with dataset colour and linetype (solid = intra, dashed = inter).
#
# Usage (from repo root):
#   Rscript code/manuscript_plotting/SI_interintra_dist.R
#   Rscript code/manuscript_plotting/SI_interintra_dist.R --perc_identity 0.90
#
# Output: manuscript/SI_interintra_dist.pdf

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(optparse)
})

stopf <- function(...) stop(sprintf(...), call. = FALSE)

# Finest to coarsest -- Species excluded
TAX_LEVELS_ORDER <- c("Genus", "Family", "Order", "Class", "Phylum")

# Strip background colours per level (matches other SI figures)
TAX_PALETTE <- c(
  Species = "#a05532",
  Genus   = "#bd863e",
  Family  = "#dab449",
  Order   = "#dcbd6c",
  Class   = "#ddc68e",
  Phylum  = "#decfb2"
)

DATASET_COLORS <- c(soil = "#8f723d", ocean = "#90afa7")

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

parse_cli_args <- function(argv = commandArgs(trailingOnly = TRUE)) {
  option_list <- list(
    optparse::make_option(
      "--perc_identity",
      type    = "character",
      default = "0.90",
      help    = "Percent identity used for CCA results [default %default]"
    ),
    optparse::make_option(
      "--output_path",
      type    = "character",
      default = NULL,
      help    = "Output directory [default: manuscript/ under repo root]"
    ),
    optparse::make_option(
      "--verbose",
      action  = "store_true",
      default = FALSE,
      help    = "Print progress messages"
    )
  )
  optparse::parse_args(optparse::OptionParser(option_list = option_list), args = argv)
}

# ---------------------------------------------------------------------------
REPO_ROOT <- get_repo_root()
args      <- parse_cli_args()

soil_dir  <- file.path(REPO_ROOT, "soil/results/reference_distances",  args$perc_identity)
ocean_dir <- file.path(REPO_ROOT, "ocean/results/reference_distances", args$perc_identity)

for (d in c(soil_dir, ocean_dir)) {
  if (!dir.exists(d)) stopf("Results directory not found: %s", d)
}

# Determine which levels are available for both datasets
present_levels <- Filter(function(lvl) {
  file.exists(file.path(soil_dir,  sprintf("dist_%s.rds", lvl))) &&
  file.exists(file.path(ocean_dir, sprintf("dist_%s.rds", lvl)))
}, TAX_LEVELS_ORDER)

if (length(present_levels) == 0) stopf("No per-level RDS files found in result directories.")
cat(sprintf("Levels found: %s\n", paste(present_levels, collapse = ", ")))

out_dir <- if (!is.null(args$output_path)) args$output_path else
  file.path(REPO_ROOT, "manuscript/SI")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
out_fn <- file.path(out_dir, "SI_interintra_dist.pdf")

panels <- vector("list", length(present_levels))
for (i in seq_along(present_levels)) {
  lvl <- present_levels[i]
  if (args$verbose) cat(sprintf("Building panel: %s\n", lvl))

  df <- dplyr::bind_rows(
    readRDS(file.path(soil_dir,  sprintf("dist_%s.rds", lvl))),
    readRDS(file.path(ocean_dir, sprintf("dist_%s.rds", lvl)))
  ) %>% mutate(
    type    = factor(type,    levels = c("intra", "inter")),
    dataset = factor(dataset, levels = c("soil", "ocean"))
  )

  panels[[i]] <- ggplot(df, aes(x = dist, colour = dataset, linetype = type)) +
    geom_density(key_glyph = "point", show.legend = c(colour = TRUE,  linetype = FALSE)) +
    geom_density(key_glyph = "path",  show.legend = c(colour = FALSE, linetype = TRUE)) +
    scale_colour_manual(values = DATASET_COLORS,
                        labels = c(soil = "soil", ocean = "ocean")) +
    scale_linetype_manual(values = c(intra = "solid", inter = "dotdash"),
                          labels = c(intra = "intra-group", inter = "inter-group")) +
    labs(x = "cophenetic distance", y = "density",
         title = lvl, colour = NULL, linetype = NULL) +
    guides(
      colour   = guide_legend(override.aes = list(shape = 16, size = 3)),
      linetype = guide_legend(override.aes = list(colour = "grey30"))
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title       = element_text(size = 10, hjust = 0.5, margin = margin(3, 2, 3, 2)),
      # plot.background  = element_rect(fill = TAX_PALETTE[lvl], colour = NA),
      # panel.background = element_rect(fill = "white", colour = NA)
    )

  rm(df); gc()
}

fig <- wrap_plots(c(panels, list(guide_area()))) +
  plot_layout(guides = "collect", ncol = 3) &
  theme(legend.text  = element_text(size = 10),
        legend.title = element_text(size = 10))

ggsave(out_fn, plot = fig, width = 6.5, height = 5)
cat(sprintf("Saved: %s\n", out_fn))
