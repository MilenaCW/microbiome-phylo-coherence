# SI_taxonomic_assignment.R
# --------------------------
# Produces a single-panel SI figure showing the percent of OTUs/ASVs assigned
# taxonomy at each rank (species -> phylum) for soil and ocean.
#
# Motivates dropping species-level analysis from the pipeline: a large
# fraction of features have no species-level assignment.
#
# Usage (from repo root):
#   Rscript code/manuscript_plotting/SI_taxonomic_assignment.R
#   Rscript code/manuscript_plotting/SI_taxonomic_assignment.R --perc_identity 0.90
#
# Output: manuscript/SI/SI_taxonomic_assignment.pdf

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(optparse)
})

stopf <- function(...) stop(sprintf(...), call. = FALSE)

# Finest to coarsest
TAX_LEVELS_ORDER <- c("Species", "Genus", "Family", "Order", "Class", "Phylum")

# =============================================================================
# COLORS (consistent with all other manuscript plots)
# =============================================================================

dataset_colors <- c(
  "soil"  = "#8f723d",
  "ocean" = "#90afa7"
)

# =============================================================================
# PATHS
# =============================================================================

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

parse_cli_args <- function(argv = commandArgs(trailingOnly = TRUE)) {
  option_list <- list(
    make_option(
      "--perc_identity",
      type    = "character",
      default = "0.90",
      help    = "Percent identity used for GG2 taxonomy (e.g. 0.90) [default %default]"
    )
  )
  parse_args(OptionParser(option_list = option_list), args = argv)
}

# =============================================================================
# DATA
# =============================================================================

compute_percent_assigned <- function(taxonomy, levels) {
  tibble(tax_level = levels) %>%
    rowwise() %>%
    mutate(
      percent_assigned = 100 * mean(taxonomy[[tax_level]] != "Unassigned", na.rm = TRUE)
    ) %>%
    ungroup()
}

main <- function() {
  args <- parse_cli_args()
  root <- get_repo_root()

  soil_tax_fn  <- file.path(root, "soil/data/processed_data/16S/GG2",  args$perc_identity, "final", "taxonomy.csv")
  ocean_tax_fn <- file.path(root, "ocean/data/processed_data/16S/GG2", args$perc_identity, "final", "taxonomy.csv")

  for (f in c(soil_tax_fn, ocean_tax_fn)) {
    if (!file.exists(f)) stopf("taxonomy.csv not found: %s", f)
  }

  soil_tax  <- read_csv(soil_tax_fn,  show_col_types = FALSE)
  ocean_tax <- read_csv(ocean_tax_fn, show_col_types = FALSE)

  all_data <- bind_rows(
    compute_percent_assigned(soil_tax,  TAX_LEVELS_ORDER) %>% mutate(dataset = "soil"),
    compute_percent_assigned(ocean_tax, TAX_LEVELS_ORDER) %>% mutate(dataset = "ocean")
  ) %>%
    mutate(tax_level = factor(tax_level, levels = TAX_LEVELS_ORDER))

  # ===========================================================================
  # PLOT
  # ===========================================================================

  p <- ggplot(all_data,
              aes(x = tax_level, y = percent_assigned, color = dataset, group = dataset)) +
    geom_line(linewidth = 0.75) +
    geom_point(shape = 19, size = 1.5) +
    scale_color_manual(values = dataset_colors) +
    scale_y_continuous(limits = c(0, 100)) +
    labs(
      x = "taxonomic level",
      y = "percent OTUs assigned (%)",
      color = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor    = element_blank(),
      axis.text           = element_text(size = 10),
      axis.title          = element_text(size = 10),
      legend.position     = "bottom",
      legend.margin       = margin(0, 0, 0, 0),
      legend.box.margin   = margin(-6, 0, 0, 0),
      plot.margin         = margin(2, 4, 2, 2)
    )

  out_dir <- file.path(root, "manuscript", "SI")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_path <- file.path(out_dir, "SI_taxonomic_assignment.pdf")

  ggsave(out_path, plot = p, width = 4, height = 3, dpi = 300)
  cat("Saved:", out_path, "\n")
}

main()
