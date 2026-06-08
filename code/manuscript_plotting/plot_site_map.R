# Plot site locations from one or more datasets on a world map.
#
# Usage (from repo root):
#   Rscript code/manuscript_plotting/plot_site_map.R --soil --ocean
#   Rscript code/manuscript_plotting/plot_site_map.R --wwtp --full
#
# Output: manuscript/site_map.pdf

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
})

stopf <- function(...) stop(sprintf(...), call. = FALSE)

# Dataset colors and shapes (edit here to change appearance)
dataset_colors <- c(
  "soil"  = "#8f723d",
  "ocean" = "#90afa7",
  "wwtp"  = "#606161"
)
dataset_shapes <- c(
  "soil"  = 16L,
  "ocean" = 16L,
  "wwtp"  = 16L
)

parse_cli_args <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (!requireNamespace("optparse", quietly = TRUE)) {
    stopf("R package 'optparse' is required. Install it with: install.packages('optparse')")
  }

  option_list <- list(
    optparse::make_option(
      "--soil",
      action = "store_true",
      default = FALSE,
      help = "Include Soil sites"
    ),
    optparse::make_option(
      "--ocean",
      action = "store_true",
      default = FALSE,
      help = "Include Ocean sites"
    ),
    optparse::make_option(
      "--wwtp",
      action = "store_true",
      default = FALSE,
      help = "Include WWTP sites"
    ),
    optparse::make_option(
      "--full",
      action = "store_true",
      default = FALSE,
      help = "Use 'full' environmental data; if not set, uses 'filtered'"
    ),
    optparse::make_option(
      "--root",
      type = "character",
      default = NA_character_,
      help = "Repo root path (optional; if omitted, uses code/setup.R)"
    )
  )

  p <- optparse::OptionParser(option_list = option_list)
  a <- optparse::parse_args(p, args = argv, positional_arguments = FALSE)

  list(
    soil = a$soil,
    ocean = a$ocean,
    wwtp = a$wwtp,
    full = a$full,
    root = a$root
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

load_site_locations <- function(root, subdir, selected) {
  out_list <- list()
  for (flag in names(selected)) {
    if (!selected[[flag]]) next
    path <- file.path(root, flag, "data", "processed_data", "environmental", subdir, "site_locations.csv")
    if (!file.exists(path)) {
      stopf("Site locations file not found: %s\nGenerate it first (e.g. run read_envdata.R for this dataset and output dir).", path)
    }
    d <- readr::read_csv(path, show_col_types = FALSE)
    d$dataset <- flag
    d$sample_id <- as.character(d$sample_id)
    d$latitude <- as.numeric(d$latitude)
    d$longitude <- as.numeric(d$longitude)
    out_list[[flag]] <- d
  }
  dplyr::bind_rows(out_list)
}

check_no_missing_coords <- function(dat) {
  missing <- is.na(dat$latitude) | is.na(dat$longitude)
  n_missing <- sum(missing, na.rm = TRUE)
  if (n_missing > 0) {
    stopf("%d rows have missing latitude or longitude. Fix or investigate the data before plotting.", n_missing)
  }
  invisible(dat)
}

main <- function() {
  args <- parse_cli_args()
  selected <- list(
    soil = args$soil,
    ocean = args$ocean,
    wwtp = args$wwtp
  )
  if (!any(unlist(selected))) {
    stop("At least one dataset must be selected (--soil, --ocean, and/or --wwtp).", call. = FALSE)
  }

  root <- args$root
  if (is.na(root) || !nzchar(trimws(root))) {
    root <- get_repo_root()
  } else {
    root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  }

  subdir <- if (args$full) "full" else "filtered"
  dat <- load_site_locations(root, subdir, selected)
  check_no_missing_coords(dat)

  # Restrict colors/shapes to datasets present
  labs <- unique(dat$dataset)
  col_scale <- dataset_colors[names(dataset_colors) %in% labs]
  shape_scale <- dataset_shapes[names(dataset_shapes) %in% labs]

  world <- ggplot2::map_data("world") %>%
    dplyr::filter(.data$region != "Antarctica")
  p <- ggplot() +
    geom_polygon(data = world, aes(.data$long, .data$lat, group = .data$group), fill = "gray90", color = "gray70", linewidth = 0.2) +
    # 'size' in geom_point() is in mm, but you can convert points (pt) to mm:
    # 1 pt = 0.3527778 mm
    # So, e.g., for 10 pt: size = 10 * 0.3527778
    geom_point(
      data = dat,
      aes(.data$longitude, .data$latitude, color = .data$dataset, shape = .data$dataset),
      size = 5 * 0.3527778,
      alpha = 0.5
    ) +
    scale_color_manual(values = col_scale, breaks = names(col_scale)) +
    scale_shape_manual(values = shape_scale, breaks = names(shape_scale)) +
    # guides(
    #   color = guide_legend(override.aes = list(size = 2)),
    #   shape  = guide_legend(override.aes = list(size = 2))
    # ) +
    guides(
      color = "none",
      shape = "none"
    ) +
    coord_quickmap() +
    theme(
      panel.background = element_rect(fill = "white"),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.margin = margin(0, 0, 0, 0),
      text = element_text(size = 10)
    )

  out_dir <- file.path(root, "manuscript")
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  out_path <- file.path(out_dir, "site_map.pdf")
  ggsave(out_path, plot = p, width = 6, height = 2.5)
  message("Saved: ", out_path)
}

main()
