# SI_GG2_performance.R
# ---------------------
# Produces a two-panel SI figure showing GG2 database performance across
# percent sequence identity thresholds for soil and ocean datasets:
#   Panel A: percent reads matched vs perc identity
#   Panel B: occupancy (taxa per GG2 feature) vs perc identity, log10 scale
#
# Usage (from repo root):
#   Rscript code/manuscript_plotting/SI_GG2_performance.R
#
# Output: manuscript/SI_GG2_performance.pdf  (6 x 3 in)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(patchwork)
})

# =============================================================================
# USER-ADJUSTABLE PARAMETERS
# =============================================================================

soil_selected_pi  <- 0.90   # percent identity used in soil analysis
ocean_selected_pi <- 0.90   # percent identity used in ocean analysis

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

root <- get_repo_root()

soil_tsv  <- file.path(root, "soil/data/processed_data/16S/GG2/summary_diagnostics/gg2_sweep_summary.tsv")
ocean_tsv <- file.path(root, "ocean/data/processed_data/16S/GG2/summary_diagnostics/gg2_sweep_summary.tsv")
out_path  <- file.path(root, "manuscript/SI/SI_GG2_performance.pdf")

# =============================================================================
# COLORS  (consistent with all other manuscript plots)
# =============================================================================

dataset_colors <- c(
  "soil"  = "#8f723d",
  "ocean" = "#90afa7"
)

# =============================================================================
# DATA
# =============================================================================

selected_pi <- c("soil" = soil_selected_pi, "ocean" = ocean_selected_pi)

all_data <- bind_rows(
  read_tsv(soil_tsv,  show_col_types = FALSE) %>% mutate(dataset = "soil"),
  read_tsv(ocean_tsv, show_col_types = FALSE) %>% mutate(dataset = "ocean")
) %>%
  mutate(
    perc_identity_pct = perc_identity * 100,
    selected = perc_identity == selected_pi[dataset]
  )

highlight_df <- filter(all_data, selected)

# =============================================================================
# PANEL A: Read mapping success
# =============================================================================

p_reads <- ggplot(all_data,
                  aes(x = perc_identity_pct, y = gg2_mapped_percent,
                      color = dataset, group = dataset)) +
  geom_line(linewidth = 0.75) +
  geom_point(shape = 19, size = 1) +
  geom_point(data = highlight_df,
             shape = 8, size = 2, stroke = 1.0) +
  scale_color_manual(values = dataset_colors) +
  scale_x_continuous(limits = c(80, 100)) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(
    x = "percent sequence identity (%)",
    y = "percent reads matched (%)",
    color = NULL,
    title = "(a)"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor       = element_blank(),
    plot.title             = element_text(hjust = 0.5, size = 10),
    plot.title.position    = "panel",
    axis.text              = element_text(size = 10),
    axis.title             = element_text(size = 10),
    plot.margin            = margin(2, 4, 2, 2)
  )

# =============================================================================
# PANEL B: Occupancy
# =============================================================================

p_occ <- ggplot(all_data,
                aes(x = perc_identity_pct, color = dataset, group = dataset)) +
  geom_line(aes(y = min_occupancy), linewidth = 0.75, linetype = "dashed") +
  geom_line(aes(y = max_occupancy), linewidth = 0.75, linetype = "dashed") +
  geom_line(aes(y = mean_occupancy), linewidth = 0.75) +
  geom_point(data = highlight_df,
             aes(y = mean_occupancy),
             shape = 8, size = 2, stroke = 1.0) +
  scale_color_manual(values = dataset_colors) +
  scale_x_continuous(limits = c(80, 100)) +
  scale_y_log10(
    breaks       = 10^(0:4),
    minor_breaks = c(outer(2:9, 10^(0:3)))
  ) +
  labs(
    x = "percent sequence identity (%)",
    y = "occupancy (no. taxa / feature)",
    color = NULL,
    title = "(b)"
  ) +
  guides(color = "none") +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor.x     = element_blank(),
    plot.title             = element_text(hjust = 0.5, size = 10),
    plot.title.position    = "panel",
    axis.text              = element_text(size = 10),
    axis.title             = element_text(size = 10),
    plot.margin            = margin(2, 4, 2, 2)
  )

# =============================================================================
# COMBINE AND SAVE
# =============================================================================

p_combined <- p_reads + p_occ +
  plot_layout(guides = "collect") &
  theme(
    legend.position   = "bottom",
    legend.margin     = margin(0, 0, 0, 0),
    legend.box.margin = margin(-6, 0, 0, 0)
  )

ggsave(out_path, plot = p_combined, width = 6, height = 3, dpi = 300)
cat("Saved:", out_path, "\n")
