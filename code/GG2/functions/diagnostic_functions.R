# Diagnostic Functions for GG2 generalized pipeline
# -------------------------------------------------
# Modular helpers for:
#   - Reading GG2 taxonomy
#   - Computing read distributions
#   - Basic occupancy summaries
#   - Plotting distributions and phylogenetic trees

suppressPackageStartupMessages({
  library(ape)
  library(ggplot2)
  library(ggtree)
  library(dplyr)
  library(readr)
  library(RColorBrewer)
  library(scales)
  library(tidyr)
})

########################################
## Data-reading helpers              ##
########################################

# Read cleaned GG2 taxonomy.csv (final)
# Expects first column 'Feature_ID' and additional taxonomic columns
read_clean_taxonomy <- function(path) {
  tax <- read_csv(path, show_col_types = FALSE)
  if (!"Feature_ID" %in% colnames(tax)) {
    stop("Expected 'Feature_ID' column in taxonomy.csv")
  }
  tax
}

########################################
## Calculations                      ##
########################################

# Calculate read distribution statistics (counts and fractions)
# og_table  : original table (Feature/ASV × sample_N)
# gg2_table : GG2 table (Feature_ID × sample_N)
calculate_read_distribution <- function(og_table, gg2_table) {
  total_reads <- og_table %>%
    select(starts_with("sample_")) %>%
    summarise(across(everything(), ~ sum(.x, na.rm = TRUE))) %>%
    pivot_longer(
      everything(),
      names_to   = "sample",
      values_to  = "total_reads",
      names_prefix = "sample_"
    ) %>%
    mutate(sample = as.numeric(sample))

  gg2_mapped_reads <- gg2_table %>%
    select(starts_with("sample_")) %>%
    summarise(across(everything(), ~ sum(.x, na.rm = TRUE))) %>%
    pivot_longer(
      everything(),
      names_to   = "sample",
      values_to  = "gg2_mapped_reads",
      names_prefix = "sample_"
    ) %>%
    mutate(sample = as.numeric(sample))

  read_distribution <- total_reads %>%
    inner_join(gg2_mapped_reads, by = "sample") %>%
    mutate(
      unmapped_reads     = total_reads - gg2_mapped_reads,
      gg2_mapped_fraction = gg2_mapped_reads / total_reads,
      unmapped_fraction   = unmapped_reads / total_reads
    )

  read_distribution
}

########################################
## Plotting                          ##
########################################

create_occupancy_distribution_plot <- function(occupancy_data, perc_identity) {
  min_occ <- min(occupancy_data$occupancy)
  max_occ <- max(occupancy_data$occupancy)
  bin_edges <- seq(min_occ - 0.5, max_occ + 0.5, by = 1)

  ggplot(occupancy_data, aes(x = occupancy)) +
    geom_histogram(breaks = bin_edges, fill = "steelblue",
                   alpha = 0.7, color = "black") +
    geom_vline(xintercept = mean(occupancy_data$occupancy),
               linetype = "dashed", color = "gray50", linewidth = 1) +
    labs(
      title = paste0("Distribution of ASV/OTU Occupancy per GG2 Feature (perc-identity: ", perc_identity, ")"),
      subtitle = paste(
        "Mean:", round(mean(occupancy_data$occupancy), 2),
        "| Max:", max(occupancy_data$occupancy),
        "| Min:", min(occupancy_data$occupancy)
      ),
      x = "Number of ASVs/OTUs per GG2 Feature",
      y = "Count"
    ) +
    theme_minimal() +
    theme(
      plot.title    = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.title    = element_text(size = 12)
    )
}

create_read_distribution_plot <- function(read_distribution, use_fraction = TRUE) {
  total_gg2_pct <- sum(read_distribution$gg2_mapped_reads) /
                   sum(read_distribution$total_reads) * 100
  total_unmapped_pct <- sum(read_distribution$unmapped_reads) /
                        sum(read_distribution$total_reads) * 100

  subtitle_text <- paste0(
    "Total: ",
    round(total_gg2_pct, 1), "% GG2 Mapped | ",
    round(total_unmapped_pct, 1), "% GG2 Unmapped"
  )

  if (use_fraction) {
    plot_data <- read_distribution %>%
      select(sample, gg2_mapped_fraction, unmapped_fraction) %>%
      pivot_longer(
        c(gg2_mapped_fraction, unmapped_fraction),
        names_to = "category",
        values_to = "fraction"
      ) %>%
      mutate(category = ifelse(category == "gg2_mapped_fraction",
                               "GG2 Mapped", "Unmapped"))

    ggplot(plot_data, aes(x = sample, y = fraction, fill = category)) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_manual(values = c("GG2 Mapped" = "steelblue",
                                   "Unmapped"   = "orange")) +
      labs(
        title    = "Read Distribution by Sample (Fractions)",
        subtitle = subtitle_text,
        x        = "Sample",
        y        = "Fraction of Reads",
        fill     = "Category"
      ) +
      theme_minimal()
  } else {
    plot_data <- read_distribution %>%
      select(sample, gg2_mapped_reads, unmapped_reads) %>%
      pivot_longer(
        c(gg2_mapped_reads, unmapped_reads),
        names_to = "category",
        values_to = "count"
      ) %>%
      mutate(category = ifelse(category == "gg2_mapped_reads",
                               "GG2 Mapped", "Unmapped"))

    ggplot(plot_data, aes(x = sample, y = count, fill = category)) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_manual(values = c("GG2 Mapped" = "steelblue",
                                   "Unmapped"   = "orange")) +
      labs(
        title    = "Read Distribution by Sample (Counts)",
        subtitle = subtitle_text,
        x        = "Sample",
        y        = "Number of Reads",
        fill     = "Category"
      ) +
      theme_minimal()
  }
}

create_read_distribution_scatter <- function(og_table, gg2_table) {
  og_reads <- og_table %>%
    select(starts_with("sample_")) %>%
    summarise(across(everything(), ~ sum(.x, na.rm = TRUE))) %>%
    pivot_longer(
      everything(),
      names_to   = "sample",
      values_to  = "reads",
      names_prefix = "sample_"
    ) %>%
    mutate(sample = as.numeric(sample), category = "Original")

  gg2_reads <- gg2_table %>%
    select(starts_with("sample_")) %>%
    summarise(across(everything(), ~ sum(.x, na.rm = TRUE))) %>%
    pivot_longer(
      everything(),
      names_to   = "sample",
      values_to  = "reads",
      names_prefix = "sample_"
    ) %>%
    mutate(sample = as.numeric(sample), category = "GG2")

  gg2_filtered <- gg2_table %>%
    mutate(across(starts_with("sample_"), ~ ifelse(. < 10, 0, .))) %>%
    filter(rowSums(select(., starts_with("sample_")) > 0) > 1) %>%
    select(starts_with("sample_")) %>%
    summarise(across(everything(), ~ sum(.x, na.rm = TRUE))) %>%
    pivot_longer(
      everything(),
      names_to   = "sample",
      values_to  = "reads",
      names_prefix = "sample_"
    ) %>%
    mutate(sample = as.numeric(sample), category = "GG2 Filtered")

  scatter_data <- bind_rows(og_reads, gg2_reads, gg2_filtered)

  ggplot(scatter_data, aes(x = sample, y = reads, color = category)) +
    geom_point(alpha = 0.7, size = 1.5) +
    scale_color_manual(values = c(
      "Original"      = "black",
      "GG2"           = "steelblue",
      "GG2 Filtered"  = "#934cbf"
    )) +
    labs(
      subtitle = "Comparison of original, GG2 mapped, and filtered GG2 reads",
      x        = "Sample ID",
      y        = "Number of Reads",
      color    = "Data Type"
    ) +
    theme_minimal() +
    theme(
      plot.title    = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.title    = element_text(size = 12),
      legend.position = "right"
    )
}

create_tree_plot <- function(tree, gg2_taxonomy_parsed, top_n = 10) {
  top_phyla <- gg2_taxonomy_parsed %>%
    count(Phylum, sort = TRUE) %>%
    head(top_n) %>%
    pull(Phylum)

  # Phylum per tip: match tree tips to taxonomy; tips not in lookup become NA
  p <- gg2_taxonomy_parsed$Phylum[match(tree$tip.label, gg2_taxonomy_parsed$Feature_ID)]
  n_unmapped <- sum(is.na(p))
  if (n_unmapped > 0) {
    warning(n_unmapped, " tip(s) from the tree are not in the taxonomy lookup; treating as \"Other\".")
  }
  p[is.na(p)] <- "Other"
  # Collapse phyla not in the top-n set to "Other"
  p[!p %in% top_phyla] <- "Other"
  phylum_df <- data.frame(
    PhylumTop = factor(p, levels = c(top_phyla, "Other")),
    row.names = tree$tip.label
  )

  phylum_cols <- c(
    scales::hue_pal(l = 70, c = 100)(length(top_phyla)),
    "grey80"
  )
  names(phylum_cols) <- c(top_phyla, "Other")

  # Use linewidth (not size) for branch thickness to avoid ggplot2 3.4+ deprecation.
  # If you see "@mapping must be <ggplot2::mapping>" or fortify() errors, either:
  # - upgrade ggtree (BiocManager::install("ggtree")) to a version compatible with ggplot2 3.4+,
  # - or pin ggplot2 to < 3.4 in your env (e.g. r-ggplot2=3.3.6). See README or env docs.
  circ <- ggtree(tree, layout = "circular", linewidth = 0.1)
  gheatmap(
    circ,
    phylum_df,
    offset = 0.02,
    width = 0.1,
    colnames = FALSE,
    color = NA
  ) +
    scale_fill_manual(
      values = phylum_cols,
      breaks = c(top_phyla, "Other"),
      name = "Phyla"
    ) +
    theme(legend.position = "right") +
    labs(
      title    = paste0("GG2 Phylogenetic Tree — Top ", top_n, " Phyla"),
      subtitle = paste("Total tips:", length(tree$tip.label))
    )
}

