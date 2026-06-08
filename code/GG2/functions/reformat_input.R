suppressPackageStartupMessages({
  library(tidyverse)
  library(Biostrings)
})

#' Reformat an input table + sequences into a BIOM-ready TSV and filtered FASTA.
#'
#' This function standardizes very different upstream formats into a common layout:
#' - First column: `Feature ID`
#' - Remaining columns: `sample_1`, `sample_2`, ..., `sample_N`
#'
#' It also creates (and returns) a sample mapping data.frame with columns:
#' - `sample_N`         : integer index (1..N)
#' - `env_sample_id`: sample identifier used in env data
#' - `sample_label`     : (TARA only) original label from the miTAG table header
#'
#' @param format One of "mitag" or "dada2"
#' @param input_table Path to the input count table
#' @param input_sequences Path to the input sequences FASTA/FNA file
#' @param output_table Path to write the BIOM-ready TSV
#' @param output_sequences Path to write the filtered sequences FASTA/FNA
#' @param sample_mapping_cfg Optional list with elements
#'   - lookup_file, sample_label_column, analysis_id_column (for TARA/mitag)
#' @param verbose Logical; print progress messages
#'
#' @return A list with:
#'   - table_path: path to BIOM-ready table
#'   - sequences_path: path to filtered sequences
#'   - sample_mapping: data.frame with sample_N/env_sample_id/(sample_label)
reformat_input <- function(format,
                           input_table,
                           input_sequences,
                           output_table,
                           output_sequences,
                           sample_mapping_cfg = NULL,
                           verbose = FALSE) {

  if (!file.exists(input_table)) {
    stop("Input table not found: ", input_table)
  }
  if (!file.exists(input_sequences)) {
    stop("Input sequences not found: ", input_sequences)
  }

  dir_out <- dirname(output_table)
  if (!dir.exists(dir_out)) {
    dir.create(dir_out, recursive = TRUE)
  }

  if (verbose) {
    message("[reformat_input] Reading table from: ", input_table)
  }

  if (format == "mitag") {
    result <- .reformat_mitag(input_table, input_sequences,
                               output_table, output_sequences,
                               sample_mapping_cfg, verbose)
  } else if (format == "dada2") {
    result <- .reformat_dada2(input_table, input_sequences,
                              output_table, output_sequences,
                              verbose)
  } else {
    stop("Unsupported format in reformat_input(): ", format,
         " (expected 'mitag' or 'dada2')")
  }

  if (verbose) {
    message("[reformat_input] Wrote BIOM-ready table to: ", result$table_path)
    message("[reformat_input] Wrote filtered sequences to: ", result$sequences_path)
  }

  return(result)
}

############################
## Internal helper: miTAG ##
############################

.reformat_mitag <- function(input_table,
                            input_sequences,
                            output_table,
                            output_sequences,
                            sample_mapping_cfg,
                            verbose = FALSE) {
  # miTAG layout:
  # Domain, Phylum, Class, Order, Family, Genus, OTU.rep, sample_1, sample_2, ...
  mitag <- read_tsv(input_table, col_types = cols(.default = "c"))

  if (!all(c("Domain", "Phylum", "Class", "Order", "Family", "Genus") %in% colnames(mitag))) {
    stop("Expected first 6 columns to be Domain, Phylum, Class, Order, Family, Genus for miTAG input.")
  }
  if (!"OTU.rep" %in% colnames(mitag)) {
    stop("Expected 'OTU.rep' column in miTAG input table.")
  }

  # Sample columns start after taxonomy + OTU.rep
  tax_cols <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
  otu_col  <- "OTU.rep"
  sample_cols <- setdiff(colnames(mitag), c(tax_cols, otu_col))

  if (length(sample_cols) == 0) {
    stop("No sample columns detected in miTAG table.")
  }

  if (verbose) {
    verbose_print(
      paste0(
        "[reformat_input:mitag] Found ",
        nrow(mitag),
        " OTU rows and ",
        length(sample_cols),
        " sample columns"
      ),
      verbose
    )
  }

  # Map sample labels to analysis IDs using the provided lookup, if available
  if (is.null(sample_mapping_cfg)) {
    stop("sample_mapping_cfg must be provided for 'mitag' format (to map sample_label -> analysis_id).")
  }
  lookup_path <- sample_mapping_cfg$lookup_file
  label_col   <- sample_mapping_cfg$sample_label_column
  id_col      <- sample_mapping_cfg$analysis_id_column

  if (!file.exists(lookup_path)) {
    stop("Sample lookup file not found: ", lookup_path)
  }
  lookup <- read_csv(lookup_path, show_col_types = FALSE)
  if (!all(c(label_col, id_col) %in% colnames(lookup))) {
    stop("Lookup file must contain columns: ", label_col, " and ", id_col)
  }

  # Ensure every sample column can be mapped
  if (!all(sample_cols %in% lookup[[label_col]])) {
    missing <- setdiff(sample_cols, lookup[[label_col]])
    stop("The following sample labels from the miTAG table were not found in the lookup file: ",
         paste(missing, collapse = ", "))
  }

  # Order by env_sample_id so sample_1, sample_2, ... follow ascending env id order
  lookup_sub <- lookup %>%
    filter(.data[[label_col]] %in% sample_cols) %>%
    distinct(.data[[label_col]], .data[[id_col]]) %>%
    arrange(.data[[id_col]])

  # Assign sample_N indices (ascending 1, 2, 3, ...) in env_sample_id order
  sample_mapping <- tibble(
    sample_N           = paste0("sample_", seq_len(nrow(lookup_sub))),
    sample_label       = lookup_sub[[label_col]],
    env_sample_id = as.character(lookup_sub[[id_col]])
  )

  # Build a named vector mapping sample_label -> sample_N colname
  new_names <- sample_mapping$sample_N
  names(new_names) <- sample_mapping$sample_label

  # Remove the taxonomy columns and then rename the OTU column to the expected format ("Feature ID")
  out <- mitag %>%
    select(-all_of(tax_cols)) %>%
    dplyr::rename(`Feature ID` = !!sym(otu_col)) %>%
    mutate(across(all_of(sample_cols), ~ suppressWarnings(as.numeric(.))))

  n_na <- out %>%
    select(all_of(sample_cols)) %>%
    summarise(across(everything(), ~sum(is.na(.)))) %>%
    sum()
  if (n_na > 0) {
    warning("[reformat_input:mitag] Number of NA values in counts: ", n_na, " (replacing with 0s)")
  }
  out <- out %>%
    mutate(across(all_of(sample_cols), ~replace(., is.na(.), 0)))

  # Rename sample columns to sample_N (preserves column order)
  out <- out %>%
    rename_with(~ unname(new_names[.]), .cols = all_of(sample_cols)) %>%
    filter(tolower(trimws(`Feature ID`)) != "unclassified") # Filter out "unclassified" OTUs (case/whitespace insensitive)

  # Order table columns by sample_N numeric value (sample_1, sample_2, ..., sample_10, not alphanumeric)
  sample_n_cols <- setdiff(colnames(out), "Feature ID")
  sample_n_cols_sorted <- sample_n_cols[order(as.integer(sub("^sample_", "", sample_n_cols)))]
  out <- out %>% select(`Feature ID`, all_of(sample_n_cols_sorted))

  # Write BIOM-ready table
  write_tsv(out, output_table)

  # Filter sequences to those present in the table
  seqs <- readDNAStringSet(input_sequences)
  keep_ids <- out$`Feature ID`
  # Match on first token of FASTA headers (headers often are "ID taxonomy..." so full name != table Feature ID)
  seq_id <- if (length(seqs) > 0L) vapply(strsplit(names(seqs), "\\s+"), `[`, 1L, FUN.VALUE = "") else character(0L)
  seqs_filt <- seqs[seq_id %in% keep_ids]
  names(seqs_filt) <- seq_id[seq_id %in% keep_ids]
  writeXStringSet(seqs_filt, output_sequences)

  if (verbose) {
    message("[reformat_input:mitag] BIOM-ready table: ", nrow(out), " features, ", ncol(out),
            " columns (1 Feature ID + ", ncol(out) - 1L, " samples)")
    message("[reformat_input:mitag] ", length(seqs_filt), " sequences saved")
  }

  list(
    table_path     = output_table,
    sequences_path = output_sequences,
    sample_mapping = sample_mapping
  )
}

#############################
## Internal helper: DADA2  ##
#############################

.reformat_dada2 <- function(input_table,
                            input_sequences,
                            output_table,
                            output_sequences,
                            verbose = FALSE) {
  # DADA2 layout (merged table after flow cells):
  # sample_id, ASV_1, ASV_2, ..., ASV_N
  asv_table <- read_csv(input_table, col_types = cols(.default = "c"))

  if (!"sample_id" %in% colnames(asv_table)) {
    stop("Expected 'sample_id' column in DADA2 input table.")
  }

  asv_cols <- setdiff(colnames(asv_table), "sample_id")
  if (length(asv_cols) == 0) {
    stop("No ASV columns found in DADA2 table.")
  }

  if (verbose) {
    message("[reformat_input:dada2] Found ", length(asv_cols),
            " ASV columns and ", nrow(asv_table), " sample rows")
  }

  # Long format: one row per (sample_id, ASV_ID)
  asv_long <- asv_table %>%
    pivot_longer(cols = all_of(asv_cols),
                 names_to = "Feature ID",
                 values_to = "reads") %>%
    mutate(reads = as.numeric(reads))

  # Count the number of NAs in 'reads', print a warning if present, then change them to 0
  num_nas <- sum(is.na(asv_long$reads))
  if (num_nas > 0) {
    warning("[reformat_input:dada2] Found ", num_nas, " NA values in read counts; setting them to 0.")
    asv_long <- asv_long %>% mutate(reads = ifelse(is.na(reads), 0, reads))
  }

  # Wide: rows = features, cols = sample_*
  asv_wide <- asv_long %>%
    pivot_wider(names_from = sample_id,
                values_from = reads,
                values_fill = 0) %>%
    arrange(`Feature ID`)

  # Build sample mapping (sample_N -> sample_id used in env data)
  sample_ids <- setdiff(colnames(asv_wide), "Feature ID")
  sample_mapping <- tibble(
    sample_N           = paste0("sample_",seq_len(length(sample_ids))),
    env_sample_id = sample_ids
  )

  # Rename sample columns to sample_N
  new_sample_names <- sample_mapping$sample_N
  names(new_sample_names) <- sample_ids
  asv_wide_renamed <- asv_wide
  asv_wide_renamed <- asv_wide_renamed %>%
    rename_with(~ new_sample_names[.], .cols = -c(`Feature ID`))

  # Write BIOM-ready table
  write_tsv(asv_wide_renamed, output_table)

  # Filter sequences to those present in the table
  seqs <- readDNAStringSet(input_sequences)
  keep_ids <- asv_wide_renamed$`Feature ID`
  seqs_filt <- seqs[names(seqs) %in% keep_ids]
  writeXStringSet(seqs_filt, output_sequences)

  if (verbose) {
    message("[reformat_input:dada2] BIOM-ready table: ", nrow(asv_wide_renamed), " features, ",
            ncol(asv_wide_renamed), " columns (1 Feature ID + ", ncol(asv_wide_renamed) - 1L, " samples)")
    message("[reformat_input:dada2] ", length(seqs_filt), " sequences saved")
  }

  list(
    table_path     = output_table,
    sequences_path = output_sequences,
    sample_mapping = sample_mapping
  )
}

