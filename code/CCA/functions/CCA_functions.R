# =============================================================================
# CCA_functions.R — Centralized CCA pipeline shared functions
# =============================================================================
# Functions for Regularized CCA: data loading (envdata.csv + composition),
# cross-validation, loadings extraction, alignment, shuffle_sample null.
# =============================================================================

library(tidyverse)
library(RCCA)
library(RColorBrewer)
library(patchwork)
library(parallel)

cov_palette <- brewer.pal(8, "Set2")

# =============================================================================
# LOADING VECTOR ALIGNMENT
# =============================================================================

align_loading_vectors <- function(loadings_df, reference_method = "first_fold") {
  canonical_directions <- sort(unique(loadings_df$canonical_direction))
  folds <- sort(unique(loadings_df$fold))

  flipped_folds <- vector("list", length(canonical_directions))
  names(flipped_folds) <- as.character(canonical_directions)

  aligned_loadings <- loadings_df

  for (CD in canonical_directions) {
    CD_data <- loadings_df %>% filter(canonical_direction == CD)

    # Build reference vector (named by var) so we can align vars safely
    if (reference_method == "first_fold") {
      ref_fold <- min(folds)
      ref_tbl <- CD_data %>%
        filter(fold == ref_fold) %>%
        select(var, value)
      ref_values <- setNames(ref_tbl$value, ref_tbl$var)
    } else if (reference_method == "consensus") {
      ref_tbl <- CD_data %>%
        group_by(var) %>%
        summarise(
          sign_direction = if_else(sum(value > 0) >= sum(value < 0), 1, -1),
          value = sign_direction * mean(abs(value)),
          .groups = "drop"
        ) %>%
        select(var, value)
      ref_values <- setNames(ref_tbl$value, ref_tbl$var)
    } else {
      stop("[align_loading_vectors] reference_method must be 'first_fold' or 'consensus'.")
    }

    to_flip <- integer(0)

    for (f in folds) {
      fold_tbl <- CD_data %>%
        filter(fold == f) %>%
        select(var, value)

      fold_values <- setNames(fold_tbl$value, fold_tbl$var)

      # Align by variable names (safer than arrange(var) if anything is missing/misaligned)
      common_vars <- intersect(names(fold_values), names(ref_values))
      if (length(common_vars) < 2L) {
        stop(
          "[align_loading_vectors] Fold ", f, " (canonical_direction ", CD, ") has too few variables in common with the reference. ",
          "Reference has ", length(ref_values), " vars; fold has ", length(fold_values), " vars; common: ", length(common_vars), ".",
          call. = FALSE
        )
      }

      r <- stats::cor(fold_values[common_vars], ref_values[common_vars])

      # flip if the correlation between current fold and the reference is negative
      if (is.finite(r) && r < 0) {
        to_flip <- c(to_flip, f)
      }
    }
    flipped_folds[[as.character(CD)]] <- to_flip

    # Apply flips for this CD in one shot (vectorized)
    if (length(to_flip) > 0) {
      aligned_loadings <- aligned_loadings %>%
        dplyr::mutate(
          value = dplyr::if_else(
            canonical_direction == CD & fold %in% to_flip,
            -value,
            value
          )
        )
    }
  }

  list(aligned_loadings = aligned_loadings, flipped_folds = flipped_folds)
}

apply_fold_flips <- function(loadings_df, flipped_folds) {
  aligned_loadings <- loadings_df
  for (CD in names(flipped_folds)) {
    CD_int <- as.integer(CD)
    for (f in flipped_folds[[CD]]) {
      aligned_loadings <- aligned_loadings %>%
        mutate(value = ifelse(canonical_direction == CD_int & fold == f, -value, value))
    }
  }
  return(aligned_loadings)
}

# =============================================================================
# NORMALIZE AND ALIGN LOADINGS (shared by steps 2 and 3)
# =============================================================================

#' L2-normalise per (fold x canonical_direction) and sign-align loading vectors.
#'
#' @param x_ld Long data frame: var, fold, canonical_direction, value (abundance loadings).
#' @param y_ld Long data frame: same columns (environment loadings).
#' @param reference_method Passed to align_loading_vectors(). Default "first_fold".
#' @return List: x_loadings and y_loadings (normalised and aligned);
#'   y_prealign (normalised, pre-alignment, for diagnostic plots); flipped_folds.
normalize_and_align_loadings <- function(x_ld, y_ld, reference_method = "first_fold") {
  norm_df <- function(df) {
    df %>%
      group_by(canonical_direction, fold) %>%
      mutate(value = value / sqrt(sum(value^2))) %>%
      ungroup()
  }
  x_norm  <- norm_df(x_ld)
  y_norm  <- norm_df(y_ld)
  y_align <- align_loading_vectors(y_norm, reference_method = reference_method)
  list(
    x_loadings    = apply_fold_flips(x_norm, y_align$flipped_folds),
    y_loadings    = y_align$aligned_loadings,
    y_prealign    = y_norm,
    flipped_folds = y_align$flipped_folds
  )
}

# =============================================================================
# DATA LOADING — env from envdata.csv only; drop NaN rows; composition as before
# =============================================================================

#' Load environment and composition data for CCA.
#'
#' Env: expect clean format after read_data (samples x features + 1)
#' > this code reads the file in, keeps all cols except sample_id, drop NaN rows.
#' Composition: expect cleaned seqtab at the end of GG2 (samples x features + 1) 
#' > this code calculates relative abundance applies a rare-taxon filter
#'
#' @param env_fn Path to envdata.csv (from read_data).
#' @param composition_fn Path to composition table (seqtab.csv).
#' @param read_threshold Minimum read count for rare-taxon filter; taxa in ≤1 sample above this are dropped.
#' @param verbose If TRUE, print dimensions and progress.
#' @return List with \code{env_df} (environmental data, samples x variables) and \code{composition_df} (compositional data, relative abundance, samples x features).
load_datasets_cca <- function(env_fn, composition_fn, read_threshold = 10, verbose = FALSE) {
  if (!file.exists(env_fn)) {
    stop("[load_datasets_cca] envdata.csv not found; run read_data first.\n  Expected: ", env_fn, call. = FALSE)
  }
  env_df <- read_csv(env_fn, show_col_types = FALSE)
  if (!"sample_id" %in% names(env_df)) {
    stop("[load_datasets_cca] envdata.csv must contain a column 'sample_id'. Run read_data first.", call. = FALSE)
  }
  if (verbose) {
    cat(paste0("[load_datasets_cca] Environmental data: ", nrow(env_df), " rows (samples), ", ncol(env_df), " columns (features + 1).\n"))
  }
  env_vars <- setdiff(names(env_df), "sample_id")
  n_samples_before <- nrow(env_df)
  env_df <- env_df %>% 
    filter(if_all(all_of(env_vars), ~ !is.na(.x))) %>%
    mutate(sample_id = as.character(sample_id)) # make sure consistent dtype for comparison with compositional data
  if (nrow(env_df) < n_samples_before) {
    warning("[load_datasets_cca] Dropped ", n_samples_before - nrow(env_df), " rows with NA; ", nrow(env_df), " samples retained.")
  }
  n_samples_env <- nrow(env_df)

  composition_df <- read_csv(composition_fn, show_col_types = FALSE)
  if (!"sample_id" %in% names(composition_df)) {
    stop("[load_datasets_cca] Expecting seqtab format: rows = samples, first column 'sample_id', then features. Check file format.", call. = FALSE)
  }
  if (verbose) {
    cat(paste0("[load_datasets_cca] Compositional data: ", nrow(composition_df), " rows (samples), ", ncol(composition_df), " columns (features + 1).\n"))
  }
  composition_df <- composition_df %>%
    mutate(sample_id = as.character(sample_id))  # make sure consistent dtype for comparison with environmental data
  n_samples_composition <- nrow(composition_df)
  
  # Find the overlapping samples between environmental and compositional data
  sample_overlap <- intersect(env_df$sample_id, composition_df$sample_id)
  if (length(sample_overlap) == 0) {
    stop("[load_datasets_cca] No overlapping sample_id found between environmental and compositional data.", call. = FALSE)
  }
  if (verbose) {
    cat(paste0("[load_datasets_cca] ", length(sample_overlap), " samples shared between environmental (", n_samples_env, " samples) and compositional (", n_samples_composition, " samples) data.\n"))
  }
  # Filter both datasets to only those samples and ensure ordering is the same
  env_df <- env_df %>%
    filter(sample_id %in% sample_overlap) %>%
    arrange(sample_id)
  composition_df <- composition_df %>%
    filter(sample_id %in% sample_overlap) %>%
    arrange(sample_id)
  if (!identical(env_df$sample_id, composition_df$sample_id)) {
    stop("[load_datasets_cca] Environmental and compositional sample_id order or set mismatch.", call. = FALSE)
  }

  # Remove rare tax (i.e. those that are present in <= 1 sample)
  rare_taxa <- composition_df %>%
    pivot_longer(cols = -sample_id, names_to = "taxon", values_to = "reads") %>%
    group_by(taxon) %>%
    summarize(nsamples = sum(reads > read_threshold, na.rm = TRUE), .groups = "drop") %>%
    filter(nsamples <= 1) %>%
    pull(taxon)
  if (length(rare_taxa) > 0) {
    if (verbose) {
      cat(paste0("[load_datasets_cca] Filtering out ", length(rare_taxa), " of ", ncol(composition_df) - 1, " total taxa.\n"))
    }
    composition_df <- composition_df %>% select(-any_of(rare_taxa))
  }
  # Calculate relative abundance
  composition_df <- composition_df %>%
    pivot_longer(cols = -sample_id, names_to = "taxon", values_to = "reads") %>%
    group_by(sample_id) %>%
    mutate(rel_abundance = reads / sum(reads)) %>%
    ungroup() %>%
    select(-reads) %>%
    pivot_wider(names_from = taxon, values_from = rel_abundance, values_fill = 0)

  if (verbose) {
    cat(paste0("[load_datasets_cca] Final:\n  ", nrow(env_df), " samples;\n  env_df ", ncol(env_df) - 1L, " variables;\n  composition_df ", ncol(composition_df) - 1L, " features (relative abundance).\n"))
  }
  list(env_df = env_df, composition_df = composition_df)
}

# =============================================================================
# FOLD TRAIN/TEST SPLITS
# =============================================================================

#' Return train and test sample indices for a given fold.
#'
#' @param fold_indices Integer vector of length n (sample index -> fold number).
#' @param fold_number Which fold is held out as test (1 to k).
#' @return List with \code{train_index} and \code{test_index} (integer vectors of row indices).
get_fold_split <- function(fold_indices, fold_number) {
  train_index <- which(fold_indices != fold_number)
  test_index <- which(fold_indices == fold_number)
  list(train_index = train_index, test_index = test_index)
}

# =============================================================================
# CROSS-VALIDATION AND LOADINGS
# =============================================================================

correlation <- function(method, X, Y, train_index, test_index) {
  X.train <- X[train_index, ]; Y.train <- Y[train_index, ]
  X.test <- X[test_index, ]; Y.test <- Y[test_index, ]
  solution <- method(X = X.train, Y = Y.train)
  num_variates <- ncol(solution$x.coefs)
  cor.train <- sapply(1:num_variates, function(i) cor(X.train %*% solution$x.coefs[, i], Y.train %*% solution$y.coefs[, i]))
  cor.test <- sapply(1:num_variates, function(i) cor(X.test %*% solution$x.coefs[, i], Y.test %*% solution$y.coefs[, i]))
  list(train = cor.train, test = cor.test)
}

run_one_rcca_cv <- function(X, Y, k, lambda1, lambda2, fold_indices) {
  fold_results_list <- list()
  for (fold in 1:k) {
    split <- get_fold_split(fold_indices, fold)
    train_index <- split$train_index
    test_index <- split$test_index
    method <- function(...) RCCA(..., lambda1 = lambda1, lambda2 = lambda2)
    cor <- correlation(method, X, Y, train_index, test_index)
    n_directions <- length(cor$train)
    fold_df <- data.frame(
      k = k, fold = fold, lambda1 = lambda1, lambda2 = lambda2,
      canonical_direction = seq_len(n_directions),
      cor.train = cor$train, cor.test = cor$test,
      stringsAsFactors = FALSE
    )
    fold_results_list[[length(fold_results_list) + 1]] <- fold_df
  }
  block <- do.call(rbind, fold_results_list)
  direction_rank <- block %>%
    group_by(canonical_direction) %>%
    summarise(mean_test = mean(cor.test, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mean_test)) %>%
    mutate(new_canonical_direction = row_number())
  canonical_direction_map <- setNames(direction_rank$new_canonical_direction, as.character(direction_rank$canonical_direction))
  block$canonical_direction <- canonical_direction_map[as.character(block$canonical_direction)]
  block
}

k_fold_cv_rcca <- function(X, Y, ks, lambda1s, lambda2s = 0, seed = 1, verbose = FALSE, n_cores = NULL) {
  set.seed(seed)
  n <- nrow(X)
  fold_indices_by_k <- list()
  for (k in ks) {
    fold_indices_by_k[[as.character(k)]] <- sample(rep(1:k, length.out = n))
  }
  tasks <- list()
  for (k in ks) {
    fold_indices <- fold_indices_by_k[[as.character(k)]]
    for (lambda1 in lambda1s) {
      for (lambda2 in lambda2s) {
        tasks[[length(tasks) + 1L]] <- list(k = k, lambda1 = lambda1, lambda2 = lambda2, fold_indices = fold_indices)
      }
    }
  }
  n_workers <- if (!is.null(n_cores)) {
    as.integer(n_cores)
  } else {
    slurm <- Sys.getenv("SLURM_CPUS_PER_TASK", "")
    if (nzchar(slurm)) {
      as.integer(slurm)
    } else {
      parallel::detectCores()
    }
  }
  if (is.na(n_workers)) n_workers <- 1L
  n_workers <- max(1L, min(n_workers, length(tasks)))
  if (verbose) {
    if (n_workers == 1L) {
      cat("[k_fold_cv_rcca] Only one worker; running in serial.\n")
    } else {
      cat("[k_fold_cv_rcca] Using", n_workers, "workers for parallelization.\n")
    }
  }
  if (n_workers == 1L) {
    results <- lapply(tasks, function(t) run_one_rcca_cv(X, Y, t$k, t$lambda1, t$lambda2, t$fold_indices))
  } else {
    results <- mclapply(tasks, function(t) run_one_rcca_cv(X, Y, t$k, t$lambda1, t$lambda2, t$fold_indices), mc.cores = n_workers)
  }
  out <- do.call(rbind, results)
  out <- out %>% arrange(k, lambda1, lambda2, fold, canonical_direction)
  if (verbose && n_workers > 1L) {
    cat("[k_fold_cv_rcca] Parallel execution complete with", n_workers, "workers.\n")
  }
  out
}

rcca_loadings <- function(X, Y, k, lambda1, lambda2, seed = 1) {
  set.seed(seed)
  n <- nrow(X)
  fold_indices <- sample(rep(1:k, length.out = n))
  x_loadings_list <- list()
  y_loadings_list <- list()
  cor_list <- list()
  proj_list <- list()
  list_idx <- 1L
  for (fold in 1:k) {
    split <- get_fold_split(fold_indices, fold)
    train_index <- split$train_index
    test_index <- split$test_index
    X_train <- X[train_index, , drop = FALSE]
    Y_train <- Y[train_index, , drop = FALSE]
    X_test <- X[test_index, , drop = FALSE]
    Y_test <- Y[test_index, , drop = FALSE]
    solution <- RCCA(X_train, Y_train, lambda1 = lambda1, lambda2 = lambda2)
    x_loadings <- solution$x.coefs
    y_loadings <- solution$y.coefs
    n_directions <- solution$n.comp
    if (is.null(rownames(x_loadings))) rownames(x_loadings) <- paste0("X_var", seq_len(nrow(x_loadings)))
    if (is.null(rownames(y_loadings))) rownames(y_loadings) <- paste0("Y_var", seq_len(nrow(y_loadings)))
    split_vec <- character(n)
    split_vec[train_index] <- "train"
    split_vec[test_index] <- "test"
    for (CD in seq_len(n_directions)) {
      x_CD <- x_loadings[, CD]
      y_CD <- y_loadings[, CD]
      proj_comp <- c(X %*% x_CD)
      proj_env <- c(Y %*% y_CD)
      train_cor <- cor(proj_comp[train_index], proj_env[train_index])
      test_cor <- cor(proj_comp[test_index], proj_env[test_index])
      proj_list[[list_idx]] <- data.frame(
        sample = seq_len(n),
        fold = fold,
        canonical_direction = CD,
        split = split_vec,
        proj_comp = proj_comp,
        proj_env = proj_env,
        stringsAsFactors = FALSE
      )
      x_loadings_list[[list_idx]] <- data.frame(
        canonical_direction = CD, var = rownames(x_loadings), fold = fold, value = x_CD, stringsAsFactors = FALSE
      )
      y_loadings_list[[list_idx]] <- data.frame(
        canonical_direction = CD, var = rownames(y_loadings), fold = fold, value = y_CD, stringsAsFactors = FALSE
      )
      cor_list[[list_idx]] <- data.frame(canonical_direction = CD, fold = fold, train = train_cor, test = test_cor, stringsAsFactors = FALSE)
      list_idx <- list_idx + 1L
    }
  }
  x_loadings_df <- do.call(rbind, x_loadings_list)
  y_loadings_df <- do.call(rbind, y_loadings_list)
  cor_df <- do.call(rbind, cor_list)
  projections_df <- do.call(rbind, proj_list)
  CD_means <- cor_df %>%
    group_by(canonical_direction) %>%
    summarize(mean_test = mean(test, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mean_test)) %>%
    mutate(new_canonical_direction = row_number())
  CD_map <- setNames(CD_means$new_canonical_direction, as.character(CD_means$canonical_direction))
  x_loadings_df$canonical_direction <- CD_map[as.character(x_loadings_df$canonical_direction)]
  y_loadings_df$canonical_direction <- CD_map[as.character(y_loadings_df$canonical_direction)]
  cor_df$canonical_direction <- CD_map[as.character(cor_df$canonical_direction)]
  projections_df$canonical_direction <- CD_map[as.character(projections_df$canonical_direction)]
  list(
    x_loadings = x_loadings_df,
    y_loadings = y_loadings_df,
    corr = cor_df,
    fold_indices = fold_indices,
    projections = projections_df
  )
}

# =============================================================================
# NULL: shuffle_sample only
# =============================================================================

shuffle_sample <- function(X, seed = 1) {
  if (!is.null(seed)) set.seed(seed)
  X[sample(nrow(X)), , drop = FALSE]
}

# =============================================================================
# PLOTTING HELPERS (slim)
# =============================================================================

#' Hyperparameter CV plot: mean ± SD test (and train) vs log10(lambda1) for all canonical directions; vertical line at best lambda1 (chosen from first direction).
plot_cv_lambda <- function(results, best_lambda1, dir_path, group_label = "") {
  summ <- results %>%
    group_by(k, lambda1, lambda2, canonical_direction) %>%
    summarise(
      mean_train = mean(cor.train, na.rm = TRUE),
      sd_train = sd(cor.train, na.rm = TRUE),
      mean_test = mean(cor.test, na.rm = TRUE),
      sd_test = sd(cor.test, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    pivot_longer(cols = c(mean_train, sd_train, mean_test, sd_test), names_to = c("stat", "type"), names_sep = "_", values_to = "value") %>%
    pivot_wider(names_from = stat, values_from = value) %>%
    mutate(type = factor(type, levels = c("test", "train")))
  if (nrow(summ) == 0) return(invisible(NULL))
  ymin <- min(summ$mean - summ$sd, 0, na.rm = TRUE)
  ymax <- max(summ$mean + summ$sd, 1, na.rm = TRUE)
  n_cov <- n_distinct(summ$canonical_direction)
  p <- ggplot(summ, aes(x = log10(lambda1), y = mean)) +
    geom_vline(xintercept = log10(best_lambda1), linetype = "dashed", color = "red", linewidth = 0.8) +
    geom_line(color = "black") +
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.15, fill = "darkblue") +
    facet_grid(rows = vars(type), cols = vars(canonical_direction),
               labeller = labeller(
                 type = as_labeller(c(test = "Test", train = "Train")),
                 canonical_direction = function(x) ifelse(as.character(x) == "1", "Direction 1 (used for selection)", paste("Direction", x))
               )) +
    labs(x = expression(log[10](lambda[1])), y = "Correlation (mean ± SD)", title = group_label) +
    ylim(ymin, ymax)
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(dir_path, "cv_plot.jpg"), plot = p, width = 4 + 2.5 * min(n_cov, 4), height = 4, dpi = 150)
  invisible(p)
}

# =============================================================================
# ALIGNMENT ASSESSMENT PLOTS
# =============================================================================

#' Plot alignment assessment: original vs aligned loadings (colored by fold).
#'
#' Expects long-format loadings with columns: canonical_direction, var, fold, value.
#' Works for any loadings (environmental or abundance).
#'
#' @param loadings_original Data frame of loadings before alignment (canonical_direction, var, fold, value).
#' @param loadings_aligned Data frame of loadings after alignment (same structure).
#' @param dir_path Directory to save the combined plot.
#' @param filename Base filename for the saved plot (default "alignment_combined.jpg").
#' @return Invisibly, the combined ggplot (patchwork) object.
plot_alignment_assessment <- function(loadings_original,
                                       loadings_aligned,
                                       dir_path,
                                       filename = "alignment_combined.jpg") {
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)

  p_original <- loadings_original %>%
    ggplot(aes(x = var, y = value, color = fold)) +
    geom_point() +
    geom_line(aes(group = factor(fold))) +
    facet_wrap(~ canonical_direction, labeller = labeller(canonical_direction = function(x) paste0("CD ", x))) +
    ylim(-1, 1) +
    labs(x = "Variable", y = "Loading", color = "Fold", title = "Original Vectors") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

  p_aligned <- loadings_aligned %>%
    ggplot(aes(x = var, y = value, color = fold)) +
    geom_point() +
    geom_line(aes(group = factor(fold))) +
    facet_wrap(~ canonical_direction, labeller = labeller(canonical_direction = function(x) paste0("CD ", x))) +
    ylim(-1, 1) +
    labs(x = "Variable", y = "Loading", color = "Fold", title = "Aligned Vectors") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

  combined <- p_original / p_aligned
  ggsave(file.path(dir_path, filename), plot = combined, width = 10, height = 12, dpi = 150)
  invisible(combined)
}
