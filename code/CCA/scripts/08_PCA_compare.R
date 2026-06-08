#!/usr/bin/env Rscript
# 08_PCA_compare.R — Compute angle and EVR comparisons between CCA canonical directions and PCA.
#
# Per tax level, saves to <dataset>/results/CCA/<perc>/<tax_level>/08_PCA_compare/:
#   angle_results.csv  (fold, cd, angle)
#   null_results.csv   (fold, cd, null_mean, null_sd)
#   evr_results.csv    (fold, j, method, split, evr)
#
# Usage:
#   Rscript code/CCA/scripts/08_PCA_compare.R \
#     --dataset soil --n_cds 4 --perc_identity 0.90
#   Rscript code/CCA/scripts/08_PCA_compare.R \
#     --dataset ocean --n_cds 3 --perc_identity 0.90

args_raw <- commandArgs(trailingOnly = FALSE)
file_arg  <- grep("^--file=", args_raw, value = TRUE)
if (length(file_arg) != 1) stop("Run with Rscript so --file= is available.", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)
script_dir  <- dirname(script_path)
setup_path  <- normalizePath(file.path(script_dir, "..", "..", "setup.R"),
                             winslash = "/", mustWork = TRUE)
source(setup_path)

suppressPackageStartupMessages({
  library(optparse)
  library(parallel)
  library(readr)
  library(dplyr)
})

TAX_LEVELS <- c("OTU", "Species", "Genus", "Family", "Order", "Class", "Phylum")

# ---------------------------------------------------------------------------
parse_cli_args <- function(argv = commandArgs(trailingOnly = TRUE)) {
  option_list <- list(
    make_option("--dataset",       type = "character", default = NA,
                help = "Dataset: soil or ocean"),
    make_option("--perc_identity", type = "character", default = "0.90",
                help = "Percent identity [default %default]"),
    make_option("--n_cds",         type = "integer",   default = 4L,
                help = "Number of significant canonical directions [default %default]"),
    make_option("--n_shuffles",    type = "integer",   default = 100L,
                help = "Shuffles per fold for angle null [default %default]"),
    make_option("--n_cores",       type = "integer",   default = 1L,
                help = "Worker processes for fold loop [default %default]"),
    make_option("--blas_threads",  type = "integer",   default = 1L,
                help = "BLAS/OpenMP threads per worker (requires RhpcBLASctl) [default %default]"),
    make_option("--root",          type = "character", default = NA,
                help = "Repo root (inferred if omitted)")
  )
  a <- parse_args(OptionParser(option_list = option_list), args = argv)
  if (is.na(a$dataset) || !nzchar(trimws(a$dataset)))
    stop("--dataset is required", call. = FALSE)
  list(dataset      = trimws(tolower(a$dataset)),
       perc         = a$perc_identity,
       n_cds        = a$n_cds,
       n_shuffles   = a$n_shuffles,
       n_cores      = a$n_cores,
       blas_threads = a$blas_threads,
       root         = a$root)
}

# ---------------------------------------------------------------------------
shuffle_rows <- function(X) t(apply(X, 1, sample))

evr_cum <- function(X_eval, Q, k_max) {
  tv <- sum(X_eval^2)
  if (tv == 0) return(rep(NA_real_, k_max))
  vapply(seq_len(k_max), function(k) {
    sum((X_eval %*% Q[, 1:k, drop = FALSE])^2) / tv
  }, numeric(1L))
}

# ---------------------------------------------------------------------------
run_tax_level <- function(tax, tax_path, n_cds, n_shuffles, n_cores, blas_threads, out_dir) {
  x_file    <- file.path(tax_path, "step0_data",    "X_matrix.csv")
  fi_file   <- file.path(tax_path, "step2_loadings","fold_indices.csv")
  load_file <- file.path(tax_path, "step2_loadings","abundance_loadings.csv")

  if (!all(file.exists(x_file, load_file))) {
    message("  Skipping ", tax, " (missing files)")
    return(invisible(NULL))
  }
  message("  Processing ", tax, " ...")

  X          <- as.matrix(read_csv(x_file,    show_col_types = FALSE))
  loads_raw  <- read_csv(load_file, show_col_types = FALSE)
  feat_names <- colnames(X)

  hyperparam <- read_csv(file.path(tax_path, "step1_hyperparam", "best_hyperparam.csv"),
                         show_col_types = FALSE)
  n_folds <- hyperparam$k[[1]]
  k_use   <- min(n_cds, max(loads_raw$canonical_direction, na.rm = TRUE))

  if (file.exists(fi_file)) {
    fold_vec <- read_csv(fi_file, show_col_types = FALSE)[[1]]
  } else {
    set.seed(1L)
    fold_vec <- sample(rep(seq_len(n_folds), length.out = nrow(X)))
  }

  run_fold <- function(f) {
    set.seed(42L + f)   # deterministic per-fold shuffle seeds across runs
    if (blas_threads > 1L) {
      if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
        RhpcBLASctl::blas_set_num_threads(blas_threads)
        RhpcBLASctl::omp_set_num_threads(blas_threads)
      } else {
        # Fallback: set env vars; honoured by OpenBLAS and most OpenMP runtimes
        Sys.setenv(OPENBLAS_NUM_THREADS = blas_threads,
                   OMP_NUM_THREADS      = blas_threads,
                   MKL_NUM_THREADS      = blas_threads)
      }
    }
    t0 <- proc.time()[["elapsed"]]

    train_idx <- which(fold_vec != f)
    test_idx  <- which(fold_vec == f)
    X_train   <- X[train_idx, , drop = FALSE]
    X_test    <- X[test_idx,  , drop = FALSE]

    pca_fit <- prcomp(X_train, center = FALSE, scale. = FALSE, rank. = k_use)
    PC1     <- pca_fit$rotation[, 1]
    Q_pca   <- pca_fit$rotation

    A_cca <- vapply(seq_len(k_use), function(j) {
      rows <- loads_raw[loads_raw$canonical_direction == j & loads_raw$fold == f,
                        c("var", "value")]
      rows$value[match(feat_names, rows$var)]
    }, numeric(length(feat_names)))
    if (!is.matrix(A_cca)) A_cca <- matrix(A_cca, ncol = 1L)
    Q_cca <- qr.Q(qr(A_cca))
    stopifnot(ncol(Q_cca) == k_use)

    null_pc1_mat <- vapply(seq_len(n_shuffles), function(s) {
      prcomp(shuffle_rows(X_train), center = FALSE, scale. = FALSE, rank. = 1L)$rotation[, 1]
    }, numeric(ncol(X)))

    message("    fold ", f, "/", n_folds,
            " — dim(X_train)=", nrow(X_train), "x", ncol(X_train),
            " — nulls done in ", round(proc.time()[["elapsed"]] - t0, 1), "s")

    angle_rows <- list(); null_rows <- list(); evr_rows <- list()

    for (j in seq_len(k_use)) {
      cd_j      <- A_cca[, j]
      cd_j_norm <- cd_j / sqrt(sum(cd_j^2))
      true_angle <- acos(pmin(1, abs(sum(cd_j_norm * PC1))))

      null_angs <- vapply(seq_len(n_shuffles), function(s) {
        acos(pmin(1, abs(sum(cd_j_norm * null_pc1_mat[, s]))))
      }, numeric(1L))

      angle_rows[[j]] <- data.frame(fold = f, cd = j, angle = true_angle)
      null_rows[[j]]  <- data.frame(fold = f, cd = j,
                                    null_mean = mean(null_angs),
                                    null_sd   = sd(null_angs))
    }

    evr_cca_tr <- evr_cum(X_train, Q_cca, k_use)
    evr_cca_te <- evr_cum(X_test,  Q_cca, k_use)
    evr_pca_tr <- evr_cum(X_train, Q_pca, k_use)
    evr_pca_te <- evr_cum(X_test,  Q_pca, k_use)

    for (j in seq_len(k_use)) {
      evr_rows[[j]] <- data.frame(
        fold   = f, j = j,
        evr    = c(evr_cca_tr[j], evr_cca_te[j], evr_pca_tr[j], evr_pca_te[j]),
        method = c("CCA", "CCA", "PCA", "PCA"),
        split  = c("train", "test", "train", "test")
      )
    }

    list(angle = bind_rows(angle_rows),
         null  = bind_rows(null_rows),
         evr   = bind_rows(evr_rows))
  }

  fold_results <- mclapply(seq_len(n_folds), run_fold, mc.cores = n_cores)

  # Surface any worker errors before attempting to bind
  errs <- which(vapply(fold_results, inherits, logical(1L), "error"))
  if (length(errs) > 0) {
    for (e in errs)
      message("  ERROR in fold ", e, ": ", conditionMessage(fold_results[[e]]))
    stop(sprintf("  %d fold(s) failed for %s — see messages above.", length(errs), tax),
         call. = FALSE)
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv(bind_rows(lapply(fold_results, `[[`, "angle")), file.path(out_dir, "angle_results.csv"))
  write_csv(bind_rows(lapply(fold_results, `[[`, "null")),  file.path(out_dir, "null_results.csv"))
  write_csv(bind_rows(lapply(fold_results, `[[`, "evr")),   file.path(out_dir, "evr_results.csv"))
  message("  Saved results for ", tax, " -> ", out_dir)
}

# ---------------------------------------------------------------------------
main <- function() {
  cli <- parse_cli_args()

  root <- if (is.na(cli$root) || !nzchar(trimws(cli$root))) {
    get("REPO_ROOT", envir = .GlobalEnv)
  } else {
    normalizePath(cli$root, winslash = "/", mustWork = TRUE)
  }

  base_path <- file.path(root, cli$dataset, "results", "CCA", cli$perc)
  message("Dataset: ", cli$dataset, " | n_cds: ", cli$n_cds,
          " | n_shuffles: ", cli$n_shuffles,
          " | n_cores: ", cli$n_cores, " | blas_threads: ", cli$blas_threads)

  for (tax in TAX_LEVELS) {
    tax_path <- file.path(base_path, tax)
    out_dir  <- file.path(tax_path, "step8_PCA_compare")
    run_tax_level(tax, tax_path, cli$n_cds, cli$n_shuffles,
                  cli$n_cores, cli$blas_threads, out_dir)
  }

  message("Done.")
}

main()
