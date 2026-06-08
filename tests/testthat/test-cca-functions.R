library(testthat)
library(dplyr)
source(file.path(Sys.getenv("REPO_ROOT", "."), "code/CCA/functions/CCA_functions.R"))

# ---- get_fold_split ----

test_that("get_fold_split returns correct test index for fold 1", {
  fold_indices <- c(1, 2, 1, 2, 1, 2)
  split <- get_fold_split(fold_indices, fold_number = 1)
  expect_equal(sort(split$test_index),  c(1L, 3L, 5L))
  expect_equal(sort(split$train_index), c(2L, 4L, 6L))
})

test_that("get_fold_split train and test cover all samples exactly once", {
  n <- 20
  fold_indices <- rep(1:5, length.out = n)
  for (f in 1:5) {
    split <- get_fold_split(fold_indices, fold_number = f)
    expect_length(intersect(split$train_index, split$test_index), 0)
    expect_length(union(split$train_index, split$test_index), n)
  }
})

# ---- align_loading_vectors ----

test_that("align_loading_vectors leaves already-aligned folds unchanged", {
  df <- data.frame(
    canonical_direction = rep(1L, 6),
    fold  = rep(1:3, each = 2),
    var   = rep(c("A", "B"), 3),
    value = c(0.5, 0.3, 0.4, 0.2, 0.6, 0.4)
  )
  result <- align_loading_vectors(df)
  expect_equal(result$aligned_loadings$value, df$value)
  expect_length(result$flipped_folds[["1"]], 0)
})

test_that("align_loading_vectors flips fold with opposite sign", {
  df <- data.frame(
    canonical_direction = rep(1L, 6),
    fold  = rep(1:3, each = 2),
    var   = rep(c("A", "B"), 3),
    value = c(0.5, 0.3,    # fold 1: reference
             -0.5, -0.3,   # fold 2: flipped relative to fold 1
              0.4, 0.2)    # fold 3: same sign as fold 1
  )
  result <- align_loading_vectors(df)
  expect_true(2L %in% result$flipped_folds[["1"]])
  expect_false(1L %in% result$flipped_folds[["1"]])
  expect_false(3L %in% result$flipped_folds[["1"]])
  aligned_fold2 <- result$aligned_loadings %>%
    filter(fold == 2L, canonical_direction == 1L) %>%
    pull(value)
  expect_true(all(aligned_fold2 > 0))
})

test_that("align_loading_vectors works with consensus reference_method", {
  df <- data.frame(
    canonical_direction = rep(1L, 6),
    fold  = rep(1:3, each = 2),
    var   = rep(c("A", "B"), 3),
    value = c(0.5, 0.3, -0.5, -0.3, 0.4, 0.2)
  )
  result <- align_loading_vectors(df, reference_method = "consensus")
  expect_true(2L %in% result$flipped_folds[["1"]])
})

test_that("apply_fold_flips produces same result as align_loading_vectors", {
  df <- data.frame(
    canonical_direction = rep(1L, 6),
    fold  = rep(1:3, each = 2),
    var   = rep(c("A", "B"), 3),
    value = c(0.5, 0.3, -0.5, -0.3, 0.4, 0.2)
  )
  aligned  <- align_loading_vectors(df)
  manually <- apply_fold_flips(df, aligned$flipped_folds)
  expect_equal(manually$value, aligned$aligned_loadings$value)
})

test_that("normalize_and_align_loadings produces unit-L2-norm per fold/canonical_direction", {
  set.seed(1)
  mk <- function(vals) data.frame(
    canonical_direction = rep(1L, 8),
    fold  = rep(1:4, each = 2),
    var   = rep(c("A", "B"), 4),
    value = vals
  )
  df     <- mk(rnorm(8))
  result <- normalize_and_align_loadings(df, df)
  norms  <- result$x_loadings |>
    dplyr::group_by(canonical_direction, fold) |>
    dplyr::summarise(norm = sqrt(sum(value^2)), .groups = "drop")
  expect_true(all(abs(norms$norm - 1) < 1e-10))
})

test_that("normalize_and_align_loadings applies y-derived flip to x loadings", {
  # Values must differ within each fold so cor() is defined after L2 normalisation.
  # (0.8, 0.6) is already unit-norm; fold 2 is exactly anti-parallel → r = -1 → flip.
  y_df <- data.frame(
    canonical_direction = rep(1L, 6),
    fold  = rep(1:3, each = 2),
    var   = rep(c("E1", "E2"), 3),
    value = c( 0.8,  0.6,   # fold 1: reference
              -0.8, -0.6,   # fold 2: anti-parallel → will be flipped
               0.8,  0.6)   # fold 3: same as reference
  )
  x_df <- data.frame(
    canonical_direction = rep(1L, 6),
    fold  = rep(1:3, each = 2),
    var   = rep(c("T1", "T2"), 3),
    value = c(0.9, 0.3, 0.9, 0.3, 0.9, 0.3)  # all positive before alignment
  )
  result  <- normalize_and_align_loadings(x_df, y_df)
  x_fold2 <- result$x_loadings |>
    dplyr::filter(fold == 2L) |>
    dplyr::pull(value)
  expect_true(all(x_fold2 < 0))
})

test_that("align_loading_vectors flips canonical directions independently", {
  df <- data.frame(
    canonical_direction = c(rep(1L, 6), rep(2L, 6)),
    fold  = c(rep(1:3, each = 2), rep(1:3, each = 2)),
    var   = c(rep(c("A", "B"), 3), rep(c("A", "B"), 3)),
    value = c( 0.5,  0.3,   # CD1 fold 1: reference
              -0.5, -0.3,   # CD1 fold 2: opposite → flip
               0.5,  0.3,   # CD1 fold 3: matches reference
              -0.5, -0.3,   # CD2 fold 1: reference (negative)
              -0.5, -0.3,   # CD2 fold 2: matches reference → no flip
              -0.5, -0.3)   # CD2 fold 3: matches reference → no flip
  )
  result <- align_loading_vectors(df)
  expect_true( 2L %in% result$flipped_folds[["1"]])
  expect_false(1L %in% result$flipped_folds[["1"]])
  expect_false(3L %in% result$flipped_folds[["1"]])
  expect_length(result$flipped_folds[["2"]], 0)
})
