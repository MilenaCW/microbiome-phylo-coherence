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

# ---- coarsen_composition ----

test_that("coarsen_composition sums relative abundance within each tax group and sample", {
  composition_df <- data.frame(
    sample_id = c("s1", "s2"),
    otu1 = c(0.1, 0.2),
    otu2 = c(0.3, 0.1),
    otu3 = c(0.6, 0.7)
  )
  taxonomy <- data.frame(
    Feature_ID = c("otu1", "otu2", "otu3"),
    Genus = c("A", "A", "B")
  )
  result <- coarsen_composition(composition_df, taxonomy, "Genus")
  result <- result[order(result$sample_id), ]
  expect_equal(sort(setdiff(names(result), "sample_id")), c("A", "B"))
  expect_equal(result$A[result$sample_id == "s1"], 0.4)
  expect_equal(result$B[result$sample_id == "s1"], 0.6)
  expect_equal(result$A[result$sample_id == "s2"], 0.3)
  expect_equal(result$B[result$sample_id == "s2"], 0.7)
})

test_that("coarsen_composition fills absent group/sample combinations with 0", {
  composition_df <- data.frame(
    sample_id = c("s1", "s2"),
    otu1 = c(1, 0),
    otu2 = c(0, 1)
  )
  taxonomy <- data.frame(
    Feature_ID = c("otu1", "otu2"),
    Genus = c("A", "B")
  )
  result <- coarsen_composition(composition_df, taxonomy, "Genus")
  result <- result[order(result$sample_id), ]
  expect_equal(result$A, c(1, 0))
  expect_equal(result$B, c(0, 1))
})

test_that("coarsen_composition errors when a feature is missing from taxonomy", {
  composition_df <- data.frame(sample_id = "s1", otu1 = 1, otu_unknown = 0)
  taxonomy <- data.frame(Feature_ID = "otu1", Genus = "A")
  expect_error(coarsen_composition(composition_df, taxonomy, "Genus"))
})

# ---- shuffle_taxonomy_labels ----

test_that("shuffle_taxonomy_labels preserves group sizes at every taxonomic level", {
  taxonomy <- data.frame(
    Feature_ID = paste0("otu", 1:6),
    Genus = c("A", "A", "A", "B", "B", "C"),
    Family = c("X", "X", "X", "X", "X", "Y"),
    stringsAsFactors = FALSE
  )
  shuffled <- shuffle_taxonomy_labels(taxonomy, seed = 1)
  expect_equal(sort(table(shuffled$Genus)), sort(table(taxonomy$Genus)))
  expect_equal(sort(table(shuffled$Family)), sort(table(taxonomy$Family)))
  expect_equal(sort(shuffled$Feature_ID), sort(taxonomy$Feature_ID))
})

test_that("shuffle_taxonomy_labels actually permutes the assignment", {
  taxonomy <- data.frame(
    Feature_ID = paste0("otu", 1:20),
    Genus = rep(c("A", "B"), 10),
    stringsAsFactors = FALSE
  )
  shuffled <- shuffle_taxonomy_labels(taxonomy, seed = 1)
  expect_false(identical(shuffled$Feature_ID, taxonomy$Feature_ID))
})

test_that("shuffle_taxonomy_labels is reproducible given a seed", {
  taxonomy <- data.frame(
    Feature_ID = paste0("otu", 1:10),
    Genus = rep(c("A", "B"), 5),
    stringsAsFactors = FALSE
  )
  s1 <- shuffle_taxonomy_labels(taxonomy, seed = 42)
  s2 <- shuffle_taxonomy_labels(taxonomy, seed = 42)
  expect_equal(s1$Feature_ID, s2$Feature_ID)
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
