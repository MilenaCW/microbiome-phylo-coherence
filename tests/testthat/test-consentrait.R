library(testthat)
library(ape)
library(dplyr)
source(file.path(Sys.getenv("REPO_ROOT", "."), "code/CCA/functions/consentrait_signed.R"))

# Balanced 4-tip tree: ((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2)
TREE <- read.tree(text = "((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);")

test_that("consentrait_signed detects two coherent clades for perfectly segregated traits", {
  traits <- c(A = 1, B = 1, C = -1, D = -1)
  result <- consentrait_signed(TREE, traits, frac_consensus = 0.9, n_shuffles = 0)
  expect_equal(sum(result$clades$type == "clade"), 2L)
})

test_that("consentrait_signed returns finite positive tau_D for conserved traits", {
  traits <- c(A = 1, B = 1, C = -1, D = -1)
  result <- consentrait_signed(TREE, traits, frac_consensus = 0.9, n_shuffles = 0)
  expect_true(is.finite(result$tau_D))
  expect_gt(result$tau_D, 0)
})

test_that("consentrait_signed returns NA tau_D when no coherent clades form", {
  # Alternating signs: no clade passes frac_consensus = 0.9
  traits <- c(A = 1, B = -1, C = 1, D = -1)
  result <- consentrait_signed(TREE, traits, frac_consensus = 0.9, n_shuffles = 0)
  expect_true(is.na(result$tau_D))
})

test_that("consentrait_signed null_tau_D has correct length", {
  traits <- c(A = 1, B = 1, C = -1, D = -1)
  result <- consentrait_signed(TREE, traits, n_shuffles = 25, seed = 42)
  expect_length(result$null_tau_D, 25L)
})

test_that("consentrait_signed clade directions are +1 and -1 for segregated traits", {
  traits <- c(A = 1, B = 1, C = -1, D = -1)
  result <- consentrait_signed(TREE, traits, n_shuffles = 0)
  clade_dirs <- sort(unique(result$clades$direction[result$clades$type == "clade"]))
  expect_equal(clade_dirs, c(-1, 1))
})

test_that("consentrait_signed errors on non-binary trait values", {
  traits_bad <- c(A = 0.5, B = 1, C = -1, D = -1)
  expect_error(consentrait_signed(TREE, traits_bad), "must contain only")
})

test_that("consentrait_signed errors when no tips match between tree and traits", {
  traits_bad <- c(X = 1, Y = -1, Z = 1, W = -1)
  expect_error(consentrait_signed(TREE, traits_bad), "No tips in common")
})

test_that("consentrait_signed tau_D equals analytic value for balanced 4-tip tree", {
  traits <- c(A = 1, B = 1, C = -1, D = -1)
  result <- consentrait_signed(TREE, traits, frac_consensus = 0.9, n_shuffles = 0)
  # AB clade: node depth=0.2, mean tip depth=0.3 → clade depth=0.1
  # CD clade: same geometry → clade depth=0.1
  # tau_D = mean(0.1, 0.1) = 0.1
  expect_equal(result$tau_D, 0.1)
})

test_that("consentrait_signed tau_D increases with within-clade branch depth", {
  traits       <- c(A = 1, B = 1, C = -1, D = -1)
  tree_shallow <- ape::read.tree(text = "((A:0.1,B:0.1):0.8,(C:0.1,D:0.1):0.8);")
  tree_deep    <- ape::read.tree(text = "((A:0.8,B:0.8):0.1,(C:0.8,D:0.8):0.1);")
  tau_shallow  <- consentrait_signed(tree_shallow, traits, n_shuffles = 0)$tau_D
  tau_deep     <- consentrait_signed(tree_deep,    traits, n_shuffles = 0)$tau_D
  expect_lt(tau_shallow, tau_deep)
})

test_that("consentrait_signed handles partial tip/trait overlap with a warning", {
  tree5  <- ape::read.tree(
    text = "((A:0.1,B:0.1):0.2,(C:0.1,(D:0.1,E:0.1):0.05):0.2);"
  )
  traits <- c(A = 1, B = 1, C = -1, D = -1)  # E absent from traits
  expect_warning(
    result <- consentrait_signed(tree5, traits, n_shuffles = 0),
    "tips in tree not in trait_values"
  )
  expect_true(is.finite(result$tau_D))
})

test_that("consentrait_signed weight_clades=TRUE shifts tau_D toward larger clades", {
  tree6      <- ape::read.tree(
    text = "((A:0.1,B:0.1):0.4,((C:0.1,D:0.1):0.1,(E:0.1,F:0.1):0.1):0.3);"
  )
  traits     <- c(A = 1, B = 1, C = -1, D = -1, E = -1, F = -1)
  unweighted <- consentrait_signed(tree6, traits, weight_clades = FALSE, n_shuffles = 0)$tau_D
  weighted   <- consentrait_signed(tree6, traits, weight_clades = TRUE,  n_shuffles = 0)$tau_D
  # AB clade: size=2, depth=0.1; CDEF clade: size=4, depth=0.2
  # Unweighted: (0.1+0.2)/2=0.15; Weighted: (2*0.1+4*0.2)/6≈0.167
  expect_gt(weighted, unweighted)
})

test_that("consentrait_signed singleton_depth_frac > 0 produces finite tau_D when no clades form", {
  traits    <- c(A = 1, B = -1, C = 1, D = -1)  # alternating: no clades at frac=0.9
  res_frac0 <- consentrait_signed(TREE, traits, singleton_depth_frac = 0, n_shuffles = 0)
  res_frac1 <- consentrait_signed(TREE, traits, singleton_depth_frac = 1, n_shuffles = 0)
  expect_true(is.na(res_frac0$tau_D))
  expect_true(is.finite(res_frac1$tau_D))
})

test_that("consentrait_signed null distribution contains NA values for balanced 4-tip tree", {
  # Only 2 of C(4,2)=6 tip assignments produce coherent clades; the rest give tau_D=NA
  traits <- c(A = 1, B = 1, C = -1, D = -1)
  result <- consentrait_signed(TREE, traits, n_shuffles = 100, seed = 42)
  expect_true(any(is.na(result$null_tau_D)))
})
