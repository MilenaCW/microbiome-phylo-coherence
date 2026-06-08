library(testthat)
source(file.path(Sys.getenv("REPO_ROOT", "."), "code/utility_functions.R"))

test_that("verbose_print emits nothing when verbose=FALSE", {
  expect_silent(verbose_print("hello", verbose = FALSE))
})

test_that("verbose_print emits message when verbose=TRUE", {
  expect_output(verbose_print("hello", verbose = TRUE), "hello")
})

test_that("hms_elapsed formats H:MM:SS correctly", {
  t1 <- as.POSIXct("2025-01-01 00:00:00", tz = "UTC")
  t2 <- as.POSIXct("2025-01-01 01:02:33", tz = "UTC")
  expect_equal(hms_elapsed(t1, t2), "01:02:33")
})

test_that("hms_elapsed handles zero elapsed time", {
  t <- as.POSIXct("2025-01-01 00:00:00", tz = "UTC")
  expect_equal(hms_elapsed(t, t), "00:00:00")
})

test_that("clean_tax_name returns Unassigned for NA", {
  expect_equal(clean_tax_name(NA_character_), "Unassigned")
})

test_that("clean_tax_name returns Unassigned for empty string", {
  expect_equal(clean_tax_name(""), "Unassigned")
})

test_that("clean_tax_name returns Unassigned for 'undef'", {
  expect_equal(clean_tax_name("undef"), "Unassigned")
})

test_that("clean_tax_name passes through normal names", {
  expect_equal(clean_tax_name("Bacteroidetes"), "Bacteroidetes")
})

test_that("clean_tax_name handles vector input", {
  input  <- c("Firmicutes", NA_character_, "undef", "", "Proteobacteria")
  expect_equal(
    clean_tax_name(input),
    c("Firmicutes", "Unassigned", "Unassigned", "Unassigned", "Proteobacteria")
  )
})

test_that("clean_tax_name returns Unassigned for 'nan', 'none', and string 'NA'", {
  expect_equal(clean_tax_name("nan"),  "Unassigned")
  expect_equal(clean_tax_name("none"), "Unassigned")
  expect_equal(clean_tax_name("NA"),   "Unassigned")  # string "NA", not R's NA_character_
})
