# Run from repo root: REPO_ROOT=$(pwd) Rscript tests/run_tests.R
library(testthat)
test_dir("tests/testthat", reporter = "progress")
