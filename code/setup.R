# setup.R
# --------
# Shared setup for R scripts under code/. Source this file using a path
# relative to the calling script (e.g. from GG2/scripts: path ../../setup.R,
# from DADA: ../setup.R), or after ensuring the working directory is the repo
# root: source("code/setup.R").
#
# Sets working directory to the repo root and assigns REPO_ROOT in the
# calling environment. Detects the repo root dynamically — no edits needed.

.repo_root <- tryCatch(
  trimws(system("git rev-parse --show-toplevel", intern = TRUE, ignore.stderr = TRUE)),
  error = function(e) character(0)
)

if (length(.repo_root) == 0 || !dir.exists(.repo_root)) {
  .this_file <- tryCatch(
    normalizePath(sys.frames()[[length(sys.frames())]]$ofile),
    error = function(e) NULL
  )
  if (!is.null(.this_file)) {
    .repo_root <- dirname(dirname(.this_file))  # code/setup.R -> repo root
  } else {
    stop("Could not determine repo root. Run from within the repository or ensure git is in PATH.")
  }
}

setwd(.repo_root)
REPO_ROOT <<- getwd()
cat("[setup.R] Working directory set to:", REPO_ROOT, "\n", file = stderr())
rm(.repo_root)
