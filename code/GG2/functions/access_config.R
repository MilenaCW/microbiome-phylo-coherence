# access_config.R
# ---------------
# Source a GG2 R config file and print a single value for a dot-separated key path.
#
# Usage:
#   Rscript access_config.R <config_R> <key.path>
#
# Examples:
#   Rscript access_config.R config/ocean.R output.directory
#   Rscript access_config.R config/ocean.R greengenes.backbone_qza
#
# Key paths use dots for nesting, e.g. output.directory, greengenes.tree_qza.

# Source repo setup (sets REPO_ROOT and working directory)
# Resolve this script's directory (works when run via Rscript)
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) != 1) stop("This script must be run with Rscript so --file= is available.")
script_path <- normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)
script_dir  <- dirname(script_path)
setup_path <- normalizePath(file.path(script_dir, "..", "..", "setup.R"),
                            winslash = "/", mustWork = TRUE)
source(setup_path)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: access_config.R <config_R> <key.path>", call. = FALSE)
}

cfg_path <- args[1]
if (!file.exists(cfg_path)) {
  stop("Config file not found: ", cfg_path, call. = FALSE)
}

cfg <- source(cfg_path, local = TRUE)$value
if (!is.list(cfg)) {
  stop("Config must evaluate to a list.", call. = FALSE)
}

key_path <- args[2]
keys <- strsplit(key_path, ".", fixed = TRUE)[[1L]]
val <- Reduce(function(x, n) x[[n]], keys, init = cfg)
if (is.null(val)) {
  stop("Config key not found or NULL: ", key_path, call. = FALSE)
}
if (!is.atomic(val) || length(val) != 1L) {
  stop("Key path must resolve to a single atomic value: ", key_path, call. = FALSE)
}
cat(as.character(val))
