# Script to do the third step of DADA2 pipeline: learn the error rates

# Source repo setup (sets REPO_ROOT and working directory)
# Resolve this script's directory (works when run via Rscript)
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) != 1) stop("This script must be run with Rscript so --file= is available.")
script_path <- normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)
script_dir  <- dirname(script_path)
setup_path <- normalizePath(file.path(script_dir, "..", "setup.R"),
                            winslash = "/", mustWork = TRUE)
source(setup_path)

library(tidyverse)
library(dada2)
library(optparse)
source('./code/utility_functions.R')

# Define command-line options
option_list <- list(
  make_option(c("--config"), action = "store", default = NULL,
              help = "Path to an R config file that evaluates to a list (required)."),
  make_option(c("--flowcell"), action = "store", default = NULL,
              help = "Learn the error model for a given flowcell (Default: NULL, all flow cells)."),
  make_option(c("--verbose"), action = "store_true", default = TRUE,
              help = "Print verbose output (Default: TRUE).")
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

load_config <- function(config_path) {
  if (is.null(config_path) || !nzchar(config_path)) {
    stop("Missing required --config <path>", call. = FALSE)
  }
  if (!file.exists(config_path)) {
    stop(paste0("Config file not found: ", config_path), call. = FALSE)
  }
  cfg <- source(config_path, local = TRUE)$value
  if (!is.list(cfg)) {
    stop("Config must evaluate to an R list", call. = FALSE)
  }
  cfg
}

cfg <- load_config(opt$config)

output_directory <- cfg$output$directory
error_cfg <- cfg$error_model
verbose <- if (!is.null(cfg$verbose)) cfg$verbose else opt$verbose

if (is.null(output_directory)) {
  stop("Config missing required field: output$directory", call. = FALSE)
}

error_method <- if (!is.null(error_cfg$method)) tolower(error_cfg$method) else "standard"
if (!error_method %in% c("standard", "loess")) {
  stop("error_model$method must be one of: 'standard', 'loess'", call. = FALSE)
}

loessErrfun_mod <- function(trans) { # from https://github.com/benjjneb/dada2/issues/1307#issuecomment-821190155
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow = 0, ncol = length(qq))
  for (nti in c("A", "C", "G", "T")) {
    for (ntj in c("A", "C", "G", "T")) {
      if (nti != ntj) {
        errs <- trans[paste0(nti, "2", ntj), ]
        tot <- colSums(trans[paste0(nti, "2", c("A", "C", "G", "T")), ])
        rlogp <- log10((errs + 1) / tot)  # 1 psuedocount for each err, but if tot=0 will give NA
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q = qq, errs = errs, tot = tot, rlogp = rlogp)

        # Gulliem Salazar's solution
        # https://github.com/benjjneb/dada2/issues/938
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot), span = 2)

        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred) > maxrli] <- pred[[maxrli]]
        pred[seq_along(pred) < minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))

  # HACKY
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est > MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est < MIN_ERROR_RATE] <- MIN_ERROR_RATE

  # enforce monotonicity
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>%
    as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)

  # Expand the err matrix with the self-transition probs
  err <- rbind(1 - colSums(est[1:3, ]), est[1:3, ],
               est[4, ], 1 - colSums(est[4:6, ]), est[5:6, ],
               est[7:8, ], 1 - colSums(est[7:9, ]), est[9, ],
               est[10:12, ], 1 - colSums(est[10:12, ]))
  rownames(err) <- paste0(rep(c("A", "C", "G", "T"), each = 4), "2", c("A", "C", "G", "T"))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}

learn_error_model <- function(flowcell, verbose = TRUE) {
  start_time <- Sys.time()

  filtered_dir <- file.path(output_directory, flowcell, "02_filterAndTrim")
  if (!dir.exists(filtered_dir)) {
    stop(paste0("Filtered directory not found: ", filtered_dir), call. = FALSE)
  }

  # create directory for the error model
  error_dir <- file.path(output_directory, flowcell, "03_learnError")
  dir.create(error_dir, recursive = TRUE, showWarnings = FALSE)
  verbose_print(paste0("Created directory for error model: ", error_dir), verbose)

  # get the filtered files
  forward_files <- list.files(filtered_dir, pattern = "_F_filt.fastq.gz", full.names = TRUE)
  reverse_files <- list.files(filtered_dir, pattern = "_R_filt.fastq.gz", full.names = TRUE)
  verbose_print(paste0("Found ", length(forward_files), " forward files and ", length(reverse_files),
                       " reverse files in ", filtered_dir), verbose)
  if (length(forward_files) == 0 || length(reverse_files) == 0) {
    stop(paste0("No filtered files found in: ", filtered_dir), call. = FALSE)
  }
  names(forward_files) <- sub("_F_filt.fastq.gz", "", basename(forward_files))
  names(reverse_files) <- sub("_R_filt.fastq.gz", "", basename(reverse_files))

  learn_args <- list(randomize = TRUE, multithread = TRUE, verbose = verbose)
  if (error_method == "loess") {
    learn_args$errorEstimationFunction <- loessErrfun_mod
    verbose_print("Using modified loess error estimation function.", verbose)
  } else {
    verbose_print("Using standard DADA2 error estimation.", verbose)
  }

  errF <- do.call(learnErrors, c(list(forward_files), learn_args))
  errR <- do.call(learnErrors, c(list(reverse_files), learn_args))

  # save the error model
  saveRDS(errF, file.path(error_dir, paste0(flowcell, "_errF.rds")))
  saveRDS(errR, file.path(error_dir, paste0(flowcell, "_errR.rds")))

  # plot the error model
  p_errF <- plotErrors(errF, nominalQ = TRUE)
  ggplot2::ggsave(file.path(error_dir, "errF.jpeg"),
                  plot = p_errF, width = 7.5, height = 7.5, dpi = 300)
  p_errR <- plotErrors(errR, nominalQ = TRUE)
  ggplot2::ggsave(file.path(error_dir, "errR.jpeg"),
                  plot = p_errR, width = 7.5, height = 7.5, dpi = 300)

  end_time <- Sys.time()
  verbose_print(paste0("Error modeling complete. Time elapsed: ", hms_elapsed(start_time, end_time)), verbose)
}

if (!is.null(opt$flowcell)) {
  learn_error_model(opt$flowcell, verbose = verbose)
} else {
  verbose_print("No flow cell specified, learning error model for all flow cells...", verbose)
  flowcells <- list.dirs(output_directory, recursive = FALSE, full.names = FALSE)
  verbose_print("Flowcell directories found in output_directory:", verbose)
  cat(paste0("  ", flowcells, collapse = "\n"), "\n")
  for (flowcell in flowcells) {
    learn_error_model(flowcell, verbose = verbose)
  }
}