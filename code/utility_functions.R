# This is a script to contain utility functions commonly used across all layers of the pipeline.

#' Print verbose messages
#' @param message Message to print
#' @param verbose Logical, whether to print the message
verbose_print <- function(message, verbose = FALSE) {
  if (verbose) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(paste0("[", timestamp, "] ", message, "\n"))
  }
}

hms_elapsed <- function(start, end) {
  dsec <- as.numeric(difftime(end, start, unit = "secs"))
  hours <- floor(dsec / 3600)
  minutes <- floor((dsec - 3600 * hours) / 60)
  seconds <- dsec - 3600*hours - 60*minutes
  paste0(
    sapply(c(hours, minutes, seconds), function(x) {
      formatC(x, width = 2, format = "d", flag = "0")
    }), collapse = ":")
}

clean_tax_name <- function(x) {
  if (length(x) > 1) {
    return(sapply(x, clean_tax_name, USE.NAMES = FALSE))
  }
  if (is.na(x) || x == "" || tolower(x) %in% c("undef", "na", "nan", "none")) {
    return("Unassigned")
  }
  as.character(x)
}