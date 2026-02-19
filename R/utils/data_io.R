# =============================================================================
# Data I/O Utilities
# =============================================================================
# Functions for reading and writing data with consistent error handling
# =============================================================================

library(here)
library(data.table)

source(here("R/utils/logging.R"))

#' Check if required files exist
#'
#' @param files.v Character vector of file paths
#' @param stop_on_missing Throw error if files missing
#' @return Logical indicating all files exist
#' @export
check_files_exist <- function(files.v, stop_on_missing = TRUE) {
  missing.v <- files.v[!file.exists(files.v)]

  if (length(missing.v) > 0) {
    msg <- sprintf(
      "Missing required file(s):\n%s",
      paste("  -", missing.v, collapse = "\n")
    )

    if (stop_on_missing) {
      log_error(msg)
      stop(msg, call. = FALSE)
    } else {
      log_warn(msg)
      return(FALSE)
    }
  }

  log_debug("All required files exist (%d files)", length(files.v))
  return(TRUE)
}

#' Read RDS file with error handling
#'
#' @param file_path Path to RDS file
#' @param description Optional description for logging
#' @return Data object from RDS file
#' @export
read_rds_safe <- function(file_path, description = NULL) {
  if (!file.exists(file_path)) {
    log_error("RDS file not found: %s", file_path)
    stop("File not found: ", file_path, call. = FALSE)
  }

  desc <- if (!is.null(description)) description else basename(file_path)
  log_debug("Reading RDS: %s", desc)

  tryCatch(
    readRDS(file_path),
    error = function(e) {
      log_error("Failed to read RDS file %s: %s", file_path, e$message)
      stop(e)
    }
  )
}

#' Write RDS file with error handling
#'
#' @param object Object to save
#' @param file_path Path to save to
#' @param description Optional description for logging
#' @export
write_rds_safe <- function(object, file_path, description = NULL) {
  # Create directory if needed
  dir_path <- dirname(file_path)
  if (!dir.exists(dir_path)) {
    log_debug("Creating directory: %s", dir_path)
    dir.create(dir_path, recursive = TRUE)
  }

  desc <- if (!is.null(description)) description else basename(file_path)
  log_debug("Writing RDS: %s", desc)

  tryCatch(
    saveRDS(object, file_path),
    error = function(e) {
      log_error("Failed to write RDS file %s: %s", file_path, e$message)
      stop(e)
    }
  )

  log_info("Saved: %s", file_path)
  invisible(file_path)
}

#' Read CSV file with error handling
#'
#' @param file_path Path to CSV file
#' @param ... Additional arguments for fread
#' @param description Optional description for logging
#' @export
read_csv_safe <- function(file_path, ..., description = NULL) {
  if (!file.exists(file_path)) {
    log_error("CSV file not found: %s", file_path)
    stop("File not found: ", file_path, call. = FALSE)
  }

  desc <- if (!is.null(description)) description else basename(file_path)
  log_debug("Reading CSV: %s", desc)

  tryCatch(
    fread(file_path, ...),
    error = function(e) {
      log_error("Failed to read CSV file %s: %s", file_path, e$message)
      stop(e)
    }
  )
}

#' Check if output file needs to be regenerated
#'
#' @param output_path Path to output file
#' @param input_paths.v Vector of input file paths
#' @param force_regenerate Force regeneration regardless of modification times
#' @return Logical indicating if regeneration is needed
#' @export
needs_regeneration <- function(output_path,
                               input_paths.v = NULL,
                               force_regenerate = FALSE) {
  if (force_regenerate) {
    log_debug("Force regeneration: %s", basename(output_path))
    return(TRUE)
  }

  if (!file.exists(output_path)) {
    log_debug("Output does not exist: %s", basename(output_path))
    return(TRUE)
  }

  if (is.null(input_paths.v)) {
    log_debug("Output exists, no input tracking: %s", basename(output_path))
    return(FALSE)
  }

  # Check modification times
  output_mtime <- file.mtime(output_path)
  existing.v <- input_paths.v[file.exists(input_paths.v)]
  input_mtimes.v <- sapply(existing.v, file.mtime)

  if (any(input_mtimes.v > output_mtime)) {
    log_debug("Inputs newer than output: %s", basename(output_path))
    return(TRUE)
  }

  log_debug("Output up-to-date: %s", basename(output_path))
  return(FALSE)
}

#' Create directory if it doesn't exist
#'
#' @param dir_path Directory path
#' @export
ensure_directory <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    log_debug("Creating directory: %s", dir_path)
    dir.create(dir_path, recursive = TRUE)
  }
  invisible(dir_path)
}

#' Join path components
#'
#' Similar to file.path but preferred for consistency with here() usage.
#' Use this when you already have a base path from get_data_path() and
#' need to append a filename or subdirectory.
#'
#' @param ... Path components to join
#' @return Character string of joined path
#' @export
path_join <- function(...) {
  components.lst <- list(...)
  # Filter out NULL and empty strings
  is_empty.fn <- function(x) is.null(x) || x == ""
  components.lst <- components.lst[!sapply(components.lst, is_empty.fn)]
  do.call(file.path, components.lst)
}

# -----------------------------------------------------------
# Manuscript adapter helpers
# -----------------------------------------------------------

#' Extract per-CVR view from merged LME results
#'
#' @param res  Full LME results list (with res$stratified)
#' @param cvr  "FRS" or "CVR_mimic"
#' @return Copy of res with $stratified narrowed to one CVR
#' @export
extract_cvr_view.fn <- function(res, cvr) {
  stopifnot(
    !is.null(res),
    !is.null(res$stratified),
    cvr %in% names(res$stratified)
  )
  out <- res
  out$stratified <- res$stratified[[cvr]]
  out
}

#' Remap LGCM mediation keys for per-sex adapter
#'
#' Converts a_path/b_path to a/b and forwards
#' fit_indices to fit.
#'
#' @param r  Single mediation result element
#' @return Modified result with renamed fields
#' @export
remap_med.fn <- function(r) {
  stopifnot(!is.null(r))
  r$a <- r$a_path
  r$b <- r$b_path
  if (!is.null(r$fit_indices))
    r$fit <- r$fit_indices
  r
}

#' Guard: fail render if any output is NULL
#'
#' @param x   Object to check
#' @param label  Human label for error messages
#' @return x, invisibly, if non-NULL
#' @export
require_output <- function(x, label) {
  if (is.null(x)) {
    stop(sprintf(
      "Output '%s' is NULL. Fix upstream.",
      label
    ))
  }
  x
}
