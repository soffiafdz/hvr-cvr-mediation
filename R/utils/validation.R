# =============================================================================
# Validation Utilities
# =============================================================================
# Functions for validating data and parameters
# =============================================================================

library(data.table)
source(here::here("R/utils/logging.R"))

#' Validate data.table has required columns
#'
#' @param data.dt data.table to validate
#' @param required_cols.v Character vector of required column names
#' @param data_name Name of data for error message
#' @export
validate_columns <- function(data.dt, required_cols.v, data_name = "data") {
  missing_cols.v <- setdiff(required_cols.v, names(data.dt))

  if (length(missing_cols.v) > 0) {
    msg <- sprintf(
      "Missing required columns in %s: %s",
      data_name,
      paste(missing_cols.v, collapse = ", ")
    )
    log_error(msg)
    stop(msg, call. = FALSE)
  }

  log_debug("All required columns present in %s", data_name)
  invisible(TRUE)
}

#' Validate data structure is not empty
#'
#' @param x Data structure to validate (data.table, data.frame, or list)
#' @param data_name Name of data for error message
#' @export
validate_not_empty <- function(x, data_name = "data") {
  if (is.data.frame(x)) {
    if (nrow(x) == 0) {
      msg <- sprintf("Dataset %s is empty", data_name)
      log_error(msg)
      stop(msg, call. = FALSE)
    }
    log_debug("Dataset %s has %d rows", data_name, nrow(x))
  } else if (is.list(x)) {
    if (length(x) == 0) {
      msg <- sprintf("List %s is empty", data_name)
      log_error(msg)
      stop(msg, call. = FALSE)
    }
    log_debug("List %s has %d elements", data_name, length(x))
  } else {
    if (length(x) == 0) {
      msg <- sprintf("Object %s is empty", data_name)
      log_error(msg)
      stop(msg, call. = FALSE)
    }
    log_debug("Object %s has length %d", data_name, length(x))
  }
  invisible(TRUE)
}

#' Validate numeric columns are in expected range
#'
#' @param data.dt data.table
#' @param col_name Column name
#' @param min_val Minimum expected value
#' @param max_val Maximum expected value
#' @param na_allowed Are NAs allowed
#' @export
validate_num_range <- function(
    data.dt, col_name, min_val = -Inf, max_val = Inf, na_allowed = TRUE) {
  if (!col_name %in% names(data.dt)) {
    msg <- sprintf("Column %s not found", col_name)
    log_error(msg)
    stop(msg, call. = FALSE)
  }

  values.v <- data.dt[[col_name]]

  # Check NAs
  n_na <- sum(is.na(values.v))
  if (n_na > 0 && !na_allowed) {
    msg <- sprintf("Column %s contains %d NAs (not allowed)", col_name, n_na)
    log_error(msg)
    stop(msg, call. = FALSE)
  }

  # Check range
  values_no_na.v <- values.v[!is.na(values.v)]
  if (length(values_no_na.v) > 0) {
    if (any(values_no_na.v < min_val) || any(values_no_na.v > max_val)) {
      actual_min <- min(values_no_na.v)
      actual_max <- max(values_no_na.v)
      msg <- sprintf(
        "Column %s out of range [%.2f, %.2f], actual [%.2f, %.2f]",
        col_name, min_val, max_val, actual_min, actual_max
      )
      log_error(msg)
      stop(msg, call. = FALSE)
    }
  }

  log_debug(
    "Column %s validated (range: [%.2f, %.2f])", col_name, min_val, max_val
  )
  invisible(TRUE)
}

#' Validate categorical column has expected levels
#'
#' @param data.dt data.table
#' @param col_name Column name
#' @param expected_levels.v Expected factor levels or unique values
#' @param allow_extra Allow extra levels not in expected
#' @export
validate_categorical <- function(
    data.dt, col_name, expected_levels.v, allow_extra = FALSE) {
  if (!col_name %in% names(data.dt)) {
    msg <- sprintf("Column %s not found", col_name)
    log_error(msg)
    stop(msg, call. = FALSE)
  }

  actual_levels.v <- unique(data.dt[[col_name]])
  actual_levels.v <- actual_levels.v[!is.na(actual_levels.v)]

  missing_levels.v <- setdiff(expected_levels.v, actual_levels.v)
  extra_levels.v <- setdiff(actual_levels.v, expected_levels.v)

  if (length(missing_levels.v) > 0) {
    log_warn(
      "Column %s missing expected levels: %s",
      col_name, paste(missing_levels.v, collapse = ", ")
    )
  }

  if (length(extra_levels.v) > 0 && !allow_extra) {
    msg <- sprintf(
      "Column %s has unexpected levels: %s",
      col_name, paste(extra_levels.v, collapse = ", ")
    )
    log_error(msg)
    stop(msg, call. = FALSE)
  }

  log_debug("Column %s validated (%d levels)", col_name,
            length(actual_levels.v))
  invisible(TRUE)
}

#' Validate age values are reasonable
#'
#' @param data.dt data.table
#' @param age_col Name of age column
#' @param min_age Minimum reasonable age
#' @param max_age Maximum reasonable age
#' @export
validate_age <- function(
    data.dt, age_col = "AGE", min_age = 18, max_age = 120, ...) {
  validate_num_range(
    data.dt, age_col,
    min_val = min_age, max_val = max_age, ...
  )
}

#' Validate IDs are unique
#'
#' @param data.dt data.table
#' @param id_cols.v Character vector of ID column names
#' @export
validate_unique_ids <- function(data.dt, id_cols.v) {
  if (!all(id_cols.v %in% names(data.dt))) {
    missing.v <- setdiff(id_cols.v, names(data.dt))
    msg <- sprintf("ID columns not found: %s",
                   paste(missing.v, collapse = ", "))
    log_error(msg)
    stop(msg, call. = FALSE)
  }

  n_before <- nrow(data.dt)
  n_unique <- nrow(unique(data.dt[, ..id_cols.v]))

  if (n_unique < n_before) {
    n_duplicates <- n_before - n_unique
    msg <- sprintf(
      "Found %d duplicate IDs in columns: %s",
      n_duplicates, paste(id_cols.v, collapse = ", ")
    )
    log_error(msg)
    stop(msg, call. = FALSE)
  }

  log_debug("IDs validated (%d unique)", n_unique)
  invisible(TRUE)
}

#' Validate required packages are installed
#'
#' @param packages.v Character vector of package names
#' @export
validate_packages <- function(packages.v) {
  missing.v <- character(0)

  for (pkg in packages.v) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing.v <- c(missing.v, pkg)
    }
  }

  if (length(missing.v) > 0) {
    install_str <- paste0("'", missing.v, "'", collapse = ", ")
    msg <- sprintf(
      "Missing packages: %s\nInstall with: install.packages(c(%s))",
      paste(missing.v, collapse = ", "),
      install_str
    )
    log_error(msg)
    stop(msg, call. = FALSE)
  }

  log_debug("All required packages installed (%d packages)",
            length(packages.v))
  invisible(TRUE)
}
