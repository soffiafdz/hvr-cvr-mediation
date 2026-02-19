# =============================================================================
# Logging Utilities
# =============================================================================
# Functions for consistent logging throughout the pipeline
# Author: Pipeline Refactoring
# Date: 2024
# =============================================================================

library(here)

# Initialize logger
.log_file <- NULL
.log_level <- "INFO"

#' Initialize logging system
#'
#' @param log_dir Directory for log files
#' @param log_level Minimum level to log (DEBUG, INFO, WARN, ERROR)
#' @param append Append to existing log file
#' @export
init_logger <- function(log_dir = here("logs"),
                       log_level = "INFO",
                       append = FALSE) {
  # Create log directory with error handling
  if (!dir.exists(log_dir)) {
    if (!dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)) {
      warning("Failed to create log directory: ", log_dir, call. = FALSE)
    }
  }

  # Set log file with timestamp
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  .log_file <<- here(log_dir, paste0("pipeline_", timestamp, ".log"))
  .log_level <<- log_level

  # Write header
  if (!append) {
    writeLines(
      c(
        paste0(rep("=", 80), collapse = ""),
        paste("Pipeline Execution Log"),
        paste("Started:", Sys.time()),
        paste("R Version:", R.version.string),
        paste0(rep("=", 80), collapse = ""),
        ""
      ),
      .log_file
    )
  }

  log_info("Logger initialized: %s", .log_file)
  invisible(.log_file)
}

#' Write log message
#'
#' @param level Log level (DEBUG, INFO, WARN, ERROR)
#' @param message Message to log
#' @param ... Arguments for sprintf formatting
#' @export
log_message <- function(level, message, ...) {
  if (is.null(.log_file)) {
    init_logger()
  }

  # Check if level should be logged
  LEVELS <- c("DEBUG" = 1, "INFO" = 2, "WARN" = 3, "ERROR" = 4)
  if (LEVELS[level] < LEVELS[.log_level]) {
    return(invisible(NULL))
  }

  # Format message
  if (length(list(...)) > 0) {
    message <- sprintf(message, ...)
  }

  # Replace glue-style placeholders
  message <- gsub("\\{([^}]+)\\}", "\\1", message)

  # Create log entry
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  entry <- sprintf("[%s] %s: %s", timestamp, level, message)

  # Write to file
  cat(entry, "\n", file = .log_file, append = TRUE)

  # Also print to console for INFO and above
  if (LEVELS[level] >= LEVELS["INFO"]) {
    message(entry)
  }

  invisible(NULL)
}

#' Log debug message
#' @export
log_debug <- function(message, ...) {
  log_message("DEBUG", message, ...)
}

#' Log info message
#' @export
log_info <- function(message, ...) {
  log_message("INFO", message, ...)
}

#' Log warning message
#' @export
log_warn <- function(message, ...) {
  log_message("WARN", message, ...)
}

#' Log error message
#' @export
log_error <- function(message, ...) {
  log_message("ERROR", message, ...)
}

#' Log script start
#' @export
log_script_start <- function(script_name) {
  log_info(paste0(rep("-", 80), collapse = ""))
  log_info("Starting script: %s", script_name)
  log_info(paste0(rep("-", 80), collapse = ""))
}

#' Log script end
#' @export
log_script_end <- function(script_name, success = TRUE) {
  status <- ifelse(success, "COMPLETED", "FAILED")
  log_info("Script %s: %s", status, script_name)
  log_info(paste0(rep("-", 80), collapse = ""))
}

#' Log section
#' @export
log_section <- function(section_name) {
  log_info("")
  log_info(">>> %s", section_name)
}

#' Time a code block and log duration
#' @export
log_time <- function(description, expr) {
  log_info("Starting: %s", description)
  start_time <- Sys.time()

  result <- tryCatch(
    expr,
    error = function(e) {
      log_error("Error in %s: %s", description, e$message)
      stop(e)
    }
  )

  duration <- difftime(Sys.time(), start_time, units = "secs")
  log_info("Completed: %s (%.2f seconds)", description, as.numeric(duration))

  invisible(result)
}
