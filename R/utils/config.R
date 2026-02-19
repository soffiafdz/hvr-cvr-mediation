# =============================================================================
# Configuration Management
# =============================================================================
# Functions for loading and accessing pipeline configuration
# =============================================================================

library(yaml)
library(here)

.config <- NULL

#' Load pipeline configuration
#'
#' @param config_file Path to YAML configuration file
#' @return List containing configuration
#' @export
load_config <- function(config_file = here("config/pipeline_config.yaml")) {
  if (!file.exists(config_file)) {
    stop("Configuration file not found: ", config_file, call. = FALSE)
  }

  .config <<- read_yaml(config_file)
  invisible(.config)
}

#' Get configuration value
#'
#' @param ... Path to configuration value (e.g., "general", "random_seed")
#' @param default Default value if not found
#' @return Configuration value
#' @export
get_config <- function(..., default = NULL) {
  if (is.null(.config)) {
    load_config()
  }

  path.lst <- list(...)
  value <- .config

  for (key in path.lst) {
    if (is.null(value[[key]])) {
      if (!is.null(default)) {
        return(default)
      }
      stop("Configuration key not found: ", paste(path.lst, collapse = " -> "),
           call. = FALSE)
    }
    value <- value[[key]]
  }

  return(value)
}

#' Get data path from configuration
#'
#' @param ... Path keys (e.g., "processed", "covars_fst")
#' @return Full path to data file
#' @export
get_data_path <- function(...) {
  rel_path <- get_config("data", ...)
  here(rel_path)
}

#' Get script setting
#'
#' @param script_name Name of script (e.g., "gamlss")
#' @param setting Setting name (e.g., "redo_plots")
#' @param default Default value
#' @export
get_script_setting <- function(..., default = NULL) {
  # Supports nested paths
  # e.g.: get_script_setting("lgcm", "bootstrap", "enabled")
  path_elements.lst <- list(...)
  do.call(get_config, c(list("scripts"), path_elements.lst,
                        list(default = default)))
}

#' Get parameter value
#'
#' @param ... Path to parameter
#' @param default Default value
#' @export
get_parameter <- function(..., default = NULL) {
  get_config("parameters", ..., default = default)
}

#' Get random seed
#'
#' @export
get_seed <- function() {
  get_config("general", "random_seed", default = 1618)
}

#' Set random seed from configuration
#'
#' @export
set_seed <- function() {

  seed <- get_seed()
  set.seed(seed)
  invisible(seed)
}

#' Validate configuration has required fields
#'
#' Checks that all required configuration sections and keys exist.
#' Call after load_config() to ensure config is complete.
#'
#' @param config Configuration list (default: loaded config)
#' @return TRUE if valid, throws error otherwise
#' @export
validate_config <- function(config = NULL) {
  if (is.null(config)) {
    if (is.null(.config)) {
      load_config()
    }
    config <- .config
  }

  # Required top-level sections
  required_sections.v <- c("general", "data", "parameters", "scripts")
  missing_sections.v <- setdiff(required_sections.v, names(config))
  if (length(missing_sections.v) > 0) {
    stop("Missing required config sections: ",
         paste(missing_sections.v, collapse = ", "), call. = FALSE)
  }

  # Required general settings
  required_general.v <- c("random_seed")
  missing_general.v <- setdiff(required_general.v, names(config$general))
  if (length(missing_general.v) > 0) {
    stop("Missing required general settings: ",
         paste(missing_general.v, collapse = ", "), call. = FALSE)
  }

  # Required data paths
  required_data.v <- c("raw", "derivatives", "models")
  missing_data.v <- setdiff(required_data.v, names(config$data))
  if (length(missing_data.v) > 0) {
    stop("Missing required data sections: ",
         paste(missing_data.v, collapse = ", "), call. = FALSE)
  }

  invisible(TRUE)
}

#' Get LGCM timepoint configuration
#'
#' Returns valid timepoint labels (T0, T1, ...) based on lgcm.n_timepoints
#' config setting.
#'
#' @param config Configuration list (optional, loads from file if NULL)
#' @return Character vector of timepoint labels (e.g., c("T0", "T1", "T5"))
#' @export
get_valid_timepoints <- function(config = NULL) {
  if (is.null(config)) {
    n_timepoints <- get_config("lgcm", "n_timepoints", default = 6)
  } else {
    n_timepoints <- config$lgcm$n_timepoints
    if (is.null(n_timepoints)) n_timepoints <- 6
  }
  paste0("T", 0:(n_timepoints - 1))
}

#' Get number of LGCM timepoints
#'
#' @return Integer number of timepoints
#' @export
get_n_timepoints <- function() {
  get_config("lgcm", "n_timepoints", default = 6)
}

#' Get standard factor levels for categorical variables
#'
#' Returns the expected factor levels for common categorical variables
#' to ensure consistency across the pipeline.
#'
#' @param variable_name Name of the variable (e.g., "SEX", "ADJ_METHOD")
#' @return Character vector of expected levels in order
#' @export
get_factor_levels <- function(variable_name) {
  # Standard factor levels used throughout the pipeline
  factor_levels.lst <- list(
    SEX = c("Female", "Male"),
    ADJ_METHOD = c("NON", "PRP", "STX", "RES"),
    ADJ_LABEL = c("Unadjusted", "Proportions", "Stereotaxic", "Residuals"),
    SIDE = c("L", "R", "LR"),
    SEGMENTATION = c("CRS", "LNG"),
    COHORT = c("ALL", "MTCH", "SENS")
  )

  if (!variable_name %in% names(factor_levels.lst)) {
    warning("Unknown variable name: ", variable_name,
            ". No standard factor levels defined.", call. = FALSE)
    return(NULL)
  }

  factor_levels.lst[[variable_name]]
}
