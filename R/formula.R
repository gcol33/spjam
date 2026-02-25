#' Create a spjam Formula Specification
#'
#' @description
#' Parses and validates a model formula for spjam.
#'
#' @param formula A formula specifying the ecological process model.
#' @param sampling Optional formula for the sampling process. If NULL,
#'   an intercept-only model is used.
#'
#' @return A \code{spjam_formula} object containing parsed formula components.
#'
#' @export
#'
#' @examples
#' # Ecological model only
#' spjam_formula(abundance ~ elevation + temperature)
#'
#' # With sampling process
#' spjam_formula(
#'   abundance ~ elevation + temperature,
#'   sampling = ~ road_distance + population
#' )
spjam_formula <- function(formula, sampling = NULL) {
  # Validate ecological formula
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula object", call. = FALSE)
  }

  # Parse sampling formula
  if (!is.null(sampling) && !inherits(sampling, "formula")) {
    stop("'sampling' must be a formula object or NULL", call. = FALSE)
  }

  structure(
    list(
      ecological = formula,
      sampling = sampling %||% ~ 1,
      terms_eco = NULL,
      terms_samp = NULL
    ),
    class = "spjam_formula"
  )
}

#' @export
print.spjam_formula <- function(x, ...) {
  cat("spjam formula specification\n")
  cat("  Ecological:", deparse(x$ecological), "\n")
  cat("  Sampling:  ", deparse(x$sampling), "\n")
  invisible(x)
}

# Null coalescing operator (canonical definition)
`%||%` <- function(x, y) if (is.null(x)) y else x
