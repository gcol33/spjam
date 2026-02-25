#' Log-Gaussian Cox Process Sampling Specification
#'
#' @description
#' Specifies a log-Gaussian Cox process (LGCP) model for the sampling
#' intensity. This is the default and most flexible option for point-referenced
#' observation data.
#'
#' @param integration_points Number of integration points for the intensity
#'   integral. Default 1000.
#' @param boundary An sf polygon defining the study region boundary.
#'
#' @return A \code{spjam_sampling} object of type "lgcp".
#'
#' @export
sampling_lgcp <- function(integration_points = 1000, boundary = NULL) {
  structure(
    list(
      type = "lgcp",
      integration_points = integration_points,
      boundary = boundary
    ),
    class = c("spjam_sampling_lgcp", "spjam_sampling")
  )
}

#' Point Process Sampling Specification
#'
#' @description
#' Specifies an inhomogeneous Poisson point process model for sampling
#' locations. Simpler than LGCP but requires explicit intensity function.
#'
#' @param intensity_formula Formula for the log-intensity.
#'
#' @return A \code{spjam_sampling} object of type "point".
#'
#' @export
sampling_point <- function(intensity_formula = ~ 1) {
  structure(
    list(
      type = "point",
      intensity_formula = intensity_formula
    ),
    class = c("spjam_sampling_point", "spjam_sampling")
  )
}

#' Binomial Sampling Specification
#'
#' @description
#' Specifies a binomial model for grid-based sampling effort, where each
#' cell is either sampled or not sampled.
#'
#' @param grid An sf or terra object defining the sampling grid.
#' @param effort Optional known effort variable name in data.
#'
#' @return A \code{spjam_sampling} object of type "binomial".
#'
#' @export
sampling_binomial <- function(grid = NULL, effort = NULL) {
  structure(
    list(
      type = "binomial",
      grid = grid,
      effort = effort
    ),
    class = c("spjam_sampling_binomial", "spjam_sampling")
  )
}

#' @export
print.spjam_sampling <- function(x, ...) {
  cat("spjam sampling specification:", x$type, "\n")
  invisible(x)
}
