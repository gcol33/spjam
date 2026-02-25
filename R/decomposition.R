#' Decompose Spatial Structure into Components
#'
#' @description
#' Decomposes the fitted spatial surface into ecological-only, sampling-only,
#' and shared components. This is the primary diagnostic for understanding
#' how much spatial structure is attributable to preferential sampling.
#'
#' @param fit A \code{spjam_fit} object.
#' @param locations Optional matrix or sf object of locations at which to
#'   compute the decomposition. Default uses observation locations.
#' @param summary Logical. If TRUE (default), returns posterior summaries.
#'   If FALSE, returns full posterior draws.
#'
#' @return A \code{spjam_decomposition} object containing:
#'   \item{ecological}{Spatial structure attributable to ecology only}
#'   \item{sampling}{Spatial structure attributable to sampling only}
#'   \item{shared}{Shared structure (the confounding)}
#'   \item{variance_partition}{Proportion of variance in each component}
#'   \item{coupling}{Posterior for beta_share}
#'
#' @export
scatter_decomposition <- function(fit, locations = NULL, summary = TRUE) {
  if (!inherits(fit, "spjam_fit")) {
    stop("'fit' must be a spjam_fit object", call. = FALSE)
  }

  # Placeholder
 stop("scatter_decomposition() not yet implemented", call. = FALSE)
}

#' @export
print.spjam_decomposition <- function(x, ...) {
  cat("spjam scatter decomposition\n\n")
  cat("Variance partitioning:\n")
  cat("  Ecological only:", sprintf("%.1f%%", x$variance_partition$ecological * 100), "\n")
  cat("  Sampling only:  ", sprintf("%.1f%%", x$variance_partition$sampling * 100), "\n")
  cat("  Shared:         ", sprintf("%.1f%%", x$variance_partition$shared * 100), "\n")
  cat("\nCoupling strength (beta_share):\n")
  cat("  Mean:", sprintf("%.3f", x$coupling$mean), "\n")
  cat("  95% CI: [", sprintf("%.3f", x$coupling$lower),
      ",", sprintf("%.3f", x$coupling$upper), "]\n")
  invisible(x)
}

#' @export
plot.spjam_decomposition <- function(x, ...) {
  # Placeholder for spatial visualization
  stop("plot.spjam_decomposition() not yet implemented", call. = FALSE)
}
