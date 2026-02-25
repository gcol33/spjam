#' Check Model Identifiability
#'
#' @description
#' Assesses whether the data can distinguish ecological signal from
#' sampling-driven structure. Reports on posterior uncertainty for
#' the coupling parameter and variance decomposition.
#'
#' @param fit A \code{spjam_fit} object.
#' @param threshold Posterior width threshold for flagging weak
#'   identifiability. Default 0.5.
#'
#' @return A \code{spjam_identifiability} object with diagnostic information.
#'
#' @export
check_identifiability <- function(fit, threshold = 0.5) {
  if (!inherits(fit, "spjam_fit")) {
    stop("'fit' must be a spjam_fit object", call. = FALSE)
  }

  # Placeholder
  stop("check_identifiability() not yet implemented", call. = FALSE)
}

#' MCMC Diagnostics
#'
#' @description
#' Computes standard MCMC diagnostics including Rhat, effective sample size,
#' and divergence counts.
#'
#' @param fit A \code{spjam_fit} object.
#'
#' @return A \code{spjam_diagnostics} object.
#'
#' @export
mcmc_diagnostics <- function(fit) {
  if (!inherits(fit, "spjam_fit")) {
    stop("'fit' must be a spjam_fit object", call. = FALSE)
  }

  # Placeholder
  stop("mcmc_diagnostics() not yet implemented", call. = FALSE)
}

#' @export
print.spjam_diagnostics <- function(x, ...) {
  cat("spjam MCMC diagnostics\n\n")
  cat("Rhat:\n")
  cat("  Max:", sprintf("%.3f", x$rhat_max), "\n")
  cat("  Parameters > 1.01:", x$rhat_bad, "\n")
  cat("\nEffective sample size:\n")
  cat("  Min:", x$ess_min, "\n")
  cat("  Parameters < 400:", x$ess_low, "\n")
  if (x$divergences > 0) {
    cat("\nWARNING:", x$divergences, "divergent transitions\n")
  }
  invisible(x)
}
