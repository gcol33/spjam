#' Specify Prior Distributions
#'
#' @description
#' Creates a prior specification object for spjam models.
#'
#' @param beta Prior for fixed effect coefficients.
#' @param intercept Prior for intercept terms.
#' @param range Prior for spatial range parameter.
#' @param sigma Prior for spatial marginal standard deviation.
#' @param coupling Prior for shared structure coupling (β_share).
#' @param phi Prior for overdispersion or mixing parameters.
#'
#' @return A \code{spjam_priors} object.
#'
#' @export
spjam_priors <- function(beta = prior_normal(0, 2.5),
                          intercept = prior_normal(0, 10),
                          range = NULL,
                          sigma = prior_half_normal(0, 1),
                          coupling = prior_normal(0, 1),
                          phi = prior_half_normal(0, 1)) {
  structure(
    list(
      beta = beta,
      intercept = intercept,
      range = range,
      sigma = sigma,
      coupling = coupling,
      phi = phi
    ),
    class = "spjam_priors"
  )
}

#' Normal Prior
#'
#' @param mean Mean of the normal distribution.
#' @param sd Standard deviation.
#'
#' @return A \code{spjam_prior} object.
#' @export
prior_normal <- function(mean = 0, sd = 1) {
  structure(
    list(type = "normal", mean = mean, sd = sd),
    class = "spjam_prior"
  )
}

#' Half-Normal Prior
#'
#' @param mean Location (typically 0).
#' @param sd Scale parameter.
#'
#' @return A \code{spjam_prior} object.
#' @export
prior_half_normal <- function(mean = 0, sd = 1) {
  structure(
    list(type = "half_normal", mean = mean, sd = sd),
    class = "spjam_prior"
  )
}

#' Half-Cauchy Prior
#'
#' @param location Location parameter (typically 0).
#' @param scale Scale parameter.
#'
#' @return A \code{spjam_prior} object.
#' @export
prior_half_cauchy <- function(location = 0, scale = 1) {
  structure(
    list(type = "half_cauchy", location = location, scale = scale),
    class = "spjam_prior"
  )
}

#' PC Prior for Spatial Range
#'
#' @description
#' Penalized complexity prior for the practical range of a spatial process.
#' Specifies P(range < U) = alpha.
#'
#' @param U Upper bound for range.
#' @param alpha Tail probability.
#'
#' @return A \code{spjam_prior} object.
#' @export
prior_pc_range <- function(U, alpha = 0.5) {
  structure(
    list(type = "pc_range", U = U, alpha = alpha),
    class = "spjam_prior"
  )
}

#' PC Prior for Marginal Standard Deviation
#'
#' @description
#' Penalized complexity prior for marginal standard deviation.
#' Specifies P(sigma > U) = alpha.
#'
#' @param U Upper bound.
#' @param alpha Tail probability.
#'
#' @return A \code{spjam_prior} object.
#' @export
prior_pc_sigma <- function(U, alpha = 0.5) {
  structure(
    list(type = "pc_sigma", U = U, alpha = alpha),
    class = "spjam_prior"
  )
}

#' @export
print.spjam_prior <- function(x, ...) {
  cat("Prior:", x$type, "\n")
  invisible(x)
}

#' @export
print.spjam_priors <- function(x, ...) {
  cat("spjam prior specification\n")
  for (nm in names(x)) {
    if (!is.null(x[[nm]])) {
      cat(" ", nm, ":", x[[nm]]$type, "\n")
    }
  }
  invisible(x)
}
