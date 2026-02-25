#' Gaussian Process Spatial Specification
#'
#' @description
#' Specifies a Gaussian process spatial effect with Matérn covariance.
#'
#' @param range_prior Prior on the practical range. See \code{\link{prior_pc_range}}.
#' @param sigma_prior Prior on the marginal standard deviation.
#' @param nu Smoothness parameter (0.5, 1.5, or 2.5). Default 1.5.
#'
#' @return A \code{spjam_spatial} object of type "gp".
#'
#' @export
spatial_gp <- function(range_prior = NULL,
                       sigma_prior = NULL,
                       nu = 1.5) {
  structure(
    list(
      type = "gp",
      range_prior = range_prior,
      sigma_prior = sigma_prior,
      nu = nu
    ),
    class = c("spjam_spatial_gp", "spjam_spatial")
  )
}

#' SPDE Spatial Approximation
#'
#' @description
#' Specifies an SPDE (stochastic partial differential equation) approximation
#' to a Gaussian process for computational efficiency with large datasets.
#'
#' @param mesh A mesh object defining the spatial discretization.
#' @param range_prior Prior on the practical range.
#' @param sigma_prior Prior on the marginal standard deviation.
#'
#' @return A \code{spjam_spatial} object of type "spde".
#'
#' @export
spatial_spde <- function(mesh = NULL,
                         range_prior = NULL,
                         sigma_prior = NULL) {
  structure(
    list(
      type = "spde",
      mesh = mesh,
      range_prior = range_prior,
      sigma_prior = sigma_prior
    ),
    class = c("spjam_spatial_spde", "spjam_spatial")
  )
}

#' CAR Spatial Specification
#'
#' @description
#' Specifies a conditional autoregressive (CAR) spatial effect for areal data.
#'
#' @param adjacency Adjacency matrix or neighborhood list.
#' @param type Type of CAR: "icar" (intrinsic) or "proper".
#'
#' @return A \code{spjam_spatial} object of type "car".
#'
#' @export
spatial_car <- function(adjacency = NULL, type = c("icar", "proper")) {
  type <- match.arg(type)
  structure(
    list(
      type = "car",
      car_type = type,
      adjacency = adjacency
    ),
    class = c("spjam_spatial_car", "spjam_spatial")
  )
}

#' BYM2 Spatial Specification
#'
#' @description
#' Specifies a BYM2 (Besag-York-Mollié 2) spatial effect combining
#' structured and unstructured components with a mixing parameter.
#'
#' @param adjacency Adjacency matrix or neighborhood list.
#' @param phi_prior Prior on the mixing parameter (proportion structured).
#'
#' @return A \code{spjam_spatial} object of type "bym2".
#'
#' @export
spatial_bym2 <- function(adjacency = NULL, phi_prior = NULL) {
  structure(
    list(
      type = "bym2",
      adjacency = adjacency,
      phi_prior = phi_prior
    ),
    class = c("spjam_spatial_bym2", "spjam_spatial")
  )
}

#' @export
print.spjam_spatial <- function(x, ...) {
  cat("spjam spatial specification:", x$type, "\n")
  invisible(x)
}
