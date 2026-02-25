#' Shared Latent Structure Specification
#'
#' @description
#' Specifies a shared latent spatial component between the ecological and
#' sampling processes. This is the core mechanism for identifying preferential
#' sampling.
#'
#' @param prior Prior on the coupling strength (β_share). Default uses a
#'   weakly informative prior centered at zero.
#' @param constrain_positive Logical. If TRUE, constrains the coupling to be
#'   positive (sampling increases where ecology is high). Default FALSE.
#'
#' @details
#' The shared component enters the model as:
#' \deqn{S_{eco}(x) = S_{shared}(x) + S_{eco-only}(x)}
#' \deqn{S_{samp}(x) = \beta_{share} \times S_{shared}(x) + S_{samp-only}(x)}
#'
#' When \eqn{\beta_{share} \neq 0}, the ecological and sampling processes
#' share spatial structure, indicating preferential sampling.
#'
#' @return A \code{spjam_latent} object of type "shared".
#'
#' @export
latent_shared <- function(prior = NULL, constrain_positive = FALSE) {
  structure(
    list(
      type = "shared",
      prior = prior,
      constrain_positive = constrain_positive
    ),
    class = c("spjam_latent_shared", "spjam_latent")
  )
}

#' Independent Latent Structure
#'
#' @description
#' Specifies independent spatial processes for ecology and sampling with
#' no shared structure. Use this as a comparison model to assess whether
#' preferential sampling correction matters.
#'
#' @return A \code{spjam_latent} object of type "independent".
#'
#' @export
latent_independent <- function() {
  structure(
    list(
      type = "independent"
    ),
    class = c("spjam_latent_independent", "spjam_latent")
  )
}

#' Correlated Latent Structure
#'
#' @description
#' Specifies correlated (but not identical) spatial processes for ecology
#' and sampling. Estimates the correlation coefficient between the two
#' latent fields.
#'
#' @param rho_prior Prior on the correlation coefficient. Default is
#'   uniform on (-1, 1).
#'
#' @return A \code{spjam_latent} object of type "correlated".
#'
#' @export
latent_correlated <- function(rho_prior = NULL) {
  structure(
    list(
      type = "correlated",
      rho_prior = rho_prior
    ),
    class = c("spjam_latent_correlated", "spjam_latent")
  )
}

#' @export
print.spjam_latent <- function(x, ...) {
  cat("spjam latent structure:", x$type, "\n")
  invisible(x)
}
