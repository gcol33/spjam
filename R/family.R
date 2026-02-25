#' Define Response Distribution Family
#'
#' @description
#' Specifies the likelihood family for the ecological response variable.
#'
#' @param family Character string specifying the family. One of:
#'   \code{"poisson"}, \code{"negbin"}, \code{"binomial"}, \code{"tweedie"},
#'   \code{"zi_poisson"}, \code{"zi_negbin"}.
#' @param link Link function. If NULL, uses the canonical link.
#'
#' @return A \code{scatr_family} object.
#'
#' @export
#'
#' @examples
#' scatr_family("poisson")
#' scatr_family("negbin")
#' scatr_family("binomial", link = "probit")
scatr_family <- function(family = c("poisson", "negbin", "binomial",
                                     "tweedie", "zi_poisson", "zi_negbin"),
                          link = NULL) {
  family <- match.arg(family)

  # Default links
  default_links <- list(
    poisson = "log",
    negbin = "log",
    binomial = "logit",
    tweedie = "log",
    zi_poisson = "log",
    zi_negbin = "log"
  )

  link <- link %||% default_links[[family]]

  structure(
    list(
      family = family,
      link = link,
      zero_inflated = grepl("^zi_", family)
    ),
    class = "scatr_family"
  )
}

#' @export
print.scatr_family <- function(x, ...) {
  cat("scatR family:", x$family, "(link =", x$link, ")\n")
  invisible(x)
}

# %||% defined in formula.R
