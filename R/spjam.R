#' Fit a Preferential Sampling Model
#'
#' @description
#' Jointly models an ecological spatial process and sampling intensity to
#' separate ecological signal from sampling-driven scatter.
#'
#' @param formula A formula or \code{spjam_formula} object specifying the
#'   ecological process model.
#' @param sampling A formula specifying the sampling process model. If NULL,
#'   uses an intercept-only LGCP.
#' @param data A data frame containing the response and covariates.
#' @param locations A matrix, sf object, or formula specifying observation
#'   locations. Required for continuous spatial models.
#' @param family A \code{spjam_family} object specifying the response
#'   distribution. Default is \code{spjam_family("poisson")}.
#' @param spatial Spatial structure specification. See \code{\link{spatial_gp}},
#'   \code{\link{spatial_spde}}, etc.
#' @param sampling_spatial Spatial structure for sampling process. If NULL,
#'   inherits from \code{spatial}.
#' @param shared Shared latent structure specification. See
#'   \code{\link{latent_shared}}.
#' @param priors Prior specifications. See \code{\link{spjam_priors}}.
#' @param backend Inference backend: \code{"hmc"} (default), \code{"laplace"},
#'   or \code{"auto"}.
#' @param chains Number of MCMC chains. Default 4.
#' @param iter Number of iterations per chain. Default 2000.
#' @param warmup Number of warmup iterations. Default \code{iter / 2}.
#' @param thin Thinning interval. Default 1.
#' @param cores Number of cores for parallel chains. Default 1.
#' @param seed Random seed for reproducibility.
#' @param init Initial values: \code{"default"}, \code{"random"}, or a list.
#' @param control A list of control parameters for the sampler.
#' @param ... Additional arguments passed to backend.
#'
#' @return An object of class \code{spjam_fit} containing:
#'   \item{draws}{Posterior samples}
#'   \item{formula}{Model specification}
#'   \item{sampling_formula}{Sampling process specification}
#'   \item{family}{Distribution family}
#'   \item{data}{Original data}
#'   \item{locations}{Observation locations}
#'   \item{backend}{Inference backend used}
#'   \item{diagnostics}{MCMC diagnostics}
#'   \item{decomposition}{Ecological/sampling/shared decomposition}
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' n <- 30
#' dat <- data.frame(
#'   abundance = rpois(n, exp(1.5)),
#'   elevation = rnorm(n),
#'   x = runif(n, 0, 5),
#'   y = runif(n, 0, 5)
#' )
#' fit <- spjam(
#'   abundance ~ elevation,
#'   data = dat,
#'   locations = ~ x + y,
#'   backend = "laplace",
#'   control = list(max_outer_iter = 5, n_samples = 50, verbose = FALSE)
#' )
#' print(fit)
#' }
spjam <- function(formula,
                  sampling = NULL,
                  data,
                  locations,
                  family = spjam_family("poisson"),
                  spatial = spatial_gp(),
                  sampling_spatial = NULL,
                  shared = latent_shared(),
                  priors = spjam_priors(),
                  backend = c("hmc", "laplace", "auto"),
                  chains = 4,
                  iter = 2000,
                  warmup = iter / 2,
                  thin = 1,
                  cores = 1,
                  seed = NULL,
                  init = "default",
                  control = list(),
                  ...) {

  # Match backend
  backend <- match.arg(backend)

  # Build spjam_formula if raw formula provided
  if (!inherits(formula, "spjam_formula")) {
    formula <- spjam_formula(formula, sampling = sampling)
  }

  # Default sampling spec if not provided
  sampling_spec <- if (inherits(sampling, "spjam_sampling")) {
    sampling
  } else {
    sampling_lgcp()
  }

  # Validate all inputs
  validate_spjam_inputs(formula, data, locations, family, spatial, shared,
                         sampling_spec)

  # Prepare model data (design matrices, locations, quadrature)
  model_data <- prepare_model_data(formula, data, locations, family,
                                    sampling_spec)

  # Prepare spatial structure (GP Cholesky, ICAR CSR, SPDE matrices)
  spatial_info <- prepare_spatial_structure(spatial, model_data)

  # For areal models, override n_spatial from adjacency matrix

  if (inherits(spatial, "spjam_spatial_car") ||
      inherits(spatial, "spjam_spatial_bym2")) {
    model_data$n_spatial <- spatial_info$n_spatial
    # Remap obs to spatial units (must be provided by user for areal data)
    # For now, assume 1:1 mapping (each observation is one area)
    if (length(model_data$obs_spatial_idx) != nrow(data)) {
      model_data$obs_spatial_idx <- seq(0L, nrow(data) - 1L)
    }
  }

  # Dispatch to backend
  if (backend == "hmc") {
    stop("HMC backend not yet implemented. Use backend = 'laplace'.",
         call. = FALSE)
  }

  if (backend == "auto") {
    backend <- "laplace"
  }

  n_samples <- control$n_samples %||% as.integer(chains * (iter - warmup) / thin)
  verbose <- control$verbose %||% TRUE

  laplace_result <- fit_laplace_joint(
    model_data = model_data,
    spatial_info = spatial_info,
    priors = priors,
    shared = shared,
    family = family,
    n_samples = n_samples,
    seed = seed,
    control = control,
    verbose = verbose
  )

  # Construct spjam_fit object
  new_spjam_fit(
    laplace_result = laplace_result,
    formula = formula,
    family = family,
    spatial = spatial,
    shared = shared,
    priors = priors,
    data = data,
    locations = locations,
    model_data = model_data,
    spatial_info = spatial_info,
    backend = backend,
    seed = seed
  )
}
