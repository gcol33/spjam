#' Prepare model data for C++ backend
#'
#' Parses formulas into design matrices, extracts locations, and generates
#' quadrature points for the LGCP sampling process.
#'
#' @keywords internal
prepare_model_data <- function(formula, data, locations, family,
                                sampling_spec) {
  # ----------------------------------------------------------------
  # 1. Extract locations
  # ----------------------------------------------------------------
  if (inherits(locations, "formula")) {
    loc_vars <- all.vars(locations)
    loc_mat <- as.matrix(data[, loc_vars, drop = FALSE])
  } else {
    loc_mat <- as.matrix(locations)
  }
  colnames(loc_mat) <- c("x", "y")

  # ----------------------------------------------------------------
  # 2. Build ecological design matrix
  # ----------------------------------------------------------------
  eco_terms <- terms(formula$ecological, data = data)
  eco_mf <- model.frame(eco_terms, data = data, na.action = na.pass)
  X_eco <- model.matrix(eco_terms, data = eco_mf)

  # Response variable
  response_name <- all.vars(formula$ecological)[1]
  y <- data[[response_name]]
  if (!is.numeric(y)) {
    stop(sprintf("response variable '%s' must be numeric", response_name),
         call. = FALSE)
  }
  y_int <- as.integer(y)

  # n_trials for binomial
  n_trials <- rep(1L, length(y_int))
  if (family$family == "binomial") {
    # Check for cbind(success, failure) or assume n=1
    resp_expr <- formula$ecological[[2]]
    if (is.call(resp_expr) && deparse(resp_expr[[1]]) == "cbind") {
      successes <- eval(resp_expr[[2]], data)
      failures <- eval(resp_expr[[3]], data)
      y_int <- as.integer(successes)
      n_trials <- as.integer(successes + failures)
    }
  }

  # ----------------------------------------------------------------
  # 3. Build sampling design matrix
  # ----------------------------------------------------------------
  samp_terms <- terms(formula$sampling, data = data)
  samp_mf <- model.frame(samp_terms, data = data, na.action = na.pass)
  Z_samp <- model.matrix(samp_terms, data = samp_mf)

  # ----------------------------------------------------------------
  # 4. Generate quadrature points for LGCP
  # ----------------------------------------------------------------
  quad <- generate_quadrature(loc_mat, sampling_spec, Z_samp)

  # ----------------------------------------------------------------
  # 5. Observation-to-spatial mapping
  # ----------------------------------------------------------------
  # For continuous spatial models (GP, SPDE), each observation IS a spatial unit
  # For areal models (CAR, BYM2), obs_spatial_idx maps obs to areas
  # For now: unique locations define spatial units
  unique_locs <- unique(loc_mat)
  n_spatial <- nrow(unique_locs)

  # Map each observation to its spatial unit (0-based for C++)
  obs_spatial_idx <- match(
    paste(loc_mat[, 1], loc_mat[, 2]),
    paste(unique_locs[, 1], unique_locs[, 2])
  ) - 1L

  list(
    y = y_int,
    n_trials = n_trials,
    X_eco = X_eco,
    Z_samp = Z_samp,
    Z_quad = quad$Z_quad,
    quad_weights = quad$weights,
    quad_spatial_idx = quad$spatial_idx,
    locations = loc_mat,
    unique_locations = unique_locs,
    n_spatial = n_spatial,
    obs_spatial_idx = obs_spatial_idx,
    p_eco = ncol(X_eco),
    p_samp = ncol(Z_samp)
  )
}

#' Generate quadrature points for LGCP integral
#' @keywords internal
generate_quadrature <- function(loc_mat, sampling_spec, Z_samp) {
  n_quad <- sampling_spec$integration_points %||% 1000L

  # Bounding box with 10% buffer
  xrange <- range(loc_mat[, 1])
  yrange <- range(loc_mat[, 2])
  dx <- diff(xrange)
  dy <- diff(yrange)
  xmin <- xrange[1] - 0.1 * dx
  xmax <- xrange[2] + 0.1 * dx
  ymin <- yrange[1] - 0.1 * dy
  ymax <- yrange[2] + 0.1 * dy

  # Regular grid
  n_side <- ceiling(sqrt(n_quad))
  x_grid <- seq(xmin, xmax, length.out = n_side)
  y_grid <- seq(ymin, ymax, length.out = n_side)
  quad_locs <- expand.grid(x = x_grid, y = y_grid)
  quad_locs <- as.matrix(quad_locs)

  # Equal weights (area per quadrature point)
  area <- (xmax - xmin) * (ymax - ymin)
  n_actual <- nrow(quad_locs)
  weights <- rep(area / n_actual, n_actual)

  # Build design matrix at quadrature points
  # For intercept-only sampling formula, Z_quad is just a column of 1s
  p_samp <- ncol(Z_samp)
  Z_quad <- matrix(0.0, nrow = n_actual, ncol = p_samp)
  Z_quad[, 1] <- 1.0  # Intercept always present

  # For non-intercept covariates, nearest-neighbor interpolation
  if (p_samp > 1) {
    for (q in seq_len(n_actual)) {
      dists <- (loc_mat[, 1] - quad_locs[q, 1])^2 +
               (loc_mat[, 2] - quad_locs[q, 2])^2
      nearest <- which.min(dists)
      Z_quad[q, ] <- Z_samp[nearest, ]
    }
  }

  # Map quadrature points to nearest spatial unit (0-based)
  unique_locs <- unique(loc_mat)
  quad_spatial_idx <- integer(n_actual)
  for (q in seq_len(n_actual)) {
    dists <- (unique_locs[, 1] - quad_locs[q, 1])^2 +
             (unique_locs[, 2] - quad_locs[q, 2])^2
    quad_spatial_idx[q] <- which.min(dists) - 1L
  }

  list(
    locations = quad_locs,
    Z_quad = Z_quad,
    weights = weights,
    spatial_idx = quad_spatial_idx
  )
}
