#' Fit joint model via Laplace approximation
#'
#' Optimizes hyperparameters via R's optim() (outer loop),
#' with C++ Newton-Raphson for the latent field (inner loop).
#'
#' @keywords internal
fit_laplace_joint <- function(model_data, spatial_info, priors, shared,
                               family, n_samples = 1000L, seed = NULL,
                               control = list(), verbose = TRUE) {
  # ----------------------------------------------------------------
  # 1. Set up hyperparameter vector
  # ----------------------------------------------------------------
  spatial_type <- spatial_info$spatial_type
  has_shared <- (shared$type == "shared")
  latent_type <- if (has_shared) 0L else 1L  # LATENT_SHARED or LATENT_INDEPENDENT

  # Hyperparameters on internal scale (log/logit transforms)
  # [log_range_eco, log_sigma_eco, log_range_samp, log_sigma_samp,
  #  log_range_shared, log_sigma_shared, beta_share, log_phi]
  theta_names <- character(0)
  theta_init <- numeric(0)
  theta_lower <- numeric(0)
  theta_upper <- numeric(0)

  # For GP: range + sigma per field
  # For ICAR: just log_tau per field
  # For BYM2: log_sigma + logit_rho per field
  # For SPDE: log_range + log_sigma per field

  if (spatial_type == 0L) {  # GP
    # Ecological field
    init_range <- median(spatial_info$dist_mat[spatial_info$dist_mat > 0]) * 0.5
    theta_names <- c(theta_names, "log_range_eco", "log_sigma_eco")
    theta_init <- c(theta_init, log(init_range), log(1.0))
    theta_lower <- c(theta_lower, log(1e-4), log(1e-4))
    theta_upper <- c(theta_upper, log(max(spatial_info$dist_mat) * 2), log(100))

    # Sampling field
    theta_names <- c(theta_names, "log_range_samp", "log_sigma_samp")
    theta_init <- c(theta_init, log(init_range), log(1.0))
    theta_lower <- c(theta_lower, log(1e-4), log(1e-4))
    theta_upper <- c(theta_upper, log(max(spatial_info$dist_mat) * 2), log(100))

    if (has_shared) {
      theta_names <- c(theta_names, "log_range_shared", "log_sigma_shared")
      theta_init <- c(theta_init, log(init_range), log(1.0))
      theta_lower <- c(theta_lower, log(1e-4), log(1e-4))
      theta_upper <- c(theta_upper, log(max(spatial_info$dist_mat) * 2), log(100))
    }

  } else if (spatial_type == 1L || spatial_type == 2L) {  # ICAR or BYM2
    theta_names <- c(theta_names, "log_tau_eco", "log_tau_samp")
    theta_init <- c(theta_init, log(1.0), log(1.0))
    theta_lower <- c(theta_lower, log(1e-4), log(1e-4))
    theta_upper <- c(theta_upper, log(1e4), log(1e4))

    if (has_shared) {
      theta_names <- c(theta_names, "log_tau_shared")
      theta_init <- c(theta_init, log(1.0))
      theta_lower <- c(theta_lower, log(1e-4))
      theta_upper <- c(theta_upper, log(1e4))
    }

  } else if (spatial_type == 3L) {  # SPDE
    init_range <- median(spatial_info$dist_mat[spatial_info$dist_mat > 0]) * 0.5
    theta_names <- c(theta_names, "log_range_eco", "log_sigma_eco")
    theta_init <- c(theta_init, log(init_range), log(1.0))
    theta_lower <- c(theta_lower, log(1e-4), log(1e-4))
    theta_upper <- c(theta_upper, log(1e4), log(100))

    theta_names <- c(theta_names, "log_range_samp", "log_sigma_samp")
    theta_init <- c(theta_init, log(init_range), log(1.0))
    theta_lower <- c(theta_lower, log(1e-4), log(1e-4))
    theta_upper <- c(theta_upper, log(1e4), log(100))

    if (has_shared) {
      theta_names <- c(theta_names, "log_range_shared", "log_sigma_shared")
      theta_init <- c(theta_init, log(init_range), log(1.0))
      theta_lower <- c(theta_lower, log(1e-4), log(1e-4))
      theta_upper <- c(theta_upper, log(1e4), log(100))
    }
  }

  # Coupling parameter beta_share
  if (has_shared) {
    theta_names <- c(theta_names, "beta_share")
    theta_init <- c(theta_init, 0.0)
    theta_lower <- c(theta_lower, -10)
    theta_upper <- c(theta_upper, 10)
  }

  # Overdispersion (negbin)
  if (family$family == "negbin") {
    theta_names <- c(theta_names, "log_phi")
    theta_init <- c(theta_init, log(1.0))
    theta_lower <- c(theta_lower, log(0.01))
    theta_upper <- c(theta_upper, log(1000))
  }

  names(theta_init) <- theta_names

  # ----------------------------------------------------------------
  # 2. Objective function (negative Laplace marginal likelihood)
  # ----------------------------------------------------------------
  tau_beta <- control$tau_beta %||% 1e-4
  max_newton <- control$max_newton %||% 50L
  newton_tol <- control$newton_tol %||% 1e-6

  eval_count <- 0L

  obj_fn <- function(theta) {
    eval_count <<- eval_count + 1L

    # Unpack hyperparameters
    hp <- unpack_hyperparams(theta, theta_names, spatial_type, has_shared,
                              family$family)

    # Build spatial precision/covariance for current hyperparameters
    spatial_args <- build_spatial_args(hp, spatial_info, model_data$n_spatial,
                                       has_shared)

    # Call C++ inner loop
    result <- cpp_joint_laplace_fit(
      y = model_data$y,
      n_trials = model_data$n_trials,
      X_eco = model_data$X_eco,
      family = family$family,
      phi_overdispersion = hp$phi,
      Z_samp = model_data$Z_samp,
      Z_quad = model_data$Z_quad,
      quad_weights = model_data$quad_weights,
      quad_spatial_idx = model_data$quad_spatial_idx,
      spatial_type = spatial_type,
      obs_spatial_idx = model_data$obs_spatial_idx,
      gp_L_shared_flat = spatial_args$gp_L_shared,
      gp_L_eco_flat = spatial_args$gp_L_eco,
      gp_L_samp_flat = spatial_args$gp_L_samp,
      adj_row_ptr = spatial_info$adj_row_ptr,
      adj_col_idx = spatial_info$adj_col_idx,
      n_neighbors_vec = spatial_info$n_neighbors,
      tau_shared = hp$tau_shared,
      tau_eco = hp$tau_eco,
      tau_samp = hp$tau_samp,
      spde_Q_row_ptr = spatial_args$spde_Q_row_ptr,
      spde_Q_col_idx = spatial_args$spde_Q_col_idx,
      spde_Q_shared_vals = spatial_args$spde_Q_shared_vals,
      spde_Q_eco_vals = spatial_args$spde_Q_eco_vals,
      spde_Q_samp_vals = spatial_args$spde_Q_samp_vals,
      latent_type = latent_type,
      beta_share = hp$beta_share,
      p_eco = model_data$p_eco,
      p_samp = model_data$p_samp,
      n_spatial = model_data$n_spatial,
      tau_beta = tau_beta,
      max_iter = max_newton,
      tol = newton_tol
    )

    if (!result$converged && verbose) {
      message(sprintf("  [eval %d] Newton did not converge (%d iters)",
                      eval_count, result$n_iter))
    }

    -result$log_marginal  # Minimize negative log-marginal
  }

  # ----------------------------------------------------------------
  # 3. Optimize hyperparameters
  # ----------------------------------------------------------------
  if (verbose) message("Optimizing hyperparameters...")

  opt <- stats::optim(
    par = theta_init,
    fn = obj_fn,
    method = "L-BFGS-B",
    lower = theta_lower,
    upper = theta_upper,
    control = list(
      maxit = control$max_outer_iter %||% 100L,
      trace = if (verbose) 1L else 0L
    )
  )

  if (verbose) {
    message(sprintf("Optimization %s after %d evaluations (value = %.2f)",
                    if (opt$convergence == 0) "converged" else "stopped",
                    eval_count, opt$value))
  }

  # ----------------------------------------------------------------
  # 4. Final fit at optimal hyperparameters
  # ----------------------------------------------------------------
  hp_opt <- unpack_hyperparams(opt$par, theta_names, spatial_type, has_shared,
                                family$family)
  spatial_args_opt <- build_spatial_args(hp_opt, spatial_info,
                                          model_data$n_spatial, has_shared)

  final_result <- cpp_joint_laplace_fit(
    y = model_data$y,
    n_trials = model_data$n_trials,
    X_eco = model_data$X_eco,
    family = family$family,
    phi_overdispersion = hp_opt$phi,
    Z_samp = model_data$Z_samp,
    Z_quad = model_data$Z_quad,
    quad_weights = model_data$quad_weights,
    quad_spatial_idx = model_data$quad_spatial_idx,
    spatial_type = spatial_type,
    obs_spatial_idx = model_data$obs_spatial_idx,
    gp_L_shared_flat = spatial_args_opt$gp_L_shared,
    gp_L_eco_flat = spatial_args_opt$gp_L_eco,
    gp_L_samp_flat = spatial_args_opt$gp_L_samp,
    adj_row_ptr = spatial_info$adj_row_ptr,
    adj_col_idx = spatial_info$adj_col_idx,
    n_neighbors_vec = spatial_info$n_neighbors,
    tau_shared = hp_opt$tau_shared,
    tau_eco = hp_opt$tau_eco,
    tau_samp = hp_opt$tau_samp,
    spde_Q_row_ptr = spatial_args_opt$spde_Q_row_ptr,
    spde_Q_col_idx = spatial_args_opt$spde_Q_col_idx,
    spde_Q_shared_vals = spatial_args_opt$spde_Q_shared_vals,
    spde_Q_eco_vals = spatial_args_opt$spde_Q_eco_vals,
    spde_Q_samp_vals = spatial_args_opt$spde_Q_samp_vals,
    latent_type = latent_type,
    beta_share = hp_opt$beta_share,
    p_eco = model_data$p_eco,
    p_samp = model_data$p_samp,
    n_spatial = model_data$n_spatial,
    tau_beta = tau_beta,
    max_iter = max_newton,
    tol = newton_tol
  )

  # ----------------------------------------------------------------
  # 5. Sample from Laplace approximation
  # ----------------------------------------------------------------
  if (verbose) message("Sampling from Laplace approximation...")

  seed_val <- seed %||% sample.int(.Machine$integer.max, 1)

  samples <- cpp_joint_laplace_sample(
    mode = final_result$mode,
    L_flat = final_result$L_flat,
    n_params = final_result$n_params,
    n_samples = n_samples,
    seed = seed_val
  )

  # Name the columns
  col_names <- make_param_names(model_data, has_shared)
  colnames(samples) <- col_names

  # ----------------------------------------------------------------
  # 6. Return results
  # ----------------------------------------------------------------
  list(
    draws = samples,
    mode = final_result$mode,
    hyperparams = hp_opt,
    hyperparams_raw = opt$par,
    log_marginal = -opt$value,
    log_lik_eco = final_result$log_lik_eco,
    log_lik_samp = final_result$log_lik_samp,
    converged = final_result$converged && (opt$convergence == 0),
    n_newton_iter = final_result$n_iter,
    n_optim_evals = eval_count,
    n_params = final_result$n_params
  )
}

#' Unpack hyperparameter vector into named list
#' @keywords internal
unpack_hyperparams <- function(theta, names, spatial_type, has_shared, family) {
  hp <- as.list(theta)
  names(hp) <- names

  out <- list(
    tau_shared = 0, tau_eco = 0, tau_samp = 0,
    range_eco = 0, sigma_eco = 0,
    range_samp = 0, sigma_samp = 0,
    range_shared = 0, sigma_shared = 0,
    beta_share = 0, phi = 1.0
  )

  if (spatial_type == 0L || spatial_type == 3L) {  # GP or SPDE
    out$range_eco <- exp(hp$log_range_eco)
    out$sigma_eco <- exp(hp$log_sigma_eco)
    out$range_samp <- exp(hp$log_range_samp)
    out$sigma_samp <- exp(hp$log_sigma_samp)
    if (has_shared) {
      out$range_shared <- exp(hp$log_range_shared)
      out$sigma_shared <- exp(hp$log_sigma_shared)
    }
  } else if (spatial_type == 1L || spatial_type == 2L) {  # ICAR/BYM2
    out$tau_eco <- exp(hp$log_tau_eco)
    out$tau_samp <- exp(hp$log_tau_samp)
    if (has_shared) {
      out$tau_shared <- exp(hp$log_tau_shared)
    }
  }

  if (has_shared && !is.null(hp$beta_share)) {
    out$beta_share <- hp$beta_share
  }

  if (family == "negbin" && !is.null(hp$log_phi)) {
    out$phi <- exp(hp$log_phi)
  }

  out
}

#' Build spatial-specific arguments for C++ call
#' @keywords internal
build_spatial_args <- function(hp, spatial_info, n_spatial, has_shared) {
  out <- list(
    gp_L_shared = numeric(0),
    gp_L_eco = numeric(0),
    gp_L_samp = numeric(0),
    spde_Q_row_ptr = integer(0),
    spde_Q_col_idx = integer(0),
    spde_Q_shared_vals = numeric(0),
    spde_Q_eco_vals = numeric(0),
    spde_Q_samp_vals = numeric(0)
  )

  if (spatial_info$spatial_type == 0L) {  # GP
    dist_mat <- spatial_info$dist_mat
    nu <- spatial_info$nu

    out$gp_L_eco <- build_gp_cholesky(dist_mat, hp$range_eco,
                                        hp$sigma_eco, nu)
    out$gp_L_samp <- build_gp_cholesky(dist_mat, hp$range_samp,
                                         hp$sigma_samp, nu)
    if (has_shared) {
      out$gp_L_shared <- build_gp_cholesky(dist_mat, hp$range_shared,
                                             hp$sigma_shared, nu)
    }

  } else if (spatial_info$spatial_type == 3L) {  # SPDE
    spde_eco <- build_spde_precision(spatial_info, hp$range_eco, hp$sigma_eco)
    spde_samp <- build_spde_precision(spatial_info, hp$range_samp, hp$sigma_samp)

    out$spde_Q_row_ptr <- spde_eco$row_ptr
    out$spde_Q_col_idx <- spde_eco$col_idx
    out$spde_Q_eco_vals <- spde_eco$vals
    out$spde_Q_samp_vals <- spde_samp$vals

    if (has_shared) {
      spde_shared <- build_spde_precision(spatial_info, hp$range_shared,
                                           hp$sigma_shared)
      out$spde_Q_shared_vals <- spde_shared$vals
    }
  }

  out
}

#' Generate parameter column names
#' @keywords internal
make_param_names <- function(model_data, has_shared) {
  eco_names <- paste0("beta_eco_", colnames(model_data$X_eco))
  samp_names <- paste0("beta_samp_", colnames(model_data$Z_samp))

  n_sp <- model_data$n_spatial
  if (has_shared) {
    spatial_names <- c(
      paste0("S_shared_", seq_len(n_sp)),
      paste0("S_eco_", seq_len(n_sp)),
      paste0("S_samp_", seq_len(n_sp))
    )
  } else {
    spatial_names <- c(
      paste0("S_eco_", seq_len(n_sp)),
      paste0("S_samp_", seq_len(n_sp))
    )
  }

  c(eco_names, samp_names, spatial_names)
}
