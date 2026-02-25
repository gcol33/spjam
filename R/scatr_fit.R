#' scatr_fit S3 class
#'
#' Output of [scatr()] containing posterior draws, model specifications,
#' and diagnostics from the Laplace approximation.
#'
#' @name scatr_fit
#' @keywords internal
NULL

#' Construct a scatr_fit object
#' @keywords internal
new_scatr_fit <- function(laplace_result, formula, family, spatial, shared,
                           priors, data, locations, model_data, spatial_info,
                           backend, seed) {
  structure(
    list(
      draws = laplace_result$draws,
      mode = laplace_result$mode,
      hyperparams = laplace_result$hyperparams,
      log_marginal = laplace_result$log_marginal,
      formula = formula,
      family = family,
      spatial = spatial,
      shared = shared,
      priors = priors,
      data = data,
      locations = locations,
      backend = backend,
      converged = laplace_result$converged,
      n_params = laplace_result$n_params,
      n_obs = nrow(data),
      n_spatial = model_data$n_spatial,
      p_eco = model_data$p_eco,
      p_samp = model_data$p_samp,
      .internal = list(
        model_data = model_data,
        spatial_info = spatial_info,
        seed = seed
      )
    ),
    class = "scatr_fit"
  )
}

#' @export
print.scatr_fit <- function(x, ...) {
  cat("scatR model fit\n")
  cat("  Family:     ", x$family$family, "(", x$family$link, ")\n")
  cat("  Spatial:    ", x$spatial$type, "\n")
  cat("  Latent:     ", x$shared$type, "\n")
  cat("  Backend:    ", x$backend, "\n")
  cat("  N obs:      ", x$n_obs, "\n")
  cat("  N spatial:  ", x$n_spatial, "\n")
  cat("  N params:   ", x$n_params, "\n")
  cat("  Converged:  ", x$converged, "\n")
  cat("  Log-marg:   ", sprintf("%.2f", x$log_marginal), "\n")

  if (x$shared$type == "shared") {
    bs <- x$hyperparams$beta_share
    cat("  beta_share: ", sprintf("%.3f", bs), "\n")
  }

  cat("\nFixed effects (posterior mode):\n")
  eco_idx <- seq_len(x$p_eco)
  samp_idx <- x$p_eco + seq_len(x$p_samp)
  eco_names <- colnames(x$.internal$model_data$X_eco)
  samp_names <- colnames(x$.internal$model_data$Z_samp)

  cat("  Ecological:\n")
  for (i in seq_along(eco_idx)) {
    cat(sprintf("    %-20s  %.4f\n", eco_names[i], x$mode[eco_idx[i]]))
  }

  cat("  Sampling:\n")
  for (i in seq_along(samp_idx)) {
    cat(sprintf("    %-20s  %.4f\n", samp_names[i], x$mode[samp_idx[i]]))
  }

  if (x$shared$type == "shared") {
    cat("\nHyperparameters:\n")
    hp <- x$hyperparams
    if (hp$range_shared > 0) {
      cat(sprintf("  range (shared): %.3f\n", hp$range_shared))
      cat(sprintf("  sigma (shared): %.3f\n", hp$sigma_shared))
    }
    if (hp$range_eco > 0) {
      cat(sprintf("  range (eco):    %.3f\n", hp$range_eco))
      cat(sprintf("  sigma (eco):    %.3f\n", hp$sigma_eco))
    }
    if (hp$range_samp > 0) {
      cat(sprintf("  range (samp):   %.3f\n", hp$range_samp))
      cat(sprintf("  sigma (samp):   %.3f\n", hp$sigma_samp))
    }
  }

  invisible(x)
}

#' @export
summary.scatr_fit <- function(object, ...) {
  # Posterior summaries from draws
  draws <- object$draws
  eco_idx <- seq_len(object$p_eco)
  samp_idx <- object$p_eco + seq_len(object$p_samp)

  summarize_params <- function(idx, names) {
    data.frame(
      parameter = names,
      mean = colMeans(draws[, idx, drop = FALSE]),
      sd = apply(draws[, idx, drop = FALSE], 2, stats::sd),
      q025 = apply(draws[, idx, drop = FALSE], 2, stats::quantile, 0.025),
      q975 = apply(draws[, idx, drop = FALSE], 2, stats::quantile, 0.975),
      row.names = NULL
    )
  }

  eco_names <- colnames(object$.internal$model_data$X_eco)
  samp_names <- colnames(object$.internal$model_data$Z_samp)

  eco_summary <- summarize_params(eco_idx, eco_names)
  samp_summary <- summarize_params(samp_idx, samp_names)

  cat("scatR model summary\n\n")
  cat("Ecological process:\n")
  print(eco_summary, digits = 4, row.names = FALSE)
  cat("\nSampling process:\n")
  print(samp_summary, digits = 4, row.names = FALSE)

  if (object$shared$type == "shared") {
    cat(sprintf("\nCoupling (beta_share): %.3f\n", object$hyperparams$beta_share))
  }

  cat(sprintf("\nLog-marginal likelihood: %.2f\n", object$log_marginal))

  invisible(list(ecological = eco_summary, sampling = samp_summary))
}

#' @export
coef.scatr_fit <- function(object, type = c("ecological", "sampling"), ...) {
  type <- match.arg(type)
  if (type == "ecological") {
    idx <- seq_len(object$p_eco)
    out <- colMeans(object$draws[, idx, drop = FALSE])
    names(out) <- colnames(object$.internal$model_data$X_eco)
  } else {
    idx <- object$p_eco + seq_len(object$p_samp)
    out <- colMeans(object$draws[, idx, drop = FALSE])
    names(out) <- colnames(object$.internal$model_data$Z_samp)
  }
  out
}

#' @export
logLik.scatr_fit <- function(object, ...) {
  val <- object$log_marginal
  attr(val, "df") <- length(object$hyperparams)
  attr(val, "nobs") <- object$n_obs
  class(val) <- "logLik"
  val
}
