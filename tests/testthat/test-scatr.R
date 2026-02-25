test_that("scatr() fits a basic Poisson GP model", {
  set.seed(123)
  n <- 30
  dat <- data.frame(
    abundance = rpois(n, exp(1.5 + 0.3 * rnorm(n))),
    elevation = rnorm(n),
    x = runif(n, 0, 5),
    y = runif(n, 0, 5)
  )

  fit <- scatr(
    abundance ~ elevation,
    data = dat,
    locations = ~ x + y,
    family = scatr_family("poisson"),
    spatial = spatial_gp(nu = 1.5),
    shared = latent_shared(),
    backend = "laplace",
    control = list(
      max_newton = 50,
      max_outer_iter = 100,
      n_samples = 100,
      verbose = FALSE
    )
  )

  expect_s3_class(fit, "scatr_fit")
  expect_equal(fit$backend, "laplace")
  expect_equal(fit$n_obs, n)
  expect_equal(nrow(fit$draws), 100)
  expect_true(fit$n_params > 0)
  # Check draws have finite values
  expect_true(all(is.finite(fit$draws[, seq_len(fit$p_eco)])))
})

test_that("scatr() works with latent_independent()", {
  set.seed(456)
  n <- 20
  dat <- data.frame(
    abundance = rpois(n, 5),
    x = runif(n, 0, 3),
    y = runif(n, 0, 3)
  )

  fit <- scatr(
    abundance ~ 1,
    data = dat,
    locations = ~ x + y,
    shared = latent_independent(),
    backend = "laplace",
    control = list(
      max_newton = 15,
      max_outer_iter = 10,
      n_samples = 50,
      verbose = FALSE
    )
  )

  expect_s3_class(fit, "scatr_fit")
  # Independent model has fewer parameters (no shared field)
  expect_equal(fit$shared$type, "independent")
})

test_that("scatr() rejects HMC backend", {
  dat <- data.frame(abundance = rpois(10, 5), x = 1:10, y = 1:10)
  expect_error(
    scatr(abundance ~ 1, data = dat, locations = ~ x + y, backend = "hmc"),
    "HMC backend not yet implemented"
  )
})

test_that("scatr() validates inputs", {
  expect_error(
    scatr(abundance ~ 1, data = "not a data frame", locations = ~ x + y),
    "must be a data.frame"
  )

  dat <- data.frame(abundance = 1:5, x = 1:5, y = 1:5)
  expect_error(
    scatr(abundance ~ missing_var, data = dat, locations = ~ x + y),
    "not found in data"
  )
})

test_that("print and summary methods work", {
  set.seed(789)
  n <- 20
  dat <- data.frame(
    abundance = rpois(n, 5),
    x = runif(n, 0, 3),
    y = runif(n, 0, 3)
  )

  fit <- scatr(
    abundance ~ 1,
    data = dat,
    locations = ~ x + y,
    backend = "laplace",
    control = list(max_outer_iter = 5, n_samples = 50, verbose = FALSE)
  )

  expect_output(print(fit), "scatR model fit")
  expect_output(summary(fit), "Ecological process")
})

test_that("coef and logLik methods work", {
  set.seed(101)
  n <- 20
  dat <- data.frame(
    abundance = rpois(n, 5),
    elevation = rnorm(n),
    x = runif(n, 0, 3),
    y = runif(n, 0, 3)
  )

  fit <- scatr(
    abundance ~ elevation,
    data = dat,
    locations = ~ x + y,
    backend = "laplace",
    control = list(max_outer_iter = 5, n_samples = 50, verbose = FALSE)
  )

  eco_coefs <- coef(fit, type = "ecological")
  expect_named(eco_coefs, c("(Intercept)", "elevation"))
  expect_true(all(is.finite(eco_coefs)))

  ll <- logLik(fit)
  expect_true(is.finite(as.numeric(ll)))
})
