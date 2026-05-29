# spjam

*the map shows where people looked*

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/gcol33/spjam/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gcol33/spjam/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Joint spatial models that estimate the ecological process and the sampling process at once, coupled through a shared latent field.**

A spatial dataset records where observers went, not where the process lives. When the two are correlated, a standard spatial model lets its random effects soak up the sampling pattern and hands you a map of effort dressed as ecology. `spjam` fits both surfaces in one likelihood, ties them with a shared latent field, and estimates how strongly sampling tracks the process, so the ecological surface is what remains after the sampling is accounted for.

```r
library(spjam)

dat <- data.frame(
  abundance = rpois(30, exp(1.5)),
  elevation = rnorm(30),
  x = runif(30, 0, 5),
  y = runif(30, 0, 5)
)

fit <- spjam(
  abundance ~ elevation,
  data      = dat,
  locations = ~ x + y,
  backend   = "laplace"
)

print(fit)
```

## Two processes, one likelihood

A standard spatial GLMM or GAM treats sampling locations as fixed and given. Its spatial random effect then absorbs whatever pattern the sampling left behind, which is fine for prediction at observed sites and misleading as a picture of the process. `spjam` writes the sampling intensity as its own spatial process and fits it jointly with the ecology:

```
S_eco(x)  = S_shared(x) + S_eco_only(x)
S_samp(x) = beta_share * S_shared(x) + S_samp_only(x)
```

`beta_share` is the coupling: how much of the sampling surface is the same field that drives the ecology. It is the quantity that says whether sampling was preferential, and by how much. Compare a `latent_shared()` fit against a `latent_independent()` fit and the difference in the ecological surface is the correction.

## How it fits

The default `"laplace"` backend runs an outer L-BFGS-B optimisation over the spatial hyperparameters and a C++ Newton-Raphson inner loop (RcppEigen) for the latent field, then draws from the Gaussian approximation at the mode. No Stan, no compilation step at fit time.

```r
fit_shared <- spjam(abundance ~ elevation, sampling = ~ road_distance,
                    data = obs, locations = ~ lon + lat,
                    shared = latent_shared(), backend = "laplace")

fit_indep  <- spjam(abundance ~ elevation, sampling = ~ road_distance,
                    data = obs, locations = ~ lon + lat,
                    shared = latent_independent(), backend = "laplace")

logLik(fit_shared)
logLik(fit_indep)
```

## Building blocks

- **`spjam()`**: fit the joint model. Families: Poisson, negative binomial, binomial, Tweedie, zero-inflated.
- **Spatial fields**: `spatial_gp()` (Matern, nu = 0.5 / 1.5 / 2.5), `spatial_spde()` for large `n`, `spatial_car()` and `spatial_bym2()` for areal data.
- **Latent coupling**: `latent_shared()` (identifies preferential sampling), `latent_independent()` (the comparison fit), `latent_correlated()` (distinct fields with estimated correlation).
- **Sampling process**: `sampling_lgcp()`, `sampling_point()`, `sampling_binomial()`.
- **Priors**: `prior_normal()`, `prior_half_normal()`, `prior_half_cauchy()`, and penalized-complexity `prior_pc_range()` / `prior_pc_sigma()` for the spatial hyperparameters.

## Standard spatial model or spjam?

| | Spatial GLMM / GAM | `spjam` |
|---|---|---|
| Sampling locations | Treated as given | Modelled as a spatial process |
| Where bias goes | Into the random effect | Estimated as `beta_share` |
| Needs a sampling model? | No | Yes |
| Returns | One spatial surface | Coupled ecology + sampling surfaces |
| Best for | Prediction at observed sites | Maps meant to read as the process |

The split mirrors how occupancy models formalise detection bias, one level up: at the choice of where to sample rather than what was detected once there.

## Installation

```r
install.packages("pak")
pak::pak("gcol33/spjam")
```

## Status

`spjam` is experimental and at an early version; the API will change. The Laplace backend is the supported path. See [ROADMAP.md](ROADMAP.md) for the plan toward a 1.0.0 release.

## Support

> "Software is like sex: it's better when it's free." — Linus Torvalds

I'm a PhD student who builds R packages in my free time because I believe good tools should be free and open. I started these projects for my own work and figured others might find them useful too.

If this package saved you some time, buying me a coffee is a nice way to say thanks. It helps with my coffee addiction.

[![Buy Me A Coffee](https://img.shields.io/badge/-Buy%20me%20a%20coffee-FFDD00?logo=buymeacoffee&logoColor=black)](https://buymeacoffee.com/gcol33)

## License

MIT

## Citation

```bibtex
@software{spjam,
  author = {Colling, Gilles},
  title = {spjam: Joint Spatial Modelling Under Preferential Sampling},
  year = {2026},
  url = {https://github.com/gcol33/spjam}
}
```
