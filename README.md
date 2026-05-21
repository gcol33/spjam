# spjam <sub>Spatial Joint Attribution Models</sub>

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/gcol33/spjam/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gcol33/spjam/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Joint Spatial Modelling of Ecological Processes Under Preferential Sampling**

Spatial ecological data rarely represent where ecological processes truly occur. Instead, they reflect where observers choose or are able to sample. When sampling locations are non-random and correlated with the ecological process, standard spatial models confound sampling structure with ecological signal, producing biased maps that appear plausible but are structurally wrong.

`spjam` addresses this by jointly modelling sampling and ecology as two coupled spatial processes, producing sampling-adjusted predictions that reflect where the process is supported to exist rather than where it was observed.

## Quick Start

```r
library(spjam)

# Fit joint model: ecological process + sampling intensity
fit <- spjam(
  abundance ~ elevation + forest_cover,
  sampling = ~ road_distance + population_density,
  data = species_obs,
  locations = ~ longitude + latitude,
  spatial = spatial_gp(),
  shared = latent_shared()
)

# Inspect results
summary(fit)
coef(fit)
```

## Statement of Need

Observed spatial patterns are *scattered* versions of the true ecological process. Sampling acts as a distortion operator, redistributing observations toward accessible regions. Standard spatial models (GLMMs, GAMs) treat sampling locations as exogenous, allowing random effects to absorb bias rather than signal. Bias covariates and pseudo-absence approaches offer partial fixes, but residual confounding persists whenever sampling correlates with the ecological surface.

`spjam` models this scattering explicitly via a Bayesian hierarchical pipeline with:

- A **latent ecological spatial process** (the signal)
- A **latent sampling intensity process** (the distortion)
- **Shared spatial structure** to identify preferential sampling
- A **joint likelihood** tying observations to both processes

The result is a sampling-adjusted ecological surface. This is analogous to how occupancy models formalize detection bias, but at the level of sampling locations rather than detection events.

| Approach | Assumption | Problem | spjam solution |
|----------|------------|---------|----------------|
| **Spatial GLMMs/GAMs** | Sampling locations exogenous | Random effects absorb bias, not signal | Models sampling explicitly |
| **Bias covariates** | Proxies capture effort | Incomplete; residual bias remains | Treats sampling as spatial process |
| **Pseudo-absences** | Background approximates availability | Sensitive to design; doesn't model mechanism | Uses observed locations; models sampling directly |
| **Post-hoc sensitivity** | Compare multiple models | Retrospective; doesn't change likelihood | Integrates sampling into single model |

## Features

### Core Modelling

- **`spjam()`**: Fit a joint spatial model with ecological and sampling processes
  - Multiple response families: Poisson, negative binomial, binomial, Tweedie, zero-inflated
  - Backends: Laplace approximation (fast) or HMC (full posterior)

### Spatial Structure

- **`spatial_gp()`**: Gaussian process with Matern covariance (nu = 0.5, 1.5, 2.5)
- **`spatial_spde()`**: SPDE approximation for computational efficiency with large datasets
- **`spatial_car()`**: Conditional autoregressive for areal data (ICAR or proper)
- **`spatial_bym2()`**: BYM2 model combining structured and unstructured spatial components

### Latent Structure

- **`latent_shared()`**: Shared latent field between ecology and sampling (identifies preferential sampling)
- **`latent_independent()`**: Independent spatial processes (comparison model)
- **`latent_correlated()`**: Correlated but distinct spatial processes with estimated correlation

### Sampling Process

- **`sampling_lgcp()`**: Log-Gaussian Cox process for point-referenced data
- **`sampling_point()`**: Inhomogeneous Poisson point process
- **`sampling_binomial()`**: Binomial model for grid-based sampling effort

### Priors

- **`prior_normal()`** / **`prior_half_normal()`** / **`prior_half_cauchy()`**: Standard priors
- **`prior_pc_range()`** / **`prior_pc_sigma()`**: Penalized complexity priors for spatial hyperparameters

## Installation

```r
install.packages("pak")
pak::pak("gcol33/spjam")
```

## Usage Examples

### Basic Joint Model

```r
library(spjam)

fit <- spjam(
  abundance ~ elevation + forest_cover,
  sampling = ~ road_distance + population_density,
  data = species_obs,
  locations = ~ longitude + latitude,
  spatial = spatial_gp(nu = 1.5),
  shared = latent_shared()
)

summary(fit)
```

### Areal Data with BYM2

```r
fit <- spjam(
  cases ~ poverty + temperature,
  sampling = ~ hospital_density,
  data = district_data,
  locations = ~ longitude + latitude,
  spatial = spatial_bym2(adjacency = adj_matrix),
  shared = latent_shared()
)
```

### Comparing Shared vs Independent Models

```r
# Joint model (corrects for preferential sampling)
fit_shared <- spjam(
  abundance ~ elevation,
  sampling = ~ road_distance,
  data = obs, locations = ~ lon + lat,
  spatial = spatial_gp(),
  shared = latent_shared()
)

# Independent model (ignores preferential sampling)
fit_indep <- spjam(
  abundance ~ elevation,
  sampling = ~ road_distance,
  data = obs, locations = ~ lon + lat,
  spatial = spatial_gp(),
  shared = latent_independent()
)

# Compare via log-likelihood
logLik(fit_shared)
logLik(fit_indep)
```

## Status

This package is **experimental**. The API may change. See [ROADMAP.md](ROADMAP.md) for the development plan from v0.0.1 to v1.0.0.

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
