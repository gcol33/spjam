# scatR

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/gcol33/scatR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gcol33/scatR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**Joint spatial modelling of ecological processes under preferential sampling**

Spatial ecological data rarely represent where ecological processes truly occur. Instead, they reflect where observers choose or are able to sample. When sampling locations are non-random and correlated with the ecological process, standard spatial models confound sampling structure with ecological signal, producing biased maps that appear plausible but are structurally wrong.

scatR addresses this by jointly modelling sampling and ecology as two coupled spatial processes.

## Installation

```r
# install.packages("pak")
pak::pak("gcol33/scatR")
```

## Core idea

Observed spatial patterns are *scattered* versions of the true ecological process. Sampling acts as a distortion operator, redistributing observations toward accessible regions. scatR explicitly models this scattering, allowing the underlying ecological surface to be inferred separately from sampling behaviour.

## Modelling framework

scatR implements a Bayesian hierarchical pipeline with:

- A **latent ecological spatial process**
- A **latent sampling intensity process**
- **Shared spatial structure** to identify preferential sampling
- A **joint likelihood** tying observations to both processes

The result is a sampling-adjusted ecological surface that reflects where the process is supported to exist, not merely where it was observed.

## Usage

```r
library(scatR)

# Fit joint model
fit <- scatr(

abundance ~ elevation + forest_cover,
  sampling = ~ road_distance + population_density,
  data = species_obs,
  locations = ~ longitude + latitude,
  spatial = spatial_gp(),
  shared = latent_shared()
)

# Check for preferential sampling
check_identifiability(fit)

# Decompose spatial structure
decomp <- scatter_decomposition(fit)
plot(decomp)

# Sampling-adjusted predictions
predict(fit, newdata = grid, type = "ecological")
```

## How scatR differs from existing approaches

| Approach | Assumption | Problem | scatR solution |
|----------|------------|---------|----------------|
| **Spatial GLMMs/GAMs** | Sampling locations exogenous | Random effects absorb bias, not signal | Models sampling explicitly |
| **Bias covariates** | Proxies capture effort | Incomplete; residual bias remains | Treats sampling as spatial process |
| **Pseudo-absences** | Background approximates availability | Sensitive to design; doesn't model mechanism | Uses observed locations; models sampling directly |
| **Post-hoc sensitivity** | Compare multiple models | Retrospective; doesn't change likelihood | Integrates sampling into single model |

## Positioning

scatR is not a general-purpose spatial modelling framework. It is a deliberately constrained pipeline for situations where sampling and ecology are entangled—analogous to how occupancy models formalize detection bias.

**spOccupancy** models imperfect detection conditional on sampled sites (separating presence from detection). **scatR** addresses a different but equally fundamental bias: non-random sampling locations. Where spOccupancy corrects for failure to detect organisms at sampled sites, scatR corrects for the fact that the sites themselves are not a random sample of space.

## Status

This package is **experimental**. The API may change. See [ROADMAP.md](ROADMAP.md) for the development plan from v0.0.1 to v1.0.0.

## License

MIT
