# scatR Development Guide

## Package Overview

Joint spatial modelling of ecological processes under preferential sampling.
Bayesian hierarchical model with Laplace approximation backend, C++ engine adapted from numdenom.

## numdenom Dependency Tracking

**Baseline version:** numdenom v1.3.0 (audited 2026-02-25)

C++ code adapted from numdenom into scatr:

| scatr file | numdenom source | What was adapted |
|---|---|---|
| `src/likelihood.h` | `src/laplace_core.cpp` lines 24-106 | Poisson, negbin, binomial log-lik + grad + neg_hess |
| `src/linalg_utils.h` | `src/linalg_fast.h` subset | safe_exp, dot_product, sparse triangular solves |
| `src/covariance.h` | `src/hmc_svc.h` covariance functions | exponential, Matérn 3/2, Matérn 5/2 |
| `src/spatial_precision.h` | `src/laplace_core.cpp` spatial sections + `src/hmc_gp.h` | ICAR, BYM2, GP NNGP, SPDE precision |
| `src/laplace_joint.h` | `src/laplace_core.h` + `src/hmc_sampler.h` (inspired) | JointModelData struct, param layout |
| `src/laplace_joint.cpp` | `src/laplace_core.cpp` Newton-Raphson sections | Extended to joint model (2 predictors, 3 fields) |

When numdenom updates beyond v1.3.0, audit the files above for:
- Bug fixes in likelihood functions
- Performance improvements in Laplace engine
- New spatial prior implementations
- Changes to NNGP or BYM2 scaling

## Build Commands

```r
Rcpp::compileAttributes()
devtools::document()
devtools::load_all()
devtools::test()
devtools::check()
```

## Architecture

```
User API: scatr() → validate → prepare_data → prepare_spatial → backend_laplace → C++ engine → scatr_fit
```

**Outer loop (R):** Hyperparameter optimization via `optim(method="L-BFGS-B")`
- Hyperparameters: range, sigma (per spatial field), β_share, overdispersion

**Inner loop (C++):** Newton-Raphson for latent field mode given hyperparameters
- Latent vector: [β_eco | β_samp | S_shared | S_eco_only | S_samp_only]

## Naming Conventions

- R functions: camelCase for exported (`scatr`, `scatrFormula`), snake_case for internal (`prepare_model_data`)
- C++ functions: snake_case (`joint_laplace_mode`)
- C++ structs: PascalCase (`JointModelData`)
- C++ constants: kPascalCase (`kMaxNewtonIter`)
