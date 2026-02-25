# spjam Implementation Roadmap

**From proof-of-concept to production-ready release**

---

## Design Philosophy

spjam is an opinionated pipeline, not a framework. Every design decision prioritizes:

1. **Identifiability over flexibility** — structure is enforced, not optional
2. **Honest inference** — flag what cannot be distinguished rather than pretend
3. **Applied use** — familiar interface, sensible defaults, clear outputs
4. **Computational pragmatism** — custom HMC backend, no Stan compilation overhead

---

## Version Summary

| Version | Focus | Key Deliverables |
|---------|-------|------------------|
| 0.0.1 | Scaffold | Package structure, API stubs, C++ skeleton |
| 0.1.0 | Core engine | HMC sampler, autodiff, basic spatial GP |
| 0.2.0 | Joint model | Coupled ecological + sampling processes |
| 0.3.0 | Shared latent | Identifiable preferential sampling structure |
| 0.4.0 | Model variants | Point/areal data, alternative spatial specs |
| 0.5.0 | Diagnostics | Identifiability checks, decomposition tools |
| 0.6.0 | Prediction | Sampling-adjusted spatial prediction |
| 0.7.0 | Inference options | Laplace approximation, mesh refinement |
| 0.8.0 | Polish | Validation, edge cases, user feedback |
| 0.9.0 | Documentation | Vignettes, examples, pkgdown site |
| 1.0.0 | Release | CRAN submission |

---

## Phase 1: Foundation (0.0.1 → 0.1.0)

### 0.0.1 — Package Scaffold

**Objective:** Establish package structure mirroring ratiod architecture.

- [x] DESCRIPTION, NAMESPACE, LICENSE
- [x] Directory structure: `R/`, `src/`, `tests/`, `inst/`, `man/`, `vignettes/`
- [x] Core function stubs with roxygen skeletons
- [ ] C++ skeleton: `src/spjam_init.cpp`, Makevars
- [ ] Basic test infrastructure with testthat
- [ ] GitHub Actions for R CMD check

**Files created:**
```
R/spjam-package.R      # Package documentation
R/spjam.R              # Main fitting function
R/formula.R            # Formula parsing (stub)
R/family.R             # Family definitions (stub)
R/priors.R             # Prior system (stub)
R/spatial.R            # Spatial specifications (stub)
R/sampling.R           # Sampling process specs (stub)
R/latent.R             # Shared structure specs (stub)
src/spjam_init.cpp     # C++ registration
src/Makevars           # Compilation flags
```

### 0.1.0 — Core Computational Engine

**Objective:** Working HMC sampler with autodiff for simple models.

**Sampler implementation:**
- Custom HMC/NUTS in C++ with Rcpp
- Dual averaging for step size adaptation
- Numerical gradient computation (stable baseline)
- Mass matrix estimation (diagonal, dense options)

**Autodiff system:**
- Scalar autodiff class with operator overloading
- Forward-mode AD for gradient computation
- Support for log, exp, sqrt, pow, trigonometric functions
- Sparse Jacobian computation for latent fields

**Spatial primitives:**
- Gaussian process with Matérn covariance
- Efficient computation via Cholesky decomposition
- Distance matrix utilities for irregular locations

**Basic inference test:**
- Fit single-process spatial GP to simulated data
- Verify posterior recovery
- Benchmark against reference implementation

**Deliverable:** `spjam(y ~ 1, spatial = spatial_gp())` produces valid posterior samples for a single spatial process.

---

## Phase 2: Joint Model Structure (0.2.0 → 0.3.0)

### 0.2.0 — Coupled Processes

**Objective:** Joint likelihood with ecological + sampling components.

**Two-process architecture:**
```
Ecological:  Y_i | S_eco(x_i), X_i ~ f(mu_eco)
Sampling:    N(A) | S_samp(A) ~ Poisson(∫_A λ(s) ds)
```

**Implementation:**
- Separate linear predictors for each process
- Separate spatial random effects (initially independent)
- Combined log-likelihood in C++
- Parameter blocking for efficient sampling

**Sampling process options:**
- `sampling_lgcp()`: Log-Gaussian Cox process for point patterns
- `sampling_binomial()`: Grid-based presence/absence of sampling effort
- `sampling_point()`: Inhomogeneous Poisson process likelihood

**Numerical integration:**
- Quadrature for LGCP intensity over domain
- Mesh-based integration with user-specified resolution
- Boundary handling for irregular domains

**Validation:**
- Simulate from known joint process
- Recover both ecological and sampling parameters
- Compare to independent-fit baseline

### 0.3.0 — Shared Latent Structure

**Objective:** Identifiable shared spatial component for preferential sampling.

**Core model structure:**
```
S_eco(x)  = S_shared(x) + S_eco_only(x)
S_samp(x) = β_share * S_shared(x) + S_samp_only(x)
```

Where `β_share` measures the coupling between sampling and ecology.

**Latent specifications:**
- `latent_shared()`: Default shared structure
- `latent_independent()`: No sharing (comparison model)
- `latent_correlated(rho)`: Correlation-based coupling

**Identifiability constraints:**
- Sum-to-zero on shared component
- Variance partitioning priors
- Prior on coupling strength with shrinkage toward zero

**Diagnostics (preliminary):**
- Posterior correlation between ecological and sampling surfaces
- Effective number of shared parameters
- Warning when shared structure dominates

**Validation:**
- Simulation with known coupling strength
- Recovery of β_share under various scenarios
- Comparison of ecological surface with/without correction

---

## Phase 3: Model Extensions (0.4.0 → 0.5.0)

### 0.4.0 — Data Types and Spatial Specifications

**Objective:** Support for common ecological data structures.

**Response types:**
- Point-referenced observations (continuous locations)
- Areal data (counts in regions)
- Presence-only records
- Abundance with zeros (zero-inflated, hurdle)

**Family extensions:**
- `spjam_family("poisson")`: Count data
- `spjam_family("binomial")`: Presence/absence
- `spjam_family("negbin")`: Overdispersed counts
- `spjam_family("tweedie")`: Biomass, density
- `spjam_family("zi_poisson")`: Zero-inflated Poisson

**Alternative spatial specifications:**
- `spatial_gp()`: Full Gaussian process (small n)
- `spatial_spde()`: SPDE approximation via mesh (large n)
- `spatial_car()`: CAR for areal data
- `spatial_bym2()`: Leroux-type structured + unstructured

**Mesh construction:**
- Integration with sf for geometry handling
- Automatic mesh generation for SPDE
- User-specified mesh refinement

### 0.5.0 — Identifiability Diagnostics

**Objective:** Explicit tools for assessing what the data can support.

**Identifiability checks:**
- `check_identifiability(fit)`: Structured assessment of separation
- Prior sensitivity analysis for shared structure
- Effective degrees of freedom decomposition

**Scatter decomposition:**
```r
decomp <- scatter_decomposition(fit)
# Returns:
#   - ecological_only: spatial structure attributable to ecology
#   - sampling_only: spatial structure attributable to sampling
#   - shared: structure common to both (the confounding)
#   - residual: unexplained variation
```

**Variance partitioning:**
- Proportion of spatial variance in each component
- Uncertainty in partitioning (posterior distributions)
- Visualization methods

**Warnings and flags:**
- Weak identifiability warning when posterior on β_share is diffuse
- Collinearity diagnostics for fixed effects
- Prior-data conflict detection

---

## Phase 4: Inference and Prediction (0.6.0 → 0.7.0)

### 0.6.0 — Sampling-Adjusted Prediction

**Objective:** Spatial prediction that reflects ecological process, not sampling.

**Prediction types:**
```r
predict(fit, newdata, type = "ecological")   # Sampling-corrected surface
predict(fit, newdata, type = "sampling")     # Sampling intensity
predict(fit, newdata, type = "observed")     # What we'd expect to observe
predict(fit, newdata, type = "shared")       # Shared component only
```

**Grid prediction:**
- Efficient batch prediction over regular grids
- Integration with terra/stars for raster output
- Credible interval surfaces

**Counterfactual surfaces:**
- "What if sampling were uniform?"
- "What if sampling followed accessibility only?"
- Comparison framework for hypothesis testing

**Fitted values:**
```r
fitted(fit, type = "ecological")  # At observed locations
fitted(fit, type = "sampling")    # Sampling intensity at locations
```

### 0.7.0 — Computational Efficiency

**Objective:** Scale to large datasets.

**Laplace approximation backend:**
- Modal approximation for ecological process
- Gaussian approximation to posterior
- Orders of magnitude faster for large n

**Backend selection:**
```r
spjam(..., backend = "auto")     # Selects based on n and model complexity
spjam(..., backend = "hmc")      # Full MCMC (default)
spjam(..., backend = "laplace")  # Fast approximate inference
```

**SPDE optimizations:**
- Sparse precision matrices
- Mesh refinement strategies
- Approximate likelihoods for very large n

**Parallelization:**
- Multi-chain parallelism
- Within-chain parallelism for large linear algebra
- GPU support (future consideration)

---

## Phase 5: Production Quality (0.8.0 → 0.9.0)

### 0.8.0 — Robustness and Edge Cases

**Objective:** Handle real-world messiness.

**Input validation:**
- Comprehensive type and range checking
- Informative error messages
- Graceful handling of missing data

**Edge cases:**
- Very sparse data
- Highly unbalanced sampling
- Near-collinear covariates
- Boundary effects in spatial prediction

**Numerical stability:**
- Scaled parameterization for HMC
- Condition number monitoring
- Fallback algorithms for ill-conditioned problems

**Convergence diagnostics:**
- Rhat, ESS, divergences (via posterior package)
- Automatic warnings for poor convergence
- Remediation suggestions

### 0.9.0 — Documentation and Examples

**Objective:** Comprehensive user guidance.

**Vignettes:**
1. `introduction.Rmd`: What is preferential sampling, why spjam?
2. `getting-started.Rmd`: First model in 10 minutes
3. `model-specification.Rmd`: Formula, family, spatial, sampling options
4. `interpretation.Rmd`: Understanding output, decomposition, prediction
5. `case-study-species.Rmd`: Full worked example with species distribution data
6. `case-study-disease.Rmd`: Epidemiological application
7. `comparison.Rmd`: spjam vs. standard spatial models, when to use each
8. `troubleshooting.Rmd`: Common issues and solutions

**Function documentation:**
- Complete roxygen documentation for all exported functions
- Examples that run in < 5 seconds
- Cross-references between related functions

**pkgdown site:**
- Themed documentation site
- Reference index with logical grouping
- Articles from vignettes

**Simulation study:**
- Reproducible code for validating the method
- Comparison with alternative approaches
- Guidance on when spjam helps/doesn't help

---

## Phase 6: Release (1.0.0)

### 1.0.0 — CRAN Submission

**Objective:** Stable, documented, tested release.

**CRAN compliance:**
- All R CMD check --as-cran passes
- Examples run in reasonable time
- No external dependencies that fail checks
- Vignettes pre-built if computationally expensive

**Testing:**
- > 80% code coverage
- Unit tests for all core functions
- Integration tests for full workflows
- Regression tests for known edge cases

**API freeze:**
- No breaking changes after 1.0.0
- Deprecation policy for future changes
- Semantic versioning commitment

**Changelog:**
- Complete NEWS.md documenting all changes
- Migration guide from any pre-release versions

**Publication:**
- JSS or Methods in Ecology and Evolution submission
- Companion paper with simulation study
- Archived simulation code

---

## API Preview

### Core Workflow

```r
library(spjam)

# 1. Fit joint model
fit <- spjam(
  # Ecological process
  abundance ~ elevation + forest_cover + (1 | region),
  # Sampling process
  sampling = ~ road_distance + population_density,
  # Data and locations
  data = species_obs,
  locations = ~ longitude + latitude,
  # Spatial structure
  spatial = spatial_spde(mesh = my_mesh),
  # Shared latent structure (what we're trying to identify)
  shared = latent_shared(prior = prior_pc(3, 0.5)),
  # Inference
  chains = 4,
  iter = 2000
)

# 2. Check identifiability
check_identifiability(fit)
# → "Preferential sampling detected (β_share posterior: 0.42 [0.18, 0.71])"
# → "Ecological and sampling processes share 34% of spatial variance"

# 3. Decompose scatter
decomp <- scatter_decomposition(fit)
plot(decomp)

# 4. Get sampling-adjusted predictions
pred_eco <- predict(fit, newdata = prediction_grid, type = "ecological")
pred_obs <- predict(fit, newdata = prediction_grid, type = "observed")

# Compare: what we'd observe vs. what actually exists
plot(pred_obs$mean, pred_eco$mean)
abline(0, 1)  # Deviation = sampling bias

# 5. Summarize
summary(fit)
# → Fixed effects (ecological): ...
# → Fixed effects (sampling): ...
# → Spatial parameters: ...
# → Shared structure: ...
# → MCMC diagnostics: all Rhat < 1.01, min ESS = 412
```

### Key S3 Methods

```r
print(fit)                           # Concise summary
summary(fit)                         # Full summary with intervals
plot(fit)                            # Diagnostic plots
plot(fit, type = "spatial")          # Spatial surface plot

fitted(fit)                          # Fitted values at observed locations
fitted(fit, type = "ecological")     # Ecological process only
fitted(fit, type = "sampling")       # Sampling intensity

predict(fit, newdata)                # Prediction at new locations
predict(fit, type = "ecological")    # Sampling-adjusted
predict(fit, type = "shared")        # Shared component

scatter_decomposition(fit)           # Variance partitioning
check_identifiability(fit)           # Identifiability assessment
mcmc_diagnostics(fit)                # Convergence diagnostics

loo(fit)                             # Leave-one-out CV
waic(fit)                            # WAIC
pp_check(fit)                        # Posterior predictive checks
```

---

## Dependencies

**Core:**
- Rcpp, RcppEigen (C++ interface)
- Matrix (sparse linear algebra)
- sf (spatial geometry)

**Suggested:**
- posterior, bayesplot, loo (MCMC diagnostics)
- ggplot2 (enhanced plotting)
- terra, stars (raster prediction)
- spdep (neighborhood structures)

**No dependency on:**
- Stan/rstan (custom HMC avoids compilation overhead)
- INLA (optional comparison, not required)
- TMB (considered, rejected for API consistency with ratiod)

---

## Milestones and Checkpoints

| Milestone | Checkpoint | Success Criteria |
|-----------|------------|------------------|
| 0.1.0 | Single-process GP | Posterior recovery in simulation |
| 0.2.0 | Joint likelihood | Both processes estimated, parameters recovered |
| 0.3.0 | Shared identification | β_share recovered, decomposition works |
| 0.5.0 | Applied ready | Works on real species distribution data |
| 0.7.0 | Scalable | 10,000+ observations in < 10 minutes |
| 0.9.0 | Documented | Complete vignettes, all functions documented |
| 1.0.0 | Released | CRAN accepted, paper submitted |

---

## Open Questions for Development

1. **SPDE vs. full GP default**: SPDE scales better but requires mesh specification. Should `spatial_gp()` auto-switch based on n?

2. **Sampling effort as covariate vs. process**: Should we support known effort layers as covariates, or always model sampling as a process?

3. **Discrete vs. continuous space**: How to handle presence-only data where "sampling locations" aren't explicit?

4. **Comparison framework**: Should spjam include built-in comparison with naive models (no preferential sampling correction)?

5. **Multiple species**: Extension to joint species distribution models with shared sampling? (Post-1.0)

---

## Success Metrics

**Technical:**
- Recovers known simulation parameters with nominal coverage
- Produces different predictions than naive spatial models when preferential sampling is present
- Scales to typical ecological dataset sizes (1,000–50,000 observations)

**Usability:**
- Fits a basic model in < 5 lines of code
- Clear diagnostic output that non-statisticians can interpret
- Error messages that suggest fixes

**Adoption:**
- Cited in applied ecological papers
- Used by researchers who previously ignored preferential sampling
- Recognized as the practical solution to a known problem

---

*Last updated: 2025-01*
