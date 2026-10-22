
# spjam

*species maps shaped by where people looked*

<!-- badges: start -->

[![R-CMD-check](https://github.com/gcol33/spjam/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gcol33/spjam/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/gcol33/spjam/graph/badge.svg)](https://app.codecov.io/gh/gcol33/spjam)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

spjam fits spatial models of an ecological process, such as abundance or
occurrence, together with a model of where observers sampled, so that
the map of the process is estimated after the sampling pattern is
accounted for.

When observers go where the species is common, sampling is preferential,
and a standard spatial model’s random effect absorbs the sampling
pattern along with the ecology. spjam writes the sampling intensity as a
second spatial process and couples the two through a shared latent
field. The coupling parameter `beta_share` measures how strongly
sampling follows the process, so you can judge whether sampling was
preferential and compare the map with and without the shared field.
spjam supports:

- **Several response families.** Poisson, negative binomial, binomial,
  Tweedie and zero-inflated responses.
- **Point and areal data.** Matérn Gaussian processes, SPDE meshes for
  many locations, and conditional autoregressive and BYM2 fields for
  areal units.
- **Choice of coupling.** A shared field that estimates preferential
  sampling, independent fields for the comparison fit, or distinct
  fields with an estimated correlation.
- **Sampling models.** Sampling as a log-Gaussian Cox process, as fixed
  points or as a binomial process.

The default Laplace backend optimises the spatial hyperparameters with
L-BFGS-B, finds the mode of the latent field with a Newton-Raphson loop
in C++ (RcppEigen), and draws from the Gaussian approximation at the
mode.

## Installation

spjam is at an early version and its interface may change; the Laplace
backend is the supported path
([roadmap](https://github.com/gcol33/spjam/blob/main/ROADMAP.md)).

The package compiles C++ code. Install the build tools for your
operating system first:

| Operating system | Build tools                                              |
|------------------|----------------------------------------------------------|
| Windows          | [Rtools](https://cran.r-project.org/bin/windows/Rtools/) |
| macOS            | Run `xcode-select --install` in Terminal                 |
| Debian / Ubuntu  | Install `r-base-dev`                                     |

Check that R can find the compiler:

``` r
pkgbuild::has_build_tools(debug = TRUE)
```

Then install spjam:

``` r
# install.packages("pak")
pak::pak("gcol33/spjam")
```

The code is compiled as C++17. The Eigen headers come from RcppEigen,
installed as a dependency, and the code links the LAPACK and BLAS
libraries that R itself uses, so no further system library is needed.

## Usage

1.  Fit the model with `spjam()`: a formula for the ecological response,
    the data, and the columns that hold the coordinates.
2.  Add a formula for the sampling process with `sampling` when
    covariates of observer effort are known.
3.  Compare a shared fit with an independent one through `shared`.

The input is a data frame with one row per observation and coordinate
columns named in `locations`.

``` r
library(spjam)

set.seed(1)
dat <- data.frame(
  abundance = rpois(30, exp(1.5)),
  elevation = rnorm(30),
  x = runif(30, 0, 5),
  y = runif(30, 0, 5)
)

fit <- spjam(abundance ~ elevation, data = dat, locations = ~ x + y,
             backend = "laplace")
```

The data are 30 simulated counts with an elevation covariate. The
printed fit gives the family, the fields, the coupling `beta_share`, and
the fixed effects of the ecological and sampling processes:

``` r
fit
#> spjam model fit
#>   Family:      poisson ( log )
#>   Spatial:     gp 
#>   Latent:      shared 
#>   Backend:     laplace 
#>   N obs:       30 
#>   N spatial:   30 
#>   N params:    93 
#>   Converged:   TRUE 
#>   Log-marg:    -100.85 
#>   beta_share:  0.139 
#> 
#> Fixed effects (posterior mode):
#>   Ecological:
#>     (Intercept)           1.5320
#>     elevation             0.0158
#>   Sampling:
#>     (Intercept)           -0.1170
#> 
#> Hyperparameters:
#>   range (shared): 11.529
#>   sigma (shared): 0.000
#>   range (eco):    11.529
#>   sigma (eco):    0.000
#>   range (samp):   11.529
#>   sigma (samp):   0.000
```

## Modelling where observers sampled

The two latent surfaces share one field:

    S_eco(x)  = S_shared(x) + S_eco_only(x)
    S_samp(x) = beta_share * S_shared(x) + S_samp_only(x)

`beta_share` is the part of the sampling surface that is the field
driving the ecology. A fit with `latent_shared()` estimates it, a fit
with `latent_independent()` sets it to zero, and the difference between
the two ecological surfaces is the correction for preferential sampling:

``` r
dat$road_distance <- rnorm(30)

fit_shared <- spjam(abundance ~ elevation, sampling = ~ road_distance,
                    data = dat, locations = ~ x + y,
                    shared = latent_shared(), backend = "laplace")
fit_indep <- spjam(abundance ~ elevation, sampling = ~ road_distance,
                   data = dat, locations = ~ x + y,
                   shared = latent_independent(), backend = "laplace")

logLik(fit_shared)
logLik(fit_indep)
```

## Model components

- Spatial fields: `spatial_gp()` (Matérn, nu = 0.5, 1.5 or 2.5),
  `spatial_spde()`, `spatial_car()` and `spatial_bym2()`.
- Coupling: `latent_shared()`, `latent_independent()` and
  `latent_correlated()`.
- Sampling process: `sampling_lgcp()`, `sampling_point()` and
  `sampling_binomial()`.
- Priors: `prior_normal()`, `prior_half_normal()`,
  `prior_half_cauchy()`, and the penalised-complexity priors
  `prior_pc_range()` and `prior_pc_sigma()` for the spatial
  hyperparameters.

## Standard spatial model or spjam?

|  | Spatial GLMM or GAM | spjam |
|----|----|----|
| Sampling locations | Treated as given | Modelled as a spatial process |
| Sampling pattern | Absorbed by the random effect | Estimated as `beta_share` |
| Needs a sampling model | No | Yes |
| Returns | One spatial surface | Coupled ecological and sampling surfaces |
| Suited to | Prediction at observed sites | Maps read as the ecological process |

The split follows how occupancy models separate detection from
occurrence, applied to the choice of where to sample instead of what was
detected there.

## Documentation

The help pages in R document every function: `?spjam` for the model and
its arguments, and `help(package = "spjam")` for the full index.

## Getting help

Questions and bug reports go to the [GitHub
issues](https://github.com/gcol33/spjam/issues). A bug report is easiest
to act on with a small reproducible example, for instance made with
[reprex](https://reprex.tidyverse.org/).

## Support

> “Software is like sex: it’s better when it’s free.” — Linus Torvalds

I’m a PhD student who builds R packages in my free time because I
believe good tools should be free and open. I started these projects for
my own work and figured others might find them useful too.

If this package saved you some time, buying me a coffee is a nice way to
say thanks. It helps with my coffee addiction.

[![Buy Me A
Coffee](https://img.shields.io/badge/-Buy%20me%20a%20coffee-FFDD00?logo=buymeacoffee&logoColor=black)](https://buymeacoffee.com/gcol33)

## License

spjam is licensed under MIT.

## Citation

To cite spjam in publications, please use:

``` text
Colling G (2026). _spjam: Joint Spatial Modelling of Ecological
Processes Under Preferential Sampling_. R package version 0.0.1,
<https://github.com/gcol33/spjam>.
```

BibTeX:

``` bibtex
@Manual{colling2026spjam,
  title = {spjam: Joint Spatial Modelling of Ecological Processes Under Preferential Sampling},
  author = {Gilles Colling},
  year = {2026},
  note = {R package version 0.0.1},
  url = {https://github.com/gcol33/spjam},
}
```
