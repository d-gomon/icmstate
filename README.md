
<!-- README.md is generated from README.Rmd. Please edit that file -->

# icmstate

<!-- <img src="man/figures/success_hex.png" align="right" width="120" /> -->
<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://lifecycle.r-lib.org/articles/stages.html)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/concept.svg)](https://www.repostatus.org/#concept)
[![GitHub Repo
stars](https://img.shields.io/github/stars/d-gomon/icmstate?style=social)](https://github.com/d-gomon/icmstate)
[![Codecov test
coverage](https://codecov.io/gh/d-gomon/icmstate/graph/badge.svg)](https://app.codecov.io/gh/d-gomon/icmstate)
[![arXiv](https://img.shields.io/badge/stat.AP-arXiv%3A2409.07176-B31B1B)](https://doi.org/10.48550/arXiv.2409.07176)
<!-- [![R-CMD-check](https://github.com/d-gomon/success/workflows/R-CMD-check/badge.svg)](https://github.com/d-gomon/success/actions/) -->
<!-- [![CRAN status](https://www.r-pkg.org/badges/version/success)](https://CRAN.R-project.org/package=success) -->
<!-- [![CRAN RStudio mirror downloads](https://cranlogs.r-pkg.org/badges/success)](https://cran.r-project.org/package=success) -->

<!-- [![Biostatistics](https://img.shields.io/badge/Biostatistics-kxac041-%23003365)](https://doi.org/10.1093/biostatistics/kxac041) -->
<!-- [![Codecov test coverage](https://codecov.io/gh/d-gomon/success/branch/main/graph/badge.svg)](https://codecov.io/gh/d-gomon/success?branch=main) -->
<!-- badges: end -->

# Interval-Censored Multi-STATE modelling (icmstate)

Using the package it is possible to non-parametrically estimate
transition intensities in interval-censored multi-state models. The
package also provides utilities for plotting panel data and transition
probabilities between states. Many of the functions from the package can
be used for visualisation as well. The package currently supports two
approaches to determine transition intensities, either using the
Multinomial likelihood approach as in [Gomon and Putter
(2024)](https://doi.org/10.48550/arXiv.2409.07176) or the latent Poisson
approach as in [Gu et al
(2023)](https://doi.org/10.1093/biomet/asad073).

## Installation

Install the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("d-gomon/icmstate")
```

## Research Code

The results of [Gomon and Putter
(2024)](https://doi.org/10.48550/arXiv.2409.07176) can be reproduced
using the code in the folder `Research_Code`.

## References

Gomon D., Putter H. (2024): [Non-parametric estimation of transition
intensities in interval censored Markov multi-state models without
loops](https://doi.org/10.48550/arXiv.2409.07176), *arXiv*

Gu Y., Zeng D., Heiss G., Lin D. Y. (2023): [Maximum likelihood
estimation for semiparametric regression models with interval-censored
multistate data](https://doi.org/10.1093/biomet/asad073), *Biometrika*
