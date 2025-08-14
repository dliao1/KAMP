# Introduction
`kperm` is an R package for efficiently estimating Ripley’s K under complete spatial randomness using a robust statistical approach called **KAMP** (K adjusted by Analytical Moments of the Permutation distribution). 

KAMP avoids the need for computationally expensive permutations while accounting for spatial inhomogeneity, making it suitable for large-scale spatial analyses such as those encountered in spatial proteomics multiplex imaging datasets.

# Overview
This package provides functions to compute both **univariate** and **bivariate** KAMP expectation and variance (`spatstat` and matrix-based implementation both included). At this time, the only edge correction methods supported are translational (`trans`) and isotropic (`iso`).

# Basic Functionality and Examples
## Installation
```r
# Install from GitHub
devtools::install_github("dliao1/kperm")
```

## Univariate
```r
library(kperm)

# Simulate a point pattern
pp <- sim_pp_data(lambda_n = 500, abundance = 0.3)

# Compute KAMP expectation
kamp_expec_univ <- kamp_expectation(pp, markvar = "immune")
print(kamp_expec_univ)

# Compute KAMP variance
kamp_var <- kamp_variance(pp, markvar = "immune")
print(kamp_var)
```

## Bivariate
```r
library(kperm)

# Simulate a point pattern
pp <- sim_pp_data_biv(lambda_n = 500, abundance = 0.3)

# Compute KAMP expectation
kamp_expec_biv <- kamp_expectation_biv(pp)
print(kamp_expec_biv)

# Compute KAMP variance
kamp_var_biv <- kamp_variance_biv(pp)
print(kamp_var_biv)
```

# Documentation
Link to documentation and vignettes: https://dliao1.github.io/KAMP/

  <!-- badges: start -->
  [![Codecov test coverage](https://codecov.io/gh/dliao1/kperm/graph/badge.svg)](https://app.codecov.io/gh/dliao1/KAMP)
  <!-- badges: end -->
