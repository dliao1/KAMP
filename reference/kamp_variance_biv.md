# KAMP bivariate variance

Computes the KAMP (K-function Adjusted for Marked Permutations) variance
for bivariate point patterns. This function calculates Ripley's K using
both the traditional Ripley's K method and the KAMP-adjusted CSR
baseline by using a matrix-based implementation.

The KAMP-adjusted CSR represents a more realistic baseline for K
(compared to traditional CSR) that accounts for spatial clustering or
inhomogeneity in a point pattern compared to the traditional CSR
assumption, while avoiding the computational burden of permuting the
point pattern.

Note: This function implements a slower, matrix-based implementation of
the KAMP variance. It is a wrapper around the `kamp_variance_biv_helper`
function that calculates the KAMP variance at one radius and maps it
over a vector of radii.

## Usage

``` r
kamp_variance_biv(
  ppp_obj,
  rvals = c(0, 0.05, 0.075, 0.1, 0.15, 0.2),
  correction = "trans",
  mark1 = "immune1",
  mark2 = "immune2"
)
```

## Arguments

- ppp_obj:

  A point pattern object from the `spatstat.geom` package.

- correction:

  Type of edge correction. Defaults to translational.

- rvec:

  A vector of radii at which to calculate the KAMP expectation. Defaults
  to c(0, 0.05, 0.075, 0.1, 0.15, 0.2).

- markvar1:

  Variable used to mark the points in the point pattern object for the
  first type. Default is "immune1".

- markvar2:

  Variable used to mark the points in the point pattern object for the
  second type. Default is "immune2".

## Value

A dataframe with the following columns:

- r:

  The radius at which K was calculated.

- k:

  The observed K value

- theo_csr:

  The theoretical K under CSR

- kamp_csr:

  The adjusted CSR representing the KAMP permuted expectation.

- var:

  Variance of K under the permutation null distribution

- z:

  Z statistic, calculated by normalizing K using the formula: (K -
  KAMP)/sqrt(var)

- pval:

  P-value, calculated using the formula: pnorm(-z)

## Details

Computes KAMP Variance for Bivariate Point Patterns
