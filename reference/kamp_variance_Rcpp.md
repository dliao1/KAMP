# KAMP univariate variance - Rcpp Version

Computes the KAMP (K-function Adjusted for Marked Permutations) variance
for a given spatial point pattern. Also returns the KAMP expectation,
z-statistic, and p-value. Calculations are performed across a vector of
radii in a single call.

Note: this a matrix-based implementation of the KAMP variance. It relies
on an Rcpp function `kamp_pair_sums` to efficiently compute the
necessary pairwise sums for the KAMP variance calculation, as well as
the `spatstat` package for handling point pattern distances and edge
corrections.

## Usage

``` r
kamp_variance_Rcpp(
  ppp_obj,
  rvals = c(0, 0.05, 0.075, 0.1, 0.15, 0.2),
  correction = "trans",
  mark1 = "immune"
)
```

## Arguments

- ppp_obj:

  A point pattern object of class "ppp" from the spatstat package.

- rvals:

  A vector of radii at which to calculate the KAMP variance.

- correction:

  Type of edge correction. Defaults to translational.

- mark1:

  The variable used to mark the points in the point pattern object.
  Default is "immune".

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

- kamp:

  The difference between observed K and KAMP CSR

- var:

  Variance of K under the permutation null distribution

- pval:

  P-value, calculated using the formula: pnorm(-z)

## Details

KAMP Variance
