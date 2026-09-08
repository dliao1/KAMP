# KAMP univariate variance - Rcpp Version

Computes the KAMP (K-function Adjusted for Marked Permutations) variance
for a given spatial point pattern. Also returns the KAMP expectation,
z-statistic, and p-value. Calculations are performed across a vector of
radii in a single call.

Note: this a matrix-based implementation of the KAMP variance. It relies
on Rcpp functions (`kamp_pair_sums_trans_iso`, `kamp_pair_sums_bord`) to
do the pairwise sums, and `spatstat` for point pattern distances/edge
corrections. Translational, isotropic, border, and no ("none") edge
corrections are all supported.

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

## Examples

``` r
win <- spatstat.geom::owin(c(0, 1), c(0, 1))
pp <- spatstat.random::rpoispp(lambda = 150, win = win)
marks <- sample(c("immune", "background"), pp$n, replace = TRUE, prob = c(0.4, 0.6))
marked_pp <- spatstat.geom::ppp(pp$x, pp$y, window = win, marks = factor(marks))

result <- kamp_variance_Rcpp(marked_pp, rvals = c(0.05, 0.1, 0.15), mark1 = "immune")
print(result)
#> # A tibble: 3 × 7
#>       r       k theo_csr kamp_csr      kamp        var pvalue
#>   <dbl>   <dbl>    <dbl>    <dbl>     <dbl>      <dbl>  <dbl>
#> 1  0.05 0.00881  0.00785  0.00790  0.000911 0.00000226  0.272
#> 2  0.1  0.0381   0.0314   0.0348   0.00336  0.0000130   0.176
#> 3  0.15 0.0753   0.0707   0.0755  -0.000179 0.0000379   0.512
```
