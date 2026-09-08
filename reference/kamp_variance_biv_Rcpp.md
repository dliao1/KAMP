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
kamp_variance_biv_Rcpp(
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

- rvals:

  A vector of radii at which to calculate the KAMP expectation. Defaults
  to c(0, 0.05, 0.075, 0.1, 0.15, 0.2).

- correction:

  Type of edge correction. Defaults to translational.

- mark1:

  Variable used to mark the points in the point pattern object for the
  first type. Default is "immune1".

- mark2:

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

## Examples

``` r
win <- spatstat.geom::owin(c(0, 1), c(0, 1))
pp <- spatstat.random::rpoispp(lambda = 150, win = win)
mark_labels <- c("immune1", "immune2", "background")
marks <- sample(mark_labels, pp$n, replace = TRUE, prob = c(0.3, 0.3, 0.4))
marked_pp <- spatstat.geom::ppp(pp$x, pp$y, window = win, marks = factor(marks))

result <- kamp_variance_biv_Rcpp(marked_pp, rvals = c(0.05, 0.1),
                                 mark1 = "immune1", mark2 = "immune2")
print(result)
#> # A tibble: 2 × 7
#>       r       k theo_csr kamp_csr      kamp        var pvalue
#>   <dbl>   <dbl>    <dbl>    <dbl>     <dbl>      <dbl>  <dbl>
#> 1  0.05 0.00876  0.00785  0.00801  0.000758 0.00000424  0.357
#> 2  0.1  0.0292   0.0314   0.0298  -0.000569 0.0000158   0.557
```
