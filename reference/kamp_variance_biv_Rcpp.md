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
# one image from the ovarian cancer data in the package
ids <- unique(ovarian_df$sample_id)
img <- ovarian_df[ovarian_df$sample_id == ids[4], ]
pp <- spatstat.geom::ppp(img$x, img$y,
                         window = spatstat.geom::owin(range(img$x), range(img$y)),
                         marks = img$phenotype)

kamp_variance_biv_Rcpp(pp, rvals = c(25, 50, 75, 100),
                       mark1 = "cytotoxic t cell", mark2 = "helper t cell")
#> # A tibble: 4 × 7
#>       r       k theo_csr kamp_csr    kamp      var pvalue
#>   <dbl>   <dbl>    <dbl>    <dbl>   <dbl>    <dbl>  <dbl>
#> 1    25  20628.    1963.    3541.  17087.   92171.      0
#> 2    50  65542.    7854.   13628.  51913.  416484.      0
#> 3    75 126400.   17671.   29464.  96936. 1136196.      0
#> 4   100 194327.   31416.   50719. 143608. 2532450.      0
```
