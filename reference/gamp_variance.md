# GAMP univariate variance

Computes the GAMP (G-function Adjusted for Marked Permutations) variance
for a given spatial point pattern: the exact variance of the uncorrected
nearest-neighbor G function of the `mark1` points under random
relabeling of marks. Also returns the GAMP expectation and p-value.

The pairwise sums are done in Rcpp (`gamp_sums`) from
[`spatstat.geom::closepairs()`](https://rdrr.io/pkg/spatstat.geom/man/closepairs.html),
so only pairs within `max(rvals)` are touched.

## Usage

``` r
gamp_variance(
  ppp_obj,
  rvals = c(0, 0.05, 0.075, 0.1, 0.15, 0.2),
  correction = "none",
  mark1 = "immune"
)
```

## Arguments

- ppp_obj:

  A point pattern object of class "ppp" from the spatstat package.

- rvals:

  A vector of radii at which to calculate the GAMP variance.

- correction:

  Type of edge correction. Only "none" (no edge correction) is currently
  supported.

- mark1:

  The variable used to mark the points in the point pattern object.
  Default is "immune".

## Value

A dataframe with the following columns:

- r:

  The radius at which G was calculated.

- g:

  The observed G value

- theo_csr:

  The theoretical G under CSR

- gamp_csr:

  The adjusted CSR representing the GAMP permuted expectation.

- gamp:

  The difference between observed G and GAMP CSR

- var:

  Variance of G under the permutation null distribution

- pvalue:

  P-value, calculated using the formula: pnorm(-z). Set to 1 when the
  variance is 0.

## Details

GAMP Variance

## Examples

``` r
# one image from the ovarian cancer data in the package
ids <- unique(ovarian_df$sample_id)
img <- ovarian_df[ovarian_df$sample_id == ids[4], ]
pp <- spatstat.geom::ppp(img$x, img$y,
                         window = spatstat.geom::owin(range(img$x), range(img$y)),
                         marks = img$immune)

gamp_variance(pp, rvals = c(10, 20, 30, 40), mark1 = "immune")
#> # A tibble: 4 × 7
#>       r     g theo_csr gamp_csr    gamp       var    pvalue
#>   <dbl> <dbl>    <dbl>    <dbl>   <dbl>     <dbl>     <dbl>
#> 1    10 0.668    0.121    0.150  0.517  0.000219  2.47e-267
#> 2    20 0.869    0.402    0.594  0.275  0.000298  1.65e- 57
#> 3    30 0.920    0.685    0.854  0.0656 0.000120  1.11e-  9
#> 4    40 0.943    0.872    0.954 -0.0108 0.0000391 9.58e-  1
```
