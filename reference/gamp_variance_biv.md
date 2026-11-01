# GAMP bivariate variance

Computes the GAMP (G-function Adjusted for Marked Permutations) variance
for bivariate point patterns: the exact variance of the uncorrected
cross-type nearest-neighbor G function (distance from each `mark1` point
to its nearest `mark2` point) under random relabeling of marks. Also
returns the GAMP expectation and p-value.

The pairwise sums are done in Rcpp (`gamp_sums`) from
[`spatstat.geom::closepairs()`](https://rdrr.io/pkg/spatstat.geom/man/closepairs.html),
so only pairs within `max(rvals)` are touched.

## Usage

``` r
gamp_variance_biv(
  ppp_obj,
  rvals = c(0, 0.05, 0.075, 0.1, 0.15, 0.2),
  correction = "none",
  mark1 = "immune1",
  mark2 = "immune2"
)
```

## Arguments

- ppp_obj:

  A point pattern object from the `spatstat.geom` package.

- rvals:

  A vector of radii at which to calculate the GAMP variance. Defaults to
  c(0, 0.05, 0.075, 0.1, 0.15, 0.2).

- correction:

  Type of edge correction. Only "none" (no edge correction) is currently
  supported.

- mark1:

  Variable used to mark the points in the point pattern object for the
  first type. Default is "immune1".

- mark2:

  Variable used to mark the points in the point pattern object for the
  second type. Default is "immune2".

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

Computes GAMP Variance for Bivariate Point Patterns

## Examples

``` r
# one image from the ovarian cancer data in the package
ids <- unique(ovarian_df$sample_id)
img <- ovarian_df[ovarian_df$sample_id == ids[4], ]
pp <- spatstat.geom::ppp(img$x, img$y,
                         window = spatstat.geom::owin(range(img$x), range(img$y)),
                         marks = img$phenotype)

gamp_variance_biv(pp, rvals = c(10, 20, 30, 40),
                  mark1 = "cytotoxic t cell", mark2 = "helper t cell")
#> # A tibble: 4 × 7
#>       r     g theo_csr gamp_csr    gamp      var    pvalue
#>   <dbl> <dbl>    <dbl>    <dbl>   <dbl>    <dbl>     <dbl>
#> 1    10 0.392   0.0449   0.0569 0.335   0.000196 1.08e-126
#> 2    20 0.582   0.168    0.276  0.306   0.000739 1.12e- 29
#> 3    30 0.667   0.339    0.511  0.156   0.000976 3.20e-  7
#> 4    40 0.711   0.520    0.705  0.00543 0.000886 4.28e-  1
```
