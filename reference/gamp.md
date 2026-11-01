# GAMP

This function computes the GAMP expectation and variance for a given
point pattern. GAMP is the nearest-neighbor G analogue of KAMP: it gives
the exact mean (and optionally variance) of the uncorrected G function
under random relabeling of marks, i.e. the quantities a permutation test
would estimate, without running any permutations.

Univariate G is the nearest-neighbor distance distribution of the
`mark1` points (matches `Gest(..., correction = "none")`). Bivariate G
is the distribution of the distance from each `mark1` point to its
nearest `mark2` point (matches `Gcross(..., correction = "none")`). Only
no edge correction is supported.

Neighbor counts are computed from
[`spatstat.geom::closepairs()`](https://rdrr.io/pkg/spatstat.geom/man/closepairs.html),
so only pairs within `max(rvals)` are touched. For variance, this
function uses the Rcpp-accelerated implementation.

## Usage

``` r
gamp(
  df,
  rvals,
  univariate = TRUE,
  mark_var = NULL,
  mark1,
  mark2 = NULL,
  variance = FALSE,
  correction = "none",
  thin = FALSE,
  p_thin = 0.5,
  ...
)
```

## Arguments

- df:

  Either a point pattern object of class `ppp` from the `spatstat`
  package, or a data.frame with `x`/`y` columns and a marks column named
  by `mark_var` (used to build a `ppp` internally via a convex-hull
  window). Should contain one point process at a time.

- rvals:

  A vector of distances at which to compute the GAMP expectation and
  variance.

- univariate:

  A logical value indicating whether to compute univariate GAMP (default
  is TRUE).

- mark_var:

  Column name in `df` containing the marks. Required when `df` is a
  data.frame (an error is thrown if it is left as `NULL`); ignored when
  `df` is already a `ppp` object.

- mark1:

  Variable used to mark the points in the point pattern object for the
  first type.

- mark2:

  Variable used to mark the points in the point pattern object for the
  second type (optional, only used if `univariate` is FALSE).

- variance:

  A logical value indicating whether to compute the variance of GAMP
  (default is FALSE).

- correction:

  Type of edge correction. Only "none" (no edge correction) is currently
  supported.

- thin:

  A logical value indicating whether to thin the point pattern before
  computing GAMP (default is FALSE).

- p_thin:

  Proportion of points to remove when `thin = TRUE`. Default is 0.5.

- ...:

  Additional arguments passed to the underlying functions.

## Value

A dataframe with the following columns:

- r:

  The radius at which G was calculated.

- g:

  The observed G value

- theo_csr:

  The theoretical G under CSR, `1 - exp(-lambda * pi * r^2)`, where
  lambda is the intensity of the `mark1` (univariate) or `mark2`
  (bivariate) points

- gamp_csr:

  The adjusted CSR representing the GAMP expectation.

- gamp:

  The difference between observed G and GAMP CSR

- var:

  If variance = TRUE, variance of G under the GAMP null distribution

- pvalue:

  If variance = TRUE, p-value for G being larger than expected
  (clustering), calculated using a normal approximation: pnorm(-z).

## Details

Computes GAMP expectation and variance

## Examples

``` r
# one image from the ovarian cancer data in the package
ids <- unique(ovarian_df$sample_id)
img <- ovarian_df[ovarian_df$sample_id == ids[4], ]
r_vals <- c(0, 10, 20, 30, 40)

# univariate with variance
gamp(img, rvals = r_vals, mark_var = "immune", mark1 = "immune", variance = TRUE)
#> We expect the dataframe to be a single point process. If you have multiple point processes, subset the dataframe by ID and please run KAMP separately for each process.
#> # A tibble: 5 × 7
#>       r     g theo_csr gamp_csr    gamp       var    pvalue
#>   <dbl> <dbl>    <dbl>    <dbl>   <dbl>     <dbl>     <dbl>
#> 1     0 0        0        0      0      0         1   e+  0
#> 2    10 0.668    0.150    0.150  0.517  0.000219  2.47e-267
#> 3    20 0.869    0.478    0.594  0.275  0.000298  1.65e- 57
#> 4    30 0.920    0.768    0.854  0.0656 0.000120  1.11e-  9
#> 5    40 0.943    0.926    0.954 -0.0108 0.0000391 9.58e-  1

# bivariate
gamp(img, rvals = r_vals, univariate = FALSE, mark_var = "phenotype",
     mark1 = "cytotoxic t cell", mark2 = "helper t cell")
#> We expect the dataframe to be a single point process. If you have multiple point processes, subset the dataframe by ID and please run KAMP separately for each process.
#> # A tibble: 5 × 5
#>       r     g theo_csr gamp_csr    gamp
#>   <dbl> <dbl>    <dbl>    <dbl>   <dbl>
#> 1     0 0       0        0      0      
#> 2    10 0.392   0.0564   0.0569 0.335  
#> 3    20 0.582   0.207    0.276  0.306  
#> 4    30 0.667   0.407    0.511  0.156  
#> 5    40 0.711   0.605    0.705  0.00543
```
