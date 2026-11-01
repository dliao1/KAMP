# KAMP

This function computes the KAMP expectation and variance for a given
point pattern. Calculates Ripley's K using both the traditional Ripley's
K method (based on `Kcross`) and the KAMP-adjusted CSR baseline (based
on `Kest`).

The KAMP-adjusted CSR represents a more robust baseline for K that
accounts for spatial clustering or inhomogeneity in a point pattern
compared to the traditional CSR assumption, while avoiding the
computational burden of permuting the point pattern.

For expectation, this function uses the `spatstat` package under the
hood, which automatically uses border correction when the number of
points in the point pattern is more than 3000.

For variance, this function uses the Rcpp-accelerated implementation.

See `?Kcross` and `?Kest` for more details on the K calculation methods.

## Usage

``` r
kamp(
  df,
  rvals,
  univariate = TRUE,
  mark_var = NULL,
  mark1,
  mark2 = NULL,
  variance = FALSE,
  correction = "trans",
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

  A vector of distances at which to compute the KAMP expectation and
  variance.

- univariate:

  A logical value indicating whether to compute univariate KAMP (default
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

  A logical value indicating whether to compute the variance of KAMP
  (default is FALSE).

- correction:

  Type of edge correction. Defaults to translational.

- thin:

  A logical value indicating whether to thin the point pattern before
  computing KAMP (default is FALSE), called KAMP-lite.

- p_thin:

  Proportion of points to remove when `thin = TRUE`. Default is 0.5.

- ...:

  Additional arguments passed to the underlying functions.

## Value

A dataframe with the following columns:

- r:

  The radius at which K was calculated.

- k:

  The observed K value

- theo_csr:

  The theoretical K under CSR

- kamp_csr:

  The adjusted CSR representing the KAMP expectation.

- kamp:

  The difference between observed K and KAMP CSR

- var:

  If variance = TRUE, variance of K under the KAMP null distribution

- pval:

  If variance = TRUE, p-value that estimates the probability of
  observing a deviation from the expected KAMP-adjusted value as large
  or larger than the one observed, under the null hypothesis of CSR).
  Calculated using a normal approximation.

## Details

Computes KAMP expectation and variance

## Examples

``` r
# one image from the ovarian cancer data in the package
ids <- unique(ovarian_df$sample_id)
img <- ovarian_df[ovarian_df$sample_id == ids[4], ]
r_vals <- c(0, 25, 50, 75, 100)

# univariate, marks column given by mark_var
kamp(img, rvals = r_vals, mark_var = "immune", mark1 = "immune")
#> We expect the dataframe to be a single point process. If you have multiple point processes, subset the dataframe by ID and please run KAMP separately for each process.
#> # A tibble: 5 × 5
#>       r       k theo_csr kamp_csr   kamp
#>   <dbl>   <dbl>    <dbl>    <dbl>  <dbl>
#> 1     0      0        0        0      0 
#> 2    25  14035.    1963.    2801. 11234.
#> 3    50  43815.    7854.   10780. 33034.
#> 4    75  82001.   17671.   23309. 58692.
#> 5   100 126106.   31416.   40133. 85973.

# with variance
kamp(img, rvals = r_vals, mark_var = "immune", mark1 = "immune", variance = TRUE)
#> We expect the dataframe to be a single point process. If you have multiple point processes, subset the dataframe by ID and please run KAMP separately for each process.
#> # A tibble: 5 × 7
#>       r       k theo_csr kamp_csr   kamp     var pvalue
#>   <dbl>   <dbl>    <dbl>    <dbl>  <dbl>   <dbl>  <dbl>
#> 1     0      0        0        0      0       0     NaN
#> 2    25  14008.    1963.    2796. 11212.  10961.      0
#> 3    50  43624.    7854.   10730. 32894.  62689.      0
#> 4    75  81395.   17671.   23120. 58275. 214741.      0
#> 5   100 124703.   31416.   39641. 85062. 562437.      0

# bivariate
kamp(img, rvals = r_vals, univariate = FALSE, mark_var = "phenotype",
     mark1 = "cytotoxic t cell", mark2 = "helper t cell")
#> We expect the dataframe to be a single point process. If you have multiple point processes, subset the dataframe by ID and please run KAMP separately for each process.
#> # A tibble: 5 × 5
#>       r       k theo_csr kamp_csr    kamp
#>   <dbl>   <dbl>    <dbl>    <dbl>   <dbl>
#> 1     0      0        0        0       0 
#> 2    25  16318.    1963.    2801.  13517.
#> 3    50  51845.    7854.   10780.  41064.
#> 4    75 100003.   17671.   23309.  76694.
#> 5   100 153760.   31416.   40133. 113626.

# KAMP-lite: randomly remove 30% of points first
kamp(img, rvals = r_vals, mark_var = "immune", mark1 = "immune",
     thin = TRUE, p_thin = 0.3)
#> We expect the dataframe to be a single point process. If you have multiple point processes, subset the dataframe by ID and please run KAMP separately for each process.
#> # A tibble: 5 × 5
#>       r       k theo_csr kamp_csr   kamp
#>   <dbl>   <dbl>    <dbl>    <dbl>  <dbl>
#> 1     0      0        0        0      0 
#> 2    25  14916.    1963.    2809. 12107.
#> 3    50  45966.    7854.   10885. 35081.
#> 4    75  84001.   17671.   23594. 60407.
#> 5   100 127497.   31416.   40569. 86927.
```
