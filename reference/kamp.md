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
  mark_var,
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

  Column name in `df` containing the marks, when `df` is a data.frame.
  Ignored when `df` is already a `ppp` object.

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

  Percentage that determines how much to thin the amount of points in
  the point pattern object. Default is 0.

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
# Loads required packages
library(spatstat.geom)
#> Loading required package: spatstat.data
#> Loading required package: spatstat.univar
#> spatstat.univar 3.2-0
#> spatstat.geom 3.8-2
library(spatstat.explore)
#> Loading required package: spatstat.random
#> spatstat.random 3.5-1
#> Loading required package: nlme
#> spatstat.explore 3.8-2

# Simulates a simple marked point pattern
set.seed(100)
x_coords <- runif(100)
y_coords <- runif(100)
marks_vec <- sample(c("immune", "background"), 100, replace = TRUE)
win <- owin(c(0,1), c(0,1))
ppp_obj <- ppp(x_coords, y_coords, window = win, marks = factor(marks_vec))

# Defines radius values for K-function estimation (must start at 0 for Kcross/Kest)
r_vals <- seq(0, 0.1, by = 0.01)

# Computes univariate KAMP expectation, passing in a ppp object directly
kamp_result <- kamp(df = ppp_obj,
                    rvals = r_vals,
                    univariate = TRUE,
                    mark1 = "immune")
#> We expect the dataframe to be a single point process. If you have multiple point processes, subset the dataframe by ID and please run KAMP separately for each process.
head(kamp_result)
#> # A tibble: 6 × 5
#>       r       k theo_csr kamp_csr      kamp
#>   <dbl>   <dbl>    <dbl>    <dbl>     <dbl>
#> 1  0    0       0        0         0       
#> 2  0.01 0       0.000314 0.000203 -0.000203
#> 3  0.02 0.00118 0.00126  0.00103   0.000153
#> 4  0.03 0.00358 0.00283  0.00249   0.00109 
#> 5  0.04 0.00478 0.00503  0.00565  -0.000872
#> 6  0.05 0.00724 0.00785  0.0102   -0.00293 

# df can also be a plain data.frame with x/y columns and a marks column,
# identified via mark_var -- kamp() builds the ppp object internally
pts_df <- data.frame(x = x_coords, y = y_coords, cell_type = marks_vec)
kamp_from_df <- kamp(df = pts_df,
                     rvals = r_vals,
                     univariate = TRUE,
                     mark_var = "cell_type",
                     mark1 = "immune")
#> We expect the dataframe to be a single point process. If you have multiple point processes, subset the dataframe by ID and please run KAMP separately for each process.
head(kamp_from_df)
#> # A tibble: 6 × 5
#>       r        k theo_csr kamp_csr      kamp
#>   <dbl>    <dbl>    <dbl>    <dbl>     <dbl>
#> 1  0    0        0        0         0       
#> 2  0.01 0        0.000314 0.000167 -0.000167
#> 3  0.02 0.000974 0.00126  0.000848  0.000126
#> 4  0.03 0.00296  0.00283  0.00206   0.000902
#> 5  0.04 0.00396  0.00503  0.00468  -0.000724
#> 6  0.05 0.00598  0.00785  0.00842  -0.00244 

# Compute univariate KAMP expectation with thinning
kamp_thin <- kamp(df = ppp_obj,
                  rvals = r_vals,
                  univariate = TRUE,
                  mark1 = "immune",
                  thin = TRUE,
                  p_thin = 0.3)
#> We expect the dataframe to be a single point process. If you have multiple point processes, subset the dataframe by ID and please run KAMP separately for each process.
head(kamp_thin)
#> # A tibble: 6 × 5
#>       r       k theo_csr kamp_csr      kamp
#>   <dbl>   <dbl>    <dbl>    <dbl>     <dbl>
#> 1  0    0       0        0         0       
#> 2  0.01 0       0.000314 0.000482 -0.000482
#> 3  0.02 0.00289 0.00126  0.00195   0.000941
#> 4  0.03 0.00586 0.00283  0.00345   0.00241 
#> 5  0.04 0.00586 0.00503  0.00646  -0.000603
#> 6  0.05 0.00888 0.00785  0.0106   -0.00168 

# Use real data from VectraPolarisData in package
data(ovarian_df)
first_sample <- unique(ovarian_df$sample_id)[1]
ov_df <- subset(ovarian_df, sample_id == first_sample)
win <- convexhull.xy(ov_df$x, ov_df$y)
ppp_real <- ppp(ov_df$x, ov_df$y, window = win, marks = ov_df$immune)
kamp_real <- kamp(df = ppp_real,
                  rvals = seq(0, 0.1, 0.01),
                  univariate = TRUE,
                  mark1 = "immune")
#> We expect the dataframe to be a single point process. If you have multiple point processes, subset the dataframe by ID and please run KAMP separately for each process.
head(kamp_real)
#> # A tibble: 6 × 5
#>       r     k theo_csr kamp_csr  kamp
#>   <dbl> <dbl>    <dbl>    <dbl> <dbl>
#> 1  0        0 0               0     0
#> 2  0.01     0 0.000314        0     0
#> 3  0.02     0 0.00126         0     0
#> 4  0.03     0 0.00283         0     0
#> 5  0.04     0 0.00503         0     0
#> 6  0.05     0 0.00785         0     0
```
