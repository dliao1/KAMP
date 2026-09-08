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

For variance, this function utlizes a matrix-based implementation.

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
  background = NULL,
  Rcpp = FALSE,
  ...
)
```

## Arguments

- rvals:

  A vector of distances at which to compute the KAMP expectation and
  variance.

- univariate:

  A logical value indicating whether to compute univariate KAMP (default
  is TRUE).

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

- background:

  Variable used to define the background for the point pattern object.

- Rcpp:

  A logical value indicating whether to use Rcpp implementation (default
  is FALSE).

- ...:

  Additional arguments passed to the underlying functions.

- ppp_obj:

  A point pattern object of class `ppp` from the `spatstat` package.

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
ppp_obj <- ppp(x_coords, y_coords, window = win, marks = marks_vec)

# Defines radius values for K-function estimation
r_vals <- seq(0.01, 0.1, by = 0.01)

# Computes univariate KAMP expectation
kamp_result <- kamp(ppp_obj = ppp_obj,
                    rvals = r_vals,
                    univariate = TRUE,
                    mark1 = "immune")
#> Error in kamp(ppp_obj = ppp_obj, rvals = r_vals, univariate = TRUE, mark1 = "immune"): argument "df" is missing, with no default
head(kamp_result)
#> Error: object 'kamp_result' not found

# Compute univariate KAMP expectation with thinning
kamp_thin <- kamp(ppp_obj = ppp_obj,
                  rvals = r_vals,
                  univariate = TRUE,
                  mark1 = "immune",
                  thin = TRUE,
                  p_thin = 0.3)
#> Error in kamp(ppp_obj = ppp_obj, rvals = r_vals, univariate = TRUE, mark1 = "immune",     thin = TRUE, p_thin = 0.3): argument "df" is missing, with no default
head(kamp_thin)
#> Error: object 'kamp_thin' not found

# Use real data from VectraPolarisData in package
data(ovarian_df)
sample_id <- unique(ovarian_df$sample_id)[1]
ov_df <- subset(ovarian_df, sample_id == sample_id)
win <- convexhull.xy(ov_df$x, ov_df$y)
ppp_real <- ppp(ov_df$x, ov_df$y, window = win, marks = ov_df$immune)
kamp_real <- kamp(ppp_obj = ppp_real,
                  rvals = seq(0.01, 0.1, 0.01),
                  univariate = TRUE,
                  mark1 = "immune")
#> Error in kamp(ppp_obj = ppp_real, rvals = seq(0.01, 0.1, 0.01), univariate = TRUE,     mark1 = "immune"): argument "df" is missing, with no default
head(kamp_real)
#> Error: object 'kamp_real' not found
```
