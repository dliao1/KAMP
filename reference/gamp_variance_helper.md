# GAMP variance helper

Exact permutation expectation and variance of G, with the pair sums done
in Rcpp (`gamp_sums`) from
[`spatstat.geom::closepairs()`](https://rdrr.io/pkg/spatstat.geom/man/closepairs.html).

## Usage

``` r
gamp_variance_helper(ppp_obj, rvals, mark1, mark2 = NULL, correction = "none")
```

## Arguments

- ppp_obj:

  A point pattern object "ppp" from the `spatstat` package.

- rvals:

  Vector of radii.

- mark1:

  Mark of the points G is computed for (e.g. immune).

- mark2:

  Mark of the nearest-neighbor points for bivariate G, or `NULL` for
  univariate G.

- correction:

  Type of edge correction. Only "none" is supported.

## Value

A dataframe with columns `r`, `g`, `theo_csr`, `gamp_csr`, `gamp`,
`var`, `pvalue`.
