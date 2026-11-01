# GAMP expectation helper

Exact permutation expectation of G, as in PermFree-NearestNeighborG:
`E(G) = 1 - (1/N) * sum_p choose(N - 1 - n_p, t) / choose(N - 1, t)`,
where n_p(r) is the number of neighbors of point p within r and t is the
number of mark2 points among the other N - 1 points (m1 - 1 univariate,
m2 bivariate).

## Usage

``` r
gamp_expectation_helper(
  ppp_obj,
  rvals,
  mark1,
  mark2 = NULL,
  correction = "none"
)
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

A dataframe with columns `r`, `g`, `theo_csr`, `gamp_csr`, `gamp`.
