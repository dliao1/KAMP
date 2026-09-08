# KAMP univariate expectation

Computes the KAMP (K-function Adjusted for Marked Permutations)
expectation for a given spatial point pattern. This function calculates
Ripley's K using both the traditional Ripley's K method (based on
`Kcross`) and the KAMP-adjusted CSR baseline (based on `Kest`).

The KAMP-adjusted CSR represents a more robust baseline for K (compared
to traditional CSR) that accounts for spatial clustering or
inhomogeneity in a point pattern compared to the traditional CSR
assumption, while avoiding the computational burden of permuting the
point pattern.

Notes:

This function uses the `spatstat` package under the hood, which
automatically uses border correction when the number of points in the
point pattern is more than 3000.

See [`?Kcross`](https://rdrr.io/pkg/spatstat.explore/man/Kcross.html)
and [`?Kest`](https://rdrr.io/pkg/spatstat.explore/man/Kest.html) for
more details on the K calculation methods.

## Usage

``` r
kamp_expectation(
  ppp_obj,
  rvals = c(0, 0.05, 0.075, 0.1, 0.15, 0.2),
  correction = "trans",
  mark1 = "immune"
)
```

## Arguments

- ppp_obj:

  A point pattern object from the `spatstat.geom` package.

- rvals:

  Vector of radii at which to calculate the KAMP expectation. Defaults
  to c(0, 0.05, 0.075, 0.1, 0.15, 0.2).

- correction:

  Type of edge correction method to be used and passed to `Kcross` and
  `Kest`. Defaults to translational edge correction.

- mark1:

  Identifies subset of marked points. Defaults to immune.

## Value

A dataframe with the following columns:

- r:

  The radius at which K was calculated.

- k:

  The observed K value from `Kcross`

- theo_csr:

  The theoretical K under CSR from `Kcross`

- kamp_csr:

  The adjusted CSR representing the permuted expectation.

- kamp:

  The difference between observed K and KAMP CSR

## Details

Compute KAMP Expectation

## Examples

``` r
win <- spatstat.geom::owin(c(0, 1), c(0, 1))
pp <- spatstat.random::rpoispp(lambda = 150, win = win)
marks <- sample(c("immune", "background"), pp$n, replace = TRUE, prob = c(0.4, 0.6))
marked_pp <- spatstat.geom::ppp(pp$x, pp$y, window = win, marks = factor(marks))

result <- kamp_expectation(marked_pp, rvals = c(0, 0.05, 0.1), mark1 = "immune")
print(result)
#> # A tibble: 3 × 5
#>       r      k theo_csr kamp_csr    kamp
#>   <dbl>  <dbl>    <dbl>    <dbl>   <dbl>
#> 1  0    0       0        0       0      
#> 2  0.05 0.0108  0.00785  0.00717 0.00365
#> 3  0.1  0.0333  0.0314   0.0314  0.00188
```
