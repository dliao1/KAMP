# KAMP bivariate expectation

Computes the KAMP (K-function Adjusted for Marked Permutations)
expectation for bivariate point patterns. This function calculates
Ripley's K using both the traditional Ripley's K method (based on
`Kcross`) and the KAMP-adjusted CSR baseline (based on `Kest`).

The KAMP-adjusted CSR represents a more robust baseline for K (compared
to traditional CSR) that accounts for spatial clustering or
inhomogeneity in a point pattern compared to the traditional CSR
assumption, while avoiding the computational burden of permuting the
point pattern.

Note: This function uses the `spatstat` package under the hood. See
[`?Kcross`](https://rdrr.io/pkg/spatstat.explore/man/Kcross.html) and
[`?Kest`](https://rdrr.io/pkg/spatstat.explore/man/Kest.html) for more
details on the K calculation methods.

## Usage

``` r
kamp_expectation_biv(
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

- correction:

  Type of edge correction. Defaults to translational.

- rvec:

  Vector of radii at which to calculate the KAMP expectation. Defaults
  to c(0, 0.05, 0.075, 0.1, 0.15, 0.2).

- markvar1:

  Variable used to mark the points in the point pattern object for the
  first type. Default is "immune1".

- markvar2:

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

- kamp:

  The difference between observed K and KAMP CSR

## Details

Computes KAMP Expectation for Bivariate Point Patterns

## Examples

``` r

if (requireNamespace("spatstat.geom", quietly = TRUE) &&
    requireNamespace("spatstat.random", quietly = TRUE)) {

  # Simulates a window
  win <- spatstat.geom::owin(c(0, 1), c(0, 1))

  # Generates 200 total points
  pp <- spatstat.random::rpoispp(lambda = 200, win = win)

  # Assigns three marks: immune1, immune2, and background
  marks <- sample(c("immune1", "immune2", "background"), pp$n, replace = TRUE, prob = c(0.3, 0.3, 0.4))

  # Creates marked point pattern
  marked_pp <- spatstat.geom::ppp(pp$x, pp$y, window = win, marks = factor(marks))

  # Computes KAMP expectation
  result <- kamp_expectation_biv(marked_pp, markvar1 = "immune1", markvar2 = "immune2")
  print(result)
}
#> Error in kamp_expectation_biv(marked_pp, markvar1 = "immune1", markvar2 = "immune2"): unused arguments (markvar1 = "immune1", markvar2 = "immune2")
```
