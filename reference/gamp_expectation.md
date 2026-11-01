# GAMP univariate expectation

Computes the GAMP (G-function Adjusted for Marked Permutations)
expectation for a given spatial point pattern: the exact mean of the
uncorrected nearest-neighbor G function of the `mark1` points under
random relabeling of marks, without running any permutations.

Observed G is the uncorrected nearest-neighbor distance distribution of
the `mark1` points, matching
`spatstat.explore::Gest(..., correction = "none")`.

## Usage

``` r
gamp_expectation(
  ppp_obj,
  rvals = c(0, 0.05, 0.075, 0.1, 0.15, 0.2),
  correction = "none",
  mark1 = "immune"
)
```

## Arguments

- ppp_obj:

  A point pattern object from the `spatstat.geom` package.

- rvals:

  Vector of radii at which to calculate the GAMP expectation. Defaults
  to c(0, 0.05, 0.075, 0.1, 0.15, 0.2).

- correction:

  Type of edge correction. Only "none" (no edge correction) is currently
  supported.

- mark1:

  Identifies subset of marked points. Defaults to immune.

## Value

A dataframe with the following columns:

- r:

  The radius at which G was calculated.

- g:

  The observed G value

- theo_csr:

  The theoretical G under CSR, `1 - exp(-lambda * pi * r^2)`

- gamp_csr:

  The adjusted CSR representing the GAMP permuted expectation.

- gamp:

  The difference between observed G and GAMP CSR

## Details

Compute GAMP Expectation

## Examples

``` r
# one image from the ovarian cancer data in the package
ids <- unique(ovarian_df$sample_id)
img <- ovarian_df[ovarian_df$sample_id == ids[4], ]
pp <- spatstat.geom::ppp(img$x, img$y,
                         window = spatstat.geom::owin(range(img$x), range(img$y)),
                         marks = img$immune)

gamp_expectation(pp, rvals = c(0, 10, 20, 30, 40), mark1 = "immune")
#> # A tibble: 5 × 5
#>       r     g theo_csr gamp_csr    gamp
#>   <dbl> <dbl>    <dbl>    <dbl>   <dbl>
#> 1     0 0        0        0      0     
#> 2    10 0.668    0.121    0.150  0.517 
#> 3    20 0.869    0.402    0.594  0.275 
#> 4    30 0.920    0.685    0.854  0.0656
#> 5    40 0.943    0.872    0.954 -0.0108
```
