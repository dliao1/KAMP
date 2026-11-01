# GAMP bivariate expectation

Computes the GAMP (G-function Adjusted for Marked Permutations)
expectation for bivariate point patterns: the exact mean of the
uncorrected cross-type nearest-neighbor G function (distance from each
`mark1` point to its nearest `mark2` point) under random relabeling of
marks, without running any permutations.

Observed G matches
`spatstat.explore::Gcross(..., i = mark1, j = mark2, correction = "none")`.

## Usage

``` r
gamp_expectation_biv(
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

  Vector of radii at which to calculate the GAMP expectation. Defaults
  to c(0, 0.05, 0.075, 0.1, 0.15, 0.2).

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

  The theoretical G under CSR, `1 - exp(-lambda2 * pi * r^2)`

- gamp_csr:

  The adjusted CSR representing the GAMP permuted expectation.

- gamp:

  The difference between observed G and GAMP CSR

## Details

Computes GAMP Expectation for Bivariate Point Patterns

## Examples

``` r
# one image from the ovarian cancer data in the package
ids <- unique(ovarian_df$sample_id)
img <- ovarian_df[ovarian_df$sample_id == ids[4], ]
pp <- spatstat.geom::ppp(img$x, img$y,
                         window = spatstat.geom::owin(range(img$x), range(img$y)),
                         marks = img$phenotype)

gamp_expectation_biv(pp, rvals = c(0, 10, 20, 30, 40),
                     mark1 = "cytotoxic t cell", mark2 = "helper t cell")
#> # A tibble: 5 × 5
#>       r     g theo_csr gamp_csr    gamp
#>   <dbl> <dbl>    <dbl>    <dbl>   <dbl>
#> 1     0 0       0        0      0      
#> 2    10 0.392   0.0449   0.0569 0.335  
#> 3    20 0.582   0.168    0.276  0.306  
#> 4    30 0.667   0.339    0.511  0.156  
#> 5    40 0.711   0.520    0.705  0.00543
```
