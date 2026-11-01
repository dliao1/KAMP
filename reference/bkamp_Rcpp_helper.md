# Compute B-KAMP across relative radii using KAMP's Rcpp functions

Calculate all requested radii once per block and combine block
statistics.

## Usage

``` r
bkamp_Rcpp_helper(
  blocks,
  rvals = c(0, 0.05, 0.075, 0.1, 0.15, 0.2),
  univariate = TRUE,
  correction = "trans",
  mark1,
  mark2 = NULL
)
```

## Arguments

- blocks:

  Nonempty list of marked rectangular `ppp` blocks.

- rvals:

  A vector of radii at which to calculate the KAMP variance. In B-KAMP,
  these are relative to each block's shorter side and are multiplied by
  that length to obtain the absolute radii used for KAMP calculations.

- univariate:

  A logical value indicating whether to compute univariate KAMP (default
  is TRUE).

- correction:

  Type of edge correction. Defaults to translational.

- mark1:

  Variable used to mark the points in the point pattern object for the
  first type.

- mark2:

  Variable used to mark the points in the point pattern object for the
  second type (optional, only used if `univariate` is FALSE).

## Value

A tibble with relative radius `r`, combined `z`, and `pvalue`.
