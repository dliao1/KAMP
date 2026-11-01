# Compute B-KAMP at one relative radius

Standardize within-block KAMP statistics and combine them.

## Usage

``` r
bkamp_helper(
  blocks,
  rvalue,
  univariate = TRUE,
  correction = "trans",
  mark1,
  mark2 = NULL
)
```

## Arguments

- blocks:

  Nonempty list of marked rectangular `ppp` blocks.

- rvalue:

  A single radius In B-KAMP, this is relative to each block's shorter
  side and is multiplied by that length to obtain the absolute radius
  used for KAMP calculations.

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

A one-row tibble with relative radius `r`, combined `z`, and `pvalue`.
