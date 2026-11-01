# Execute an adaptive block identification phase

Repeatedly extract valid rectangles under one aspect ratio bound.

## Usage

``` r
bkamp_phase(
  max_aspect_ratio,
  grid,
  n_min,
  success_blocks,
  assigned,
  ref_wd,
  ref_hi,
  base_tol,
  univariate = TRUE
)
```

## Arguments

- max_aspect_ratio:

  Aspect ratio bound for the current phase.

- grid:

  List of grid cells carrying temporary blocking marks.

- n_min:

  Minimum counts required in each candidate block.

- success_blocks:

  Existing block records, each containing `ppp` and cell `ids`.

- assigned:

  Logical vector identifying cells already assigned to blocks.

- ref_wd, ref_hi:

  Width and height of a grid cell.

- base_tol:

  Numerical tolerance for matching grid boundaries.

- univariate:

  A logical value indicating whether to compute univariate KAMP (default
  is TRUE).

## Value

A list containing updated `success_blocks` and `assigned`.
