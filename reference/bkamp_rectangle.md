# Identify a valid rectangle using row-wise histograms

Find a count-valid rectangle with the smallest point count using the
points in the available cells.

## Usage

``` r
bkamp_rectangle(
  grid_mat,
  max_aspect_ratio,
  ppp_obj,
  n_min,
  x_coords,
  y_coords,
  ref_wd,
  ref_hi,
  univariate = TRUE
)
```

## Arguments

- grid_mat:

  Binary availability matrix; rows run from top to bottom.

- max_aspect_ratio:

  Maximum longer/shorter side ratio of a candidate.

- ppp_obj:

  A point pattern object of class "ppp" from the spatstat package.

- n_min:

  Target/background minimum counts, or type1/type2/background counts.

- x_coords, y_coords:

  Increasing grid boundary coordinates.

- ref_wd, ref_hi:

  Width and height of a grid cell.

- univariate:

  A logical value indicating whether to compute univariate KAMP (default
  is TRUE).

## Value

The selected rectangle bounds, or NULL if no candidate is valid.
