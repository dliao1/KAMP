# Construct the adaptive regular grid satisfying the maximum cell count

Refine a regular grid over the rectangular observation window until
every cell satisfies `split_max`.

## Usage

``` r
bkamp_grid(ppp_obj, split_max, max_grid_dim = 100L)
```

## Arguments

- ppp_obj:

  A marked `ppp` with a rectangular observation window.

- split_max:

  Maximum point count per cell, controlling grid refinement.

- max_grid_dim:

  Maximum number of cells per axis to search. Defaults to 100.

## Value

A list of marked `ppp` grid cells with distinct point ownership.
