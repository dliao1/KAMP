# Complement retained blocks with adjacent residual cells

Greedily expand blocks while keeping them rectangular and disjoint.

## Usage

``` r
bkamp_complement(grid, success_blocks, assigned, base_tol)
```

## Arguments

- grid:

  List of grid cells carrying temporary blocking marks.

- success_blocks:

  Retained block records, each containing `ppp` and cell `ids`.

- assigned:

  Logical vector identifying cells already assigned to blocks.

- base_tol:

  Numerical tolerance for adjacency and rectangle checks.

## Value

A list containing updated `success_blocks` and `assigned`.
