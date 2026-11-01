# Adaptive spatial blocking for B-KAMP

Partition a marked pattern into rectangular blocks using grid
representation, two block identification phases, and block complement.

## Usage

``` r
bkamp_block(
  ppp_obj,
  univariate = TRUE,
  mark1,
  mark2 = NULL,
  split_max = NULL,
  n_min = NULL,
  rho1 = NULL,
  rho2 = Inf,
  max_grid_dim = 100L
)
```

## Arguments

- ppp_obj:

  A marked `ppp` with a rectangular observation window.

- univariate:

  Logical; TRUE for one target type, FALSE for two.

- mark1:

  Variable used to mark the points in the point pattern object for the
  first type.

- mark2:

  Variable used to mark the points in the point pattern object for the
  second type (optional, only used if `univariate` is FALSE).

- split_max:

  Maximum point count per initial grid cell. NULL uses `sqrt(n)`, where
  n denotes the total number of points in `ppp_obj`.

- n_min:

  Minimum target/background counts, or type1/type2/background counts.
  NULL uses the paper's square-root thresholds based on full-pattern
  counts.

- rho1:

  Phase I bound on longer/shorter side ratio. NULL uses the full-window
  ratio.

- rho2:

  Phase II aspect ratio bound. Defaults to Inf; must be \>= rho1.

- max_grid_dim:

  Maximum number of grid cells per axis to search. Defaults to 100.

## Value

A list whose elements are marked `ppp` blocks.

## Details

Adaptive spatial blocking

## Note

The observation window must be rectangular. Returned blocks retain their
original marks. Unassigned points are excluded with a warning.

## Examples

``` r
# one image from the ovarian cancer data in the package
ids <- unique(ovarian_df$sample_id)
img <- ovarian_df[ovarian_df$sample_id == ids[4], ]
pp <- spatstat.geom::ppp(img$x, img$y,
                         window = spatstat.geom::owin(range(img$x), range(img$y)),
                         marks = img$immune)

blocks <- bkamp_block(pp, mark1 = "immune")
length(blocks)
#> [1] 21
```
