# B-KAMP - Rcpp version

Identify adaptive blocks and combine block-level KAMP statistics
calculated using
[`kamp_variance_Rcpp()`](https://dliao1.github.io/KAMP/reference/kamp_variance_Rcpp.md)
or
[`kamp_variance_biv_Rcpp()`](https://dliao1.github.io/KAMP/reference/kamp_variance_biv_Rcpp.md).

## Usage

``` r
bkamp_Rcpp(
  ppp_obj,
  rvals = c(0, 0.05, 0.075, 0.1, 0.15, 0.2),
  univariate = TRUE,
  correction = "trans",
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

- rvals:

  Relative radii. Each value is multiplied by each block's shorter side
  to obtain the absolute radius used for within-block K calculations.

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

A tibble with the following columns:

- r:

  The relative radius at which K was calculated within each block.

- z:

  Weighted statistic combining block-level KAMP z-statistics.

- pvalue:

  right-sided normal p-value, `pnorm(-z)`.

## Details

B-KAMP for scalable clustering inference using Rcpp

## Note

Blocks are identified once and reused at all radii. Undefined block
statistics give combined NA results; radius 0 is typically undefined.

## See also

[`bkamp_block()`](https://dliao1.github.io/KAMP/reference/bkamp_block.md),
[`bkamp()`](https://dliao1.github.io/KAMP/reference/bkamp.md)

## Examples

``` r
# one image from the ovarian cancer data in the package
ids <- unique(ovarian_df$sample_id)
img <- ovarian_df[ovarian_df$sample_id == ids[4], ]
pp <- spatstat.geom::ppp(img$x, img$y,
                         window = spatstat.geom::owin(range(img$x), range(img$y)),
                         marks = img$immune)

bkamp_Rcpp(pp, rvals = c(0.05, 0.1, 0.15), mark1 = "immune")
#> # A tibble: 3 × 3
#>       r     z    pvalue
#>   <dbl> <dbl>     <dbl>
#> 1  0.05  27.7 1.74e-169
#> 2  0.1   36.5 3.25e-291
#> 3  0.15  37.5 2.25e-308
```
