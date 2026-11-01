# Close pairs and edge correction weights for the Rcpp KAMP variance

Finds all close pairs up to `max(rvals)`, builds their edge correction
weights, and sorts them by distance, as expected by the
`kamp_pair_sums_*` Rcpp functions.

## Usage

``` r
kamp_pair_prep(ppp_obj, rvals, correction)
```

## Arguments

- ppp_obj:

  A point pattern object of class "ppp".

- rvals:

  A vector of radii.

- correction:

  One of "trans", "iso", "border", or "none". Border pairs get weight 1;
  eligibility is handled in `kamp_pair_sums_bord`.

## Value

A list with sorted pair indices `i`, `j`, distances `d`, and weights
`w`.
