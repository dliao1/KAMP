# kamp_variance_helper

Helper function to calculate the KAMP variance for a point pattern
object and single radius.

## Usage

``` r
kamp_variance_helper(ppp_obj, rvalue, correction = "trans", mark1 = "immune")
```

## Arguments

- ppp_obj:

  A point pattern object "ppp" from the `spatstat` package.

- rvalue:

  A single radius

- correction:

  Type of edge correction. Defaults to translational.

- mark1:

  Value used to mark the points in the point pattern object. Default is
  "immune".

## Value

A single-row dataframe with the following columns:

- r:

  The current radius at which K was calculated.

- k:

  The observed K value

- theo_csr:

  The theoretical K under CSR

- kamp_csr:

  The adjusted CSR representing the KAMP permuted expectation.

- kamp:

  The difference between observed K and KAMP CSR

- var:

  Variance of K under the permutation null distribution

- pval:

  P-value

## Details

Helper function for KAMP Variance
