# kamp_variance_biv_helper

Helper function to calculate the KAMP variance for bivariate point
patterns.

## Usage

``` r
kamp_variance_biv_helper(
  ppp_obj,
  rval,
  correction = "trans",
  mark1 = "immune1",
  mark2 = "immune2"
)
```

## Arguments

- ppp_obj:

  A point pattern object "ppp" from the `spatstat` package.

- rval:

  A single radius

- correction:

  Type of edge correction. Defaults to translational.

- mark1:

  Variable used to mark the points in the point pattern object for the
  first type. Default is "immune1".

- mark2:

  Variable used to mark the points in the point pattern object for the
  second type. Default is "immune2".

## Value

A single-row dataframe with the following columns:

- r:

  The radius at which K was calculated.

- k:

  The observed K value

- theo_csr:

  The theoretical K under CSR

- kamp_csr:

  The adjusted CSR representing the KAMP permuted expectation.

- var:

  Variance of K under the permutation null distribution

- pval:

  P-value, calculated using the formula: pnorm(-z)

## Details

Helper function for bivariate KAMP Variance
