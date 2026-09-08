# Checks inputs for KAMP functions

Checks inputs for KAMP functions

## Usage

``` r
check_inputs(
  df,
  rvals,
  univariate,
  correction,
  mark_var,
  mark1,
  mark2,
  variance,
  thin,
  p_thin,
  ...
)
```

## Arguments

- df:

  A dataframe containing the point pattern data. Will be converted into
  a `ppp` object.

- rvals:

  Vector of radius values at which to compute the KAMP expectation.

- univariate:

  Logical indicating whether to compute univariate KAMP expectation.
  Defaults to TRUE.

- mark_var:

  Column name in `df` that contains the marks for the point pattern
  object.

- mark1:

  Value used to mark the points in the point pattern object.

- mark2:

  Value used to mark the points in the point pattern object for the
  second type (optional, only used if `univariate` is FALSE).

- variance:

  Logical indicating whether to compute the variance of KAMP (default is
  FALSE).

- thin:

  Logical indicating whether to thin the point pattern before computing
  KAMP (default is FALSE), called KAMP-lite.

- p_thin:

  Percentage that determines how much to thin

- ...:

  Additional arguments (currently unused).

## Value

The `ppp` point pattern object (built from `df` if it was a data.frame)
if all inputs are valid, otherwise throws an error with a descriptive
message.
