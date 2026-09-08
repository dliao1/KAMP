# Ovarian Cancer Example Dataframe

A processed dataframe derived from the VectraPolaris ovarian cancer
dataset. It is a snapshot of 5 images from the `HumanOvarianCancerVP()`
dataset in the `VectraPolarisData` package

## Usage

``` r
ovarian_df
```

## Format

A dataframe with x and y coordinates and marks indicating immune status.

- cell_id:

  The ID of the current cell in the current sample/image

- sample_id:

  The ID of the sample/image

- x:

  X coordinate of the current cell

- y:

  Y coordinate of the current cell

- immune:

  A factor indicating whether the cell is immune or background

- phenotype:

  A factor indicating the cell type (e.g., b cell, cytotoxic t cell,
  helper t cell, macrophage, tumor, other)

## Source

Subset of `HumanOvarianCancerVP()` from the VectraPolarisData package.
