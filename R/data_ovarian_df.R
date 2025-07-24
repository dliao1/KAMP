#' Ovarian Cancer Example Dataframe
#'
#' A processed dataframe derived from the VectraPolaris ovarian cancer dataset.
#' It is a snapshot of 5 images from the `HumanOvarianCancerVP()` dataset in the `VectraPolarisData` package
#'
#' @format A dataframe with x and y coordinates and marks indicating immune status.
#' \describe{
#'   \item{cell_id}{The ID of the current cell in the current sample/image}
#'   \item{sample_id}{The ID of the sample/image}
#'   \item{x}{X coordinate of the current cell}
#'   \item{y}{Y coordinate of the current cell}
#'   \item{immune}{A factor indicating whether the cell is immune or background}
#'   \item{phenotype}{A factor indicating the cell type (e.g., b cell, cytotoxic t cell, helper t cell, macrophage, tumor, other)}
#' }
#'
#' @source Subset of `HumanOvarianCancerVP()` from the VectraPolarisData package.
#'
"ovarian_df"
