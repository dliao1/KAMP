#' Ovarian Cancer Example Point Pattern
#'
#' A processed spatial point pattern object derived from the VectraPolaris ovarian cancer dataset.
#' This object contains immune and background cells within the tumor region of one image.
#'
#' @format A `ppp` object (from the `spatstat.geom` package) with marks indicating immune status.
#' \describe{
#'   \item{x}{X}
#'   \item{y}{Y}
#'   \item{marks}{Cell type: "immune" or "background"}
#' }
#'
#' @source Subset of `HumanOvarianCancerVP()` from the VectraPolarisData package.
#'
"ovarian_df"
