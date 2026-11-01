#' Computes GAMP Variance for Bivariate Point Patterns
#'
#' @title GAMP bivariate variance
#' @description
#' Computes the GAMP (G-function Adjusted for Marked Permutations) variance
#' for bivariate point patterns: the exact variance of the uncorrected
#' cross-type nearest-neighbor G function (distance from each `mark1` point to
#' its nearest `mark2` point) under random relabeling of marks. Also returns
#' the GAMP expectation and p-value.
#'
#' The pairwise sums are done in Rcpp (`gamp_sums`) from
#' `spatstat.geom::closepairs()`, so only pairs within `max(rvals)` are touched.
#'
#' @param ppp_obj A point pattern object from the `spatstat.geom` package.
#' @param rvals A vector of radii at which to calculate the GAMP variance. Defaults to c(0, 0.05, 0.075, 0.1, 0.15, 0.2).
#' @param correction Type of edge correction. Only "none" (no edge correction) is currently supported.
#' @param mark1 Variable used to mark the points in the point pattern object for the first type. Default is "immune1".
#' @param mark2 Variable used to mark the points in the point pattern object for the second type. Default is "immune2".
#'
#' @returns
#' A dataframe with the following columns:
#' \describe{
#' \item{r}{The radius at which G was calculated.}
#' \item{g}{The observed G value}
#' \item{theo_csr}{The theoretical G under CSR}
#' \item{gamp_csr}{The adjusted CSR representing the GAMP permuted expectation.}
#' \item{gamp}{The difference between observed G and GAMP CSR}
#' \item{var}{Variance of G under the permutation null distribution}
#' \item{pvalue}{P-value, calculated using the formula: pnorm(-z). Set to 1 when the variance is 0.}
#' }
#'
#' @export
#' @useDynLib KAMP, .registration = TRUE
#' @examples
#' # one image from the ovarian cancer data in the package
#' ids <- unique(ovarian_df$sample_id)
#' img <- ovarian_df[ovarian_df$sample_id == ids[4], ]
#' pp <- spatstat.geom::ppp(img$x, img$y,
#'                          window = spatstat.geom::owin(range(img$x), range(img$y)),
#'                          marks = img$phenotype)
#'
#' gamp_variance_biv(pp, rvals = c(10, 20, 30, 40),
#'                   mark1 = "cytotoxic t cell", mark2 = "helper t cell")
gamp_variance_biv <- function(ppp_obj,
                              rvals = c(0, .05, .075, .1, .15, .2),
                              correction = "none",
                              mark1 = "immune1",
                              mark2 = "immune2") {
  gamp_variance_helper(ppp_obj = ppp_obj,
                       rvals = rvals,
                       correction = correction,
                       mark1 = mark1,
                       mark2 = mark2)
}
