#' Computes KAMP Variance for Bivariate Point Patterns
#'
#' @title KAMP bivariate variance
#' @description
#' Computes the KAMP (K-function Adjusted for Marked Permutations) variance
#' for bivariate point patterns. This function calculates Ripley's K
#' using both the traditional Ripley's K method and the KAMP-adjusted CSR
#' baseline by using a matrix-based implementation.
#'
#' The KAMP-adjusted CSR represents a more realistic baseline for K (compared
#' to traditional CSR) that accounts for spatial clustering or inhomogeneity
#' in a point pattern compared to the traditional CSR assumption, while
#' avoiding the computational burden of permuting the point pattern.
#'
#' Note: This function implements a slower, matrix-based implementation
#' of the KAMP variance. It is a wrapper around the `kamp_variance_biv_helper`
#' function that calculates the KAMP variance at one radius and maps it
#' over a vector of radii.
#'
#' @param ppp_obj A point pattern object from the `spatstat.geom` package.
#' @param rvals A vector of radii at which to calculate the KAMP expectation. Defaults to c(0, 0.05, 0.075, 0.1, 0.15, 0.2).
#' @param correction Type of edge correction. Defaults to translational.
#' @param mark1 Variable used to mark the points in the point pattern object for the first type. Default is "immune1".
#' @param mark2 Variable used to mark the points in the point pattern object for the second type. Default is "immune2".
#'
#' @importFrom spatstat.explore Kcross Kest edge.Trans edge.Ripley
#' @importFrom spatstat.geom area.owin ppp as.owin npoints Window
#' @importFrom dplyr mutate select
#' @importFrom tibble as_tibble
#' @importFrom magrittr %>%
#' @importFrom purrr map_dfr
#' @importFrom stats dist pnorm
#' @importFrom tibble tibble
#'
#' @returns
#' A dataframe with the following columns:
#' \describe{
#' \item{r}{The radius at which K was calculated.}
#' \item{k}{The observed K value}
#' \item{theo_csr}{The theoretical K under CSR}
#' \item{kamp_csr}{The adjusted CSR representing the KAMP permuted expectation.}
#' \item{var}{Variance of K under the permutation null distribution}
#' \item{z}{Z statistic, calculated by normalizing K using the formula: (K - KAMP)/sqrt(var)}
#' \item{pval}{P-value, calculated using the formula: pnorm(-z)}
#' }
#'
#' @export
#' @examples
#' # one image from the ovarian cancer data in the package
#' ids <- unique(ovarian_df$sample_id)
#' img <- ovarian_df[ovarian_df$sample_id == ids[4], ]
#' pp <- spatstat.geom::ppp(img$x, img$y,
#'                          window = spatstat.geom::owin(range(img$x), range(img$y)),
#'                          marks = img$phenotype)
#'
#' kamp_variance_biv(pp, rvals = c(25, 50),
#'                   mark1 = "cytotoxic t cell", mark2 = "helper t cell")
kamp_variance_biv <- function(ppp_obj,
                                 rvals = c(0, .05, .075, .1, .15, .2),
                                 correction = "trans",
                                 mark1 = "immune1",
                                 mark2 = "immune2") {

  map_dfr(rvals,
          ~kamp_variance_biv_helper(ppp_obj = ppp_obj,
                                    rval = .x,
                                    correction = correction,
                                    mark1 = mark1,
                                    mark2 = mark2),
          .progress = TRUE)
}
