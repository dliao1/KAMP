#' KAMP Variance
#' @title KAMP univariate variance
#'
#' @description
#' Computes the KAMP (K-function Adjusted for Marked Permutations) variance
#' for a given spatial point pattern. Also returns the KAMP expectation,
#' z-statistic, and p-value.
#'
#' Note: this a matrix-based implementation of the KAMP variance that
#' does not use the `spatstat` package. It is a wrapper around the
#' `kamp_variance_helper` function that calculates the KAMP variance at
#' one radius and maps it over a vector of radii.
#'
#' @param ppp_obj A point pattern object of class "ppp" from the spatstat package.
#' @param rvals A vector of radii at which to calculate the KAMP variance.
#' @param correction Type of edge correction. Defaults to translational.
#' @param mark1 The variable used to mark the points in the point pattern object. Default is "immune".
#'
#' @returns
#' A dataframe with the following columns:
#' \describe{
#'   \item{r}{The radius at which K was calculated.}
#'   \item{k}{The observed K value}
#'   \item{theo_csr}{The theoretical K under CSR}
#'   \item{kamp_csr}{The adjusted CSR representing the KAMP permuted expectation.}
#'   \item{kamp}{The difference between observed K and KAMP CSR}
#'   \item{var}{Variance of K under the permutation null distribution}
#'   \item{pval}{P-value, calculated using the formula: pnorm(-z)}
#' }
#'
#' @export
#' @importFrom purrr map_dfr
#' @importFrom dplyr %>%
#' @importFrom tibble as_tibble
kamp_variance = function(ppp_obj,
                         rvals = c(0, .05, .075, .1, .15, .2),
                         correction = "trans",
                         mark1 = "immune") {


  results <- map_dfr(rvals,
          ~kamp_variance_helper(ppp_obj = ppp_obj,
                                rvalue = .x,
                                correction = correction,
                                mark1 = mark1),
          .progress = TRUE)

  return(results)
}
