#' GAMP Variance
#' @title GAMP univariate variance
#'
#' @description
#' Computes the GAMP (G-function Adjusted for Marked Permutations) variance
#' for a given spatial point pattern: the exact variance of the uncorrected
#' nearest-neighbor G function of the `mark1` points under random relabeling
#' of marks. Also returns the GAMP expectation and p-value.
#'
#' The pairwise sums are done in Rcpp (`gamp_sums`) from
#' `spatstat.geom::closepairs()`, so only pairs within `max(rvals)` are touched.
#'
#' @param ppp_obj A point pattern object of class "ppp" from the spatstat package.
#' @param rvals A vector of radii at which to calculate the GAMP variance.
#' @param correction Type of edge correction. Only "none" (no edge correction) is currently supported.
#' @param mark1 The variable used to mark the points in the point pattern object. Default is "immune".
#'
#' @returns
#' A dataframe with the following columns:
#' \describe{
#'   \item{r}{The radius at which G was calculated.}
#'   \item{g}{The observed G value}
#'   \item{theo_csr}{The theoretical G under CSR}
#'   \item{gamp_csr}{The adjusted CSR representing the GAMP permuted expectation.}
#'   \item{gamp}{The difference between observed G and GAMP CSR}
#'   \item{var}{Variance of G under the permutation null distribution}
#'   \item{pvalue}{P-value, calculated using the formula: pnorm(-z). Set to 1 when the variance is 0.}
#' }
#'
#' @export
#' @useDynLib KAMP, .registration = TRUE
#' @examples
#' win <- spatstat.geom::owin(c(0, 1), c(0, 1))
#' pp <- spatstat.random::rpoispp(lambda = 150, win = win)
#' marks <- sample(c("immune", "background"), pp$n, replace = TRUE, prob = c(0.4, 0.6))
#' marked_pp <- spatstat.geom::ppp(pp$x, pp$y, window = win, marks = factor(marks))
#'
#' result <- gamp_variance(marked_pp, rvals = c(0.05, 0.1), mark1 = "immune")
#' print(result)
gamp_variance <- function(ppp_obj,
                          rvals = c(0, .05, .075, .1, .15, .2),
                          correction = "none",
                          mark1 = "immune") {
  gamp_variance_helper(ppp_obj = ppp_obj,
                       rvals = rvals,
                       correction = correction,
                       mark1 = mark1)
}
