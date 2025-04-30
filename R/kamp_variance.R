#' KAMP Variance
#' @title Computes univariate KAMP variance
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
#' @param rvec A vector of radii at which to calculate the KAMP variance.
#' @param correction Type of edge correction. Defaults to translational.
#' @param markvar The variable used to mark the points in the point pattern object. Default is "immune".
#' @param thin_pct Percentage that determines how much to thin the amount of points in the point pattern object.
#'
#' @returns
#' A dataframe with the following columns:
#' \describe{
#'   \item{r}{The radius at which K was calculated.}
#'   \item{k}{The observed K value}
#'   \item{theo_csr}{The theoretical K under CSR}
#'   \item{kamp_csr}{The adjusted CSR representing the KAMP permuted expectation.}
#'   \item{kamp_fundiff}{The difference between observed K and KAMP CSR}
#'   \item{var}{Variance of K under the permutation null distribution}
#'   \item{z}{Z statistic, calculated by normalizing K using the formula: (K - KAMP)/sqrt(var)}
#'   \item{pval}{P-value, calculated using the formula: pnorm(-z)}
#' }
#'
#' @export
#'
#' @examples
#' if (requireNamespace("spatstat.geom", quietly = TRUE)) {
#'   # simulates a simple spatial point pattern with two types
#'   win <- spatstat.geom::owin(c(0, 1), c(0, 1))
#'   pp <- spatstat.random::rpoispp(lambda = 100, win = win)
#'   marks <- sample(c("immune", "background"), pp$n, replace = TRUE)
#'   marked_pp <- spatstat.geom::ppp(pp$x, pp$y, window = win, marks = factor(marks))
#'
#'   # computes KAMP variance
#'   kamp_result <- kamp_variance(marked_pp, markvar = "immune")
#'   print(kamp_result)
#' }
kamp_variance = function(ppp_obj,
                         rvec = c(0, .05, .075, .1, .15, .2),
                         correction = "trans",
                         markvar = "immune",
                         thin_pct = 0) {

  check_valid_inputs_univ(ppp_obj = ppp_obj,
                          rvec = rvec,
                          correction = correction,
                          markvar = markvar,
                          thin_pct = thin_pct)


  if (thin_pct != 0) {
    ppp_obj = rthin(ppp_obj, 1 - thin_pct)
  }

  map_dfr(rvec,
          ~kamp_variance_helper(ppp_obj = ppp_obj,
                                rvalue = .x,
                                correction = correction,
                                markvar = markvar),
          .progress = TRUE)
}
