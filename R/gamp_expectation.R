#' Compute GAMP Expectation
#'
#' @title GAMP univariate expectation
#'
#' @description
#' Computes the GAMP (G-function Adjusted for Marked Permutations) expectation
#' for a given spatial point pattern: the exact mean of the uncorrected
#' nearest-neighbor G function of the `mark1` points under random relabeling
#' of marks, without running any permutations.
#'
#' Observed G is the uncorrected nearest-neighbor distance distribution of the
#' `mark1` points, matching `spatstat.explore::Gest(..., correction = "none")`.
#'
#' @param ppp_obj A point pattern object from the `spatstat.geom` package.
#' @param rvals Vector of radii at which to calculate the GAMP expectation. Defaults to c(0, 0.05, 0.075, 0.1, 0.15, 0.2).
#' @param correction Type of edge correction. Only "none" (no edge correction) is currently supported.
#' @param mark1 Identifies subset of marked points. Defaults to immune.
#'
#' @returns
#' A dataframe with the following columns:
#' \describe{
#'   \item{r}{The radius at which G was calculated.}
#'   \item{g}{The observed G value}
#'   \item{theo_csr}{The theoretical G under CSR, `1 - exp(-lambda * pi * r^2)`}
#'   \item{gamp_csr}{The adjusted CSR representing the GAMP permuted expectation.}
#'   \item{gamp}{The difference between observed G and GAMP CSR}
#' }
#'
#' @export
#'
#' @examples
#' win <- spatstat.geom::owin(c(0, 1), c(0, 1))
#' pp <- spatstat.random::rpoispp(lambda = 150, win = win)
#' marks <- sample(c("immune", "background"), pp$n, replace = TRUE, prob = c(0.4, 0.6))
#' marked_pp <- spatstat.geom::ppp(pp$x, pp$y, window = win, marks = factor(marks))
#'
#' result <- gamp_expectation(marked_pp, rvals = c(0, 0.05, 0.1), mark1 = "immune")
#' print(result)
gamp_expectation <- function(ppp_obj,
                             rvals = c(0, .05, .075, .1, .15, .2),
                             correction = "none",
                             mark1 = "immune") {
  gamp_expectation_helper(ppp_obj = ppp_obj,
                          rvals = rvals,
                          correction = correction,
                          mark1 = mark1)
}
