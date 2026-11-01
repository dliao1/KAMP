#' Computes GAMP Expectation for Bivariate Point Patterns
#'
#' @title GAMP bivariate expectation
#'
#' @description
#' Computes the GAMP (G-function Adjusted for Marked Permutations) expectation
#' for bivariate point patterns: the exact mean of the uncorrected cross-type
#' nearest-neighbor G function (distance from each `mark1` point to its nearest
#' `mark2` point) under random relabeling of marks, without running any permutations.
#'
#' Observed G matches `spatstat.explore::Gcross(..., i = mark1, j = mark2, correction = "none")`.
#'
#' @param ppp_obj A point pattern object from the `spatstat.geom` package.
#' @param rvals Vector of radii at which to calculate the GAMP expectation. Defaults to c(0, 0.05, 0.075, 0.1, 0.15, 0.2).
#' @param correction Type of edge correction. Only "none" (no edge correction) is currently supported.
#' @param mark1 Variable used to mark the points in the point pattern object for the first type. Default is "immune1".
#' @param mark2 Variable used to mark the points in the point pattern object for the second type. Default is "immune2".
#'
#' @returns
#' A dataframe with the following columns:
#' \describe{
#'  \item{r}{The radius at which G was calculated.}
#'  \item{g}{The observed G value}
#'  \item{theo_csr}{The theoretical G under CSR, `1 - exp(-lambda2 * pi * r^2)`}
#'  \item{gamp_csr}{The adjusted CSR representing the GAMP permuted expectation.}
#'  \item{gamp}{The difference between observed G and GAMP CSR}
#' }
#'
#' @export
#'
#' @examples
#' win <- spatstat.geom::owin(c(0, 1), c(0, 1))
#' pp <- spatstat.random::rpoispp(lambda = 150, win = win)
#' mark_labels <- c("immune1", "immune2", "background")
#' marks <- sample(mark_labels, pp$n, replace = TRUE, prob = c(0.3, 0.3, 0.4))
#' marked_pp <- spatstat.geom::ppp(pp$x, pp$y, window = win, marks = factor(marks))
#'
#' result <- gamp_expectation_biv(marked_pp, mark1 = "immune1", mark2 = "immune2")
#' print(result)
gamp_expectation_biv <- function(ppp_obj,
                                 rvals = c(0, .05, .075, .1, .15, .2),
                                 correction = "none",
                                 mark1 = "immune1",
                                 mark2 = "immune2") {
  gamp_expectation_helper(ppp_obj = ppp_obj,
                          rvals = rvals,
                          correction = correction,
                          mark1 = mark1,
                          mark2 = mark2)
}
