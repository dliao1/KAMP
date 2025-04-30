#' Computes KAMP Variance for Bivariate Point Patterns
#'
#' @title Computes bivariate KAMP variance
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
#' @param rvec A vector of radii at which to calculate the KAMP expectation. Defaults to c(0, 0.05, 0.075, 0.1, 0.15, 0.2).
#' @param correction Type of edge correction. Defaults to translational.
#' @param markvar1 Variable used to mark the points in the point pattern object for the first type. Default is "immune1".
#' @param markvar2 Variable used to mark the points in the point pattern object for the second type. Default is "immune2".
#' @param thin_pct Percentage that determines how much to thin the amount of points in the point pattern object. Default is 0.
#'
#' @importFrom spatstat.explore Kcross Kest edge.Trans edge.Ripley
#' @importFrom spatstat.geom area.owin ppp as.owin npoints Window
#' @importFrom dplyr mutate select
#' @importFrom tibble as_tibble
#' @importFrom magrittr %>%
#' @importFrom purrr map_dfr
#' @importFrom tictoc tic toc
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
#'
#' @examples
#'
#' if (requireNamespace("spatstat.geom", quietly = TRUE) &&
#'     requireNamespace("spatstat.random", quietly = TRUE)) {
#'
#'   # Simulates a window
#'   win <- spatstat.geom::owin(c(0, 1), c(0, 1))
#'
#'   # Generates 200 total points
#'   pp <- spatstat.random::rpoispp(lambda = 200, win = win)
#'
#'   # Assigns three marks: immune1, immune2, and background
#'   marks <- sample(c("immune1", "immune2", "background"), pp$n, replace = TRUE, prob = c(0.3, 0.3, 0.4))
#'
#'   # Creates marked point pattern
#'   marked_pp <- spatstat.geom::ppp(pp$x, pp$y, window = win, marks = factor(marks))
#'
#'   # Computes KAMP variance
#'   result <- kamp_variance_biv(marked_pp, markvar1 = "immune1", markvar2 = "immune2")
#'   print(result)
#' }
kamp_variance_biv <- function(ppp_obj,
                                 rvec = c(0, .05, .075, .1, .15, .2),
                                 correction = "trans",
                                 markvar1 = "immune1",
                                 markvar2 = "immune2",
                                 thin_pct = 0) {

  check_valid_inputs_biv(ppp_obj = ppp_obj,
                        rvec = rvec,
                        correction = correction,
                        markvar1 = markvar1,
                        markvar2 = markvar2,
                        thin_pct = thin_pct)


  if (thin_pct != 0) {
    ppp_obj = rthin(ppp_obj, 1 - thin_pct)
  }

  map_dfr(rvec,
          ~kamp_variance_biv_helper(ppp_obj = ppp_obj,
                                    rvalue = .x,
                                    correction = correction,
                                    markvar1 = markvar1,
                                    markvar2 = markvar2),
          .progress = TRUE)
}
