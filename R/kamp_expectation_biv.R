#' Computes KAMP Expectation for Bivariate Point Patterns
#'
#' @title KAMP bivariate Expectation
#'
#' @description Computes the KAMP (K-function Adjusted for Marked Permutations) expectation
#' for bivariate point patterns. This function calculates Ripley's K
#' using both the traditional Ripley's K method (based on `Kcross`)
#' and the KAMP-adjusted CSR baseline (based on `Kest`).
#'
#' The KAMP-adjusted CSR represents a more robust baseline for K (compared
#' to traditional CSR) that accounts for spatial clustering or inhomogeneity
#' in a point pattern compared to the traditional CSR assumption, while
#' avoiding the computational burden of permuting the point pattern.
#'
#' Note: This function uses the `spatstat` package under the hood.
#' See `?Kcross` and `?Kest` for more details on the K calculation methods.
#'
#' See `kamp_expectation_biv_mat` for the matrix-based implementation of the KAMP
#' bivariate expectation.
#'
#' @param ppp_obj A point pattern object from the `spatstat.geom` package.
#' @param rvec Vector of radii at which to calculate the KAMP expectation. Defaults to c(0, 0.05, 0.075, 0.1, 0.15, 0.2).
#' @param correction Type of edge correction. Defaults to translational.
#' @param markvar1 Variable used to mark the points in the point pattern object for the first type. Default is "immune1".
#' @param markvar2 Variable used to mark the points in the point pattern object for the second type. Default is "immune2".
#' @param thin_pct Percentage that determines how much to thin the amount of points in the point pattern object. Default is 0.
#'
#' @returns
#' A dataframe with the following columns:
#' \describe{
#'  \item{r}{The radius at which K was calculated.}
#'  \item{k}{The observed K value}
#'  \item{theo_csr}{The theoretical K under CSR}
#'  \item{kamp_csr}{The adjusted CSR representing the KAMP permuted expectation.}
#'  \item{kamp_fundiff}{The difference between observed K and KAMP CSR}
#' }
#'
#' @importFrom spatstat.explore Kcross Kest
#' @importFrom spatstat.geom area.owin ppp as.owin
#' @importFrom spatstat.random rthin
#' @importFrom dplyr mutate select
#' @importFrom tibble as_tibble
#' @importFrom magrittr %>%
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
#'   # Computes KAMP expectation
#'   result <- kamp_expectation_biv(marked_pp, markvar1 = "immune1", markvar2 = "immune2")
#'   print(result)
#' }
kamp_expectation_biv <- function(ppp_obj,
                             rvec = c(0, .05, .075, .1, .15, .2),
                             correction = "trans",
                             markvar1 = "immune1",
                             markvar2 = "immune2",
                             thin_pct = 0) {

  if (!is.numeric(thin_pct)) {
    stop("thin_pct must be numeric.")
  }

  if (thin_pct < 0 || thin_pct > 1) {
    stop("thin_pct must be between 0 and 1.")
  }

  if (markvar1 == markvar2) {
    stop("markvar1 and markvar2 must be different.")
  }

  if (thin_pct != 0) {
    ppp_obj = rthin(ppp_obj, 1 - thin_pct)
  }

  k_orig = Kcross(ppp_obj, i = markvar1, j = markvar2,
             r = rvec,
             correction = correction)

  kamp = Kest(ppp_obj,
              r = rvec,
              correction = correction) %>%
    as_tibble() %>%
    mutate(kamp_csr = trans,
           k = k_orig$trans,
           theo_csr = k_orig$theo,
           kamp_fundiff = k - kamp_csr) %>%
    select(r, k, theo_csr, kamp_csr, kamp_fundiff)

  return(kamp)
}

#' Bivariate KAMP Expectation (Matrix Implementation)
#'
#' @title KAMP bivariate expectation (Matrix Implementation)
#' @description
#' Computes the KAMP (K-function Adjusted for Marked Permutations) expectation
#' for bivariate point patterns using a matrix-based approach. Note that this
#' is slower.
#'
#' See `kamp_expectation_biv` for the `spatstat`-based implementation of the KAMP
#' bivariate expectation.
#'
#' @param ppp_obj A point pattern object from the `spatstat.geom` package.
#' @param rvec Vector of radii at which to calculate the KAMP expectation. Defaults to c(0, 0.05, 0.075, 0.1, 0.15, 0.2).
#' @param correction Type of edge correction method. Defaults to translational.
#' @param markvar1 Variable used to mark the points in the point pattern object for the first type. Default is "immune1".
#' @param markvar2 Variable used to mark the points in the point pattern object for the second type. Default is "immune2".
#' @param thin_pct Percentage that determines how much to thin the amount of points in the point pattern object. Default is 0.
#'
#' @returns
#' A dataframe with the following columns:
#' \describe{
#'  \item{r}{The radius at which K was calculated.}
#'  \item{k}{The observed K value}
#'  \item{theo_csr}{The theoretical K under CSR}
#'  \item{kamp_csr}{The adjusted CSR representing the KAMP permuted expectation.}
#'  \item{kamp_fundiff}{The difference between observed K and KAMP CSR}
#' }
#'
#'
#' @importFrom spatstat.explore Kcross Kest
#' @importFrom spatstat.geom area.owin ppp as.owin
#' @importFrom spatstat.random rthin
#' @importFrom dplyr mutate select
#' @importFrom tibble as_tibble
#' @importFrom magrittr %>%
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
#'   # Computes KAMP expectation
#'   result <- kamp_expectation_biv_mat(marked_pp, markvar1 = "immune1", markvar2 = "immune2")
#'   print(result)
#' }
kamp_expectation_biv_mat <- function(ppp_obj,
                                     rvec = c(0, .05, .075, .1, .15, .2),
                                     correction = "trans",
                                     markvar1 = "immune1",
                                     markvar2 = "immune2",
                                     thin_pct = 0) {

  if (!is.numeric(thin_pct)) {
    stop("thin_pct must be numeric.")
  }

  if (thin_pct < 0 || thin_pct > 1) {
    stop("thin_pct must be between 0 and 1.")
  }

  if (markvar1 == markvar2) {
    stop("markvar1 and markvar2 must be different.")
  }

  if (thin_pct != 0) {
    ppp_obj = rthin(ppp_obj, 1 - thin_pct)
  }

  map_dfr(rvec,
          ~kamp_expectation_biv_mat_helper(ppp_obj = ppp_obj,
                                           rvalue = .x,
                                           correction = correction,
                                           markvar1 = markvar1,
                                           markvar2 = markvar2),
          .progress = TRUE)
}
