#' Compute KAMP Expectation
#'
#' @title get_kamp_expectation
#'
#' @description
#' Computes the KAMP (K-function Adjusted for Marked Permutations) expectation
#' for a given spatial point pattern. This function calculates a Ripley's K
#' using both the traditional theoretical expectation (based on `Kcross`)
#' and the KAMP-adjusted CSR baseline (based on `Kest`)
#
#' @param ppp_obj A point pattern object of class from the `spatstat.geom` package.
#' @param rvec Vector of radii at which to calculate the KAMP expectation. Defaults to c(0, 0.05, 0.075, 0.1, 0.15, 0.2).
#' @param correction Specifies the edge correction method to be used and passed to `Kcross` and `Kest`.
#' @param markvar Identifies subset of marked points, defaults to immune
#' @param thin_pct Percentage that determines how much to thin the amount of points in the point pattern object.
#'
#' @return
#' A tibble with the following columns:
#' \describe{
#'   \item{r}{The radius at which the K-function was calculated.}
#'   \item{theo_csr}{The theoretical K under CSR from `Kcross`}
#'   \item{kamp_csr}{The adjusted CSR representing the permuted expectation.}
#'   \item{k}{The observed K-function value from `Kcross`}
#'   \item{kamp_fundiff}{The difference between observed K and KAMP CSR}
#' }
#'
#' @importFrom spatstat.explore Kcross Kest
#' @importFrom spatstat.geom area.owin ppp as.owin
#' @importFrom dplyr mutate select
#' @importFrom tibble as_tibble
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @examples
#' if (requireNamespace("spatstat.geom", quietly = TRUE)) {
#'   # simulate a simple spatial point pattern with two types
#'   win <- spatstat.geom::owin(c(0, 1), c(0, 1))
#'   pp <- spatstat.random::rpoispp(lambda = 100, win = win)
#'   marks <- sample(c("immune", "background"), pp$n, replace = TRUE)
#'   marked_pp <- spatstat.geom::ppp(pp$x, pp$y, window = win, marks = factor(marks))
#'
#'   # compute KAMP expectation
#'   kamp_result <- get_kamp_expectation(marked_pp, markvar = "immune")
#'   print(kamp_result)
#' }
get_kamp_expectation <- function(ppp_obj,
                                 rvec = c(0, .05, .075, .1, .15, .2),
                                 correction = "trans",
                                 markvar = "immune",
                                 thin_pct = 0) {

  # Pre-existing code that uses spatstat
  # Gets original K using translational correction
  k_orig = Kcross(ppp_obj, i = markvar, j = markvar,
             r = rvec,
             correction = correction)


  kamp = Kest(ppp_obj,
              r = rvec,
              correction = correction) %>%
    as_tibble() %>%
    mutate(kamp_csr = trans,
           k = k_orig$trans, # takes calculated K from Kcross
           theo_csr = k_orig$theo,
           kamp_fundiff = k - kamp_csr) %>% # takes calculated theoretical K from Kcross
    select(r, theo_csr, kamp_csr, k, kamp_fundiff) # gets r, kamp CSR, K statistic, method

  return(kamp)
}


#' Title
#' @title get_kamp_expectation_mat
#' @description Wrapper function for mapping across multiple radii to calculate the KAMP expectation using matrices.
#' @param ppp_obj  A point pattern object of class "ppp" from the spatstat package.
#' @param rvec A vector of radii at which to calculate the KAMP expectation.
#' @param markvar The variable used to mark the points in the point pattern object. Default is "immune".
#'
#' @returns A tibble containing the radius, expectation, and time taken to compute the expectation.
#' @export
#'
#' @examples
get_kamp_expectation_mat = function(ppp_obj,
                                    rvec,
                                    markvar = "immune") {
  return(map_dfr(rvec, get_kamp_expectation_mat_helper, ppp_obj = ppp_obj, .progress = TRUE))
}
