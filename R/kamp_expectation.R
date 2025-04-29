#' Compute KAMP Expectation
#'
#' @title kamp_expectation
#'
#' @description
#' Computes the KAMP (K-function Adjusted for Marked Permutations) expectation
#' for a given spatial point pattern. This function calculates Ripley's K
#' using both the traditional Ripley's K method (based on `Kcross`)
#' and the KAMP-adjusted CSR baseline (based on `Kest`).
#'
#' The KAMP-adjusted CSR represents a more realistic baseline for K (compared
#' to traditional CSR) that accounts for spatial clustering or inhomogeneity
#' in a point pattern compared to the traditional CSR assumption, while
#' avoiding the computational burden of permuting the point pattern.
#'
#' See `?Kcross` and `?Kest` for more details on the K calculation methods.
#' See `get_kamp_expectation_mat` for the matrix-based implementation.
#'
#
#' @param ppp_obj A point pattern object from the `spatstat.geom` package.
#' @param rvec Vector of radii at which to calculate the KAMP expectation. Defaults to c(0, 0.05, 0.075, 0.1, 0.15, 0.2).
#' @param correction Specifies the edge correction method to be used and passed to `Kcross` and `Kest`.
#' @param markvar Identifies subset of marked points, defaults to immune
#' @param thin_pct Percentage that determines how much to thin the amount of points in the point pattern object.
#'
#' @return
#' A dataframe with the following columns:
#' \describe{
#'   \item{r}{The radius at which K was calculated.}
#'   \item{k}{The observed K value from `Kcross`}
#'   \item{theo_csr}{The theoretical K under CSR from `Kcross`}
#'   \item{kamp_csr}{The adjusted CSR representing the permuted expectation.}
#'   \item{kamp_fundiff}{The difference between observed K and KAMP CSR}
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
#' if (requireNamespace("spatstat.geom", quietly = TRUE)) {
#'   # simulates a simple spatial point pattern with two types
#'   win <- spatstat.geom::owin(c(0, 1), c(0, 1))
#'   pp <- spatstat.random::rpoispp(lambda = 100, win = win)
#'   marks <- sample(c("immune", "background"), pp$n, replace = TRUE)
#'   marked_pp <- spatstat.geom::ppp(pp$x, pp$y, window = win, marks = factor(marks))
#'
#'   # computes KAMP expectation
#'   kamp_result <- kamp_expectation(marked_pp, markvar = "immune")
#'   print(kamp_result)
#' }
kamp_expectation <- function(ppp_obj,
                             rvec = c(0, .05, .075, .1, .15, .2),
                             correction = "trans",
                             markvar = "immune",
                             thin_pct = 0) {

  if (thin_pct != 0) {
    # Thin the point pattern object
    #npts = npoints(ppp_obj)
    #npts_thin = round(npts * (1 - thin_pct))
    #ppp_obj = ppp_obj[sample(1:npts, npts_thin), ]
    ppp_obj = rthin(ppp_obj, 1 - thin_pct)
  }

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
    select(r, k, theo_csr, kamp_csr, kamp_fundiff) # gets r, kamp CSR, K statistic, method

  return(kamp)
}


#' Computes KAMP Expectation using matrices
#' @title kamp_expectation_mat
#'
#' @description
#' Computes the KAMP (K-function Adjusted for Marked Permutations) expectation
#' for a given spatial point pattern. This function calculates Ripley's K
#' using both the traditional Ripley's K method (based on `Kcross`)
#' and the KAMP-adjusted CSR baseline (based on `Kest`).
#'
#' The KAMP-adjusted CSR represents a more realistic baseline for K (compared
#' to traditional CSR) that accounts for spatial clustering or inhomogeneity
#' in a point pattern compared to the traditional CSR assumption, while
#' avoiding the computational burden of permuting the point pattern.
#'
#' Note: this is a slower, matrix-based implementation of the KAMP expectation
#'
#' See `get_kamp_expectation` for the faster implementation that utilizes
#' functions in `spatstat`
#'
#' @param ppp_obj  A point pattern object of class "ppp" from the spatstat package.
#' @param rvec A vector of radii at which to calculate the KAMP expectation.
#' @param correction Type of border correction (can either be translational or border)
#' @param markvar The variable used to mark the points in the point pattern object. Default is "immune".
#'
#' @returns
#' A dataframe with the following columns:
#' \describe{
#'   \item{r}{The radius at which K was calculated.}
#'   \item{k}{The observed K value calculated using matrices}
#'   \item{theo_csr}{The theoretical K under CSR}
#'   \item{kamp_csr}{The adjusted CSR representing the permuted expectation.}
#'   \item{kamp_fundiff}{The difference between observed K and KAMP CSR}
#' }
#'
#' @export
#'
#' @examples
#' if (requireNamespace("spatstat.geom", quietly = TRUE)) {
#' # simulate a simple spatial point pattern with two types
#' win <- spatstat.geom::owin(c(0, 1), c(0, 1))
#' pp <- spatstat.random::rpoispp(lambda = 100, win = win)
#' marks <- sample(c("immune", "background"), pp$n, replace = TRUE)
#' marked_pp <- spatstat.geom::ppp(pp$x, pp$y, window = win, marks = factor(marks))
#' # compute KAMP expectation using matrix method
#' kamp_result_mat <- kamp_expectation_mat(marked_pp, markvar = "immune")
#' print(kamp_result_mat)
#' }
kamp_expectation_mat = function(ppp_obj,
                                rvec = c(0, .05, .075, .1, .15, .2),
                                correction = "trans",
                                markvar = "immune") {
  map_dfr(rvec,
          ~kamp_expectation_mat_helper(ppp_obj = ppp_obj,
                                       rvalue = .x,
                                       correction = correction,
                                       markvar = markvar),
          .progress = TRUE)
}

