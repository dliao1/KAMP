#' Compute KAMP Expectation
#'
#' @title KAMP univariate expectation
#'
#' @description
#' Computes the KAMP (K-function Adjusted for Marked Permutations) expectation
#' for a given spatial point pattern. This function calculates Ripley's K
#' using both the traditional Ripley's K method (based on `Kcross`)
#' and the KAMP-adjusted CSR baseline (based on `Kest`).
#'
#' The KAMP-adjusted CSR represents a more robust baseline for K (compared
#' to traditional CSR) that accounts for spatial clustering or inhomogeneity
#' in a point pattern compared to the traditional CSR assumption, while
#' avoiding the computational burden of permuting the point pattern.
#'
#' Notes:
#'
#' This function uses the `spatstat` package under the hood, which
#' automatically uses border correction when the number of
#' points in the point pattern is more than 3000.
#'
#' See `?Kcross` and `?Kest` for more details on the K calculation methods.
#'
#' See `kamp_expectation_mat` for the matrix-based implementation.
#'
#'
#'
#
#' @param ppp_obj A point pattern object from the `spatstat.geom` package.
#' @param rvals Vector of radii at which to calculate the KAMP expectation. Defaults to c(0, 0.05, 0.075, 0.1, 0.15, 0.2).
#' @param correction Type of edge correction method to be used and passed to `Kcross` and `Kest`. Defaults to translational edge correction.
#' @param mark1 Identifies subset of marked points. Defaults to immune.
#'
#' @returns
#' A dataframe with the following columns:
#' \describe{
#'   \item{r}{The radius at which K was calculated.}
#'   \item{k}{The observed K value from `Kcross`}
#'   \item{theo_csr}{The theoretical K under CSR from `Kcross`}
#'   \item{kamp_csr}{The adjusted CSR representing the permuted expectation.}
#'   \item{kamp}{The difference between observed K and KAMP CSR}
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
kamp_expectation <- function(ppp_obj,
                             rvals = c(0, .05, .075, .1, .15, .2),
                             correction = "trans",
                             mark1 = "immune") {
  # gets number of points
  npts  <- npoints(ppp_obj)

  if (npts > 10000) {
    warning("Point pattern has more than 10,000 points. Switching to border correction")
    correction = "border"
  }
  # Pre-existing code that uses spatstat
  # Gets original K using translational correction
  k_orig = Kcross(ppp_obj, i = mark1, j = mark1,
             r = rvals,
             correction = correction)

  # Gets KAMP CSR using Kest
  kamp_df = Kest(ppp_obj,
              r = rvals,
              correction = correction) %>%
    as_tibble()

  if (correction == "trans") {
    kamp_df = kamp_df %>%
      mutate(kamp_csr = trans,
             k = k_orig$trans, # takes calculated K from Kcross
             theo_csr = k_orig$theo,
             kamp = k - kamp_csr) %>% # difference between K and KAMP CSR
      select(r, k, theo_csr, kamp_csr, kamp)
  } else if (correction == "iso") {
    kamp_df = kamp_df %>%
      mutate(kamp_csr = iso,
             k = k_orig$iso, # takes calculated K from Kcross
             theo_csr = k_orig$theo,
             kamp = k - kamp_csr) %>%
      select(r, k, theo_csr, kamp_csr, kamp)
  } else if (correction == "border") {
    kamp_df = kamp_df %>%
      mutate(kamp_csr = border,
             k = k_orig$border, # takes calculated K from Kcross
             theo_csr = k_orig$theo,
             kamp = k - kamp_csr) %>%
      select(r, k, theo_csr, kamp_csr, kamp)
  }

  return(kamp_df)
}

