#' Computes KAMP Expectation for Bivariate Point Patterns
#'
#' @title KAMP bivariate expectation
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
#' @param ppp_obj A point pattern object from the `spatstat.geom` package.
#' @param rvals Vector of radii at which to calculate the KAMP expectation. Defaults to c(0, 0.05, 0.075, 0.1, 0.15, 0.2).
#' @param correction Type of edge correction. Defaults to translational.
#' @param mark1 Variable used to mark the points in the point pattern object for the first type. Default is "immune1".
#' @param mark2 Variable used to mark the points in the point pattern object for the second type. Default is "immune2".
#'
#' @returns
#' A dataframe with the following columns:
#' \describe{
#'  \item{r}{The radius at which K was calculated.}
#'  \item{k}{The observed K value}
#'  \item{theo_csr}{The theoretical K under CSR}
#'  \item{kamp_csr}{The adjusted CSR representing the KAMP permuted expectation.}
#'  \item{kamp}{The difference between observed K and KAMP CSR}
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
#' # one image from the ovarian cancer data in the package
#' ids <- unique(ovarian_df$sample_id)
#' img <- ovarian_df[ovarian_df$sample_id == ids[4], ]
#' pp <- spatstat.geom::ppp(img$x, img$y,
#'                          window = spatstat.geom::owin(range(img$x), range(img$y)),
#'                          marks = img$phenotype)
#'
#' kamp_expectation_biv(pp, rvals = c(0, 25, 50, 75, 100),
#'                      mark1 = "cytotoxic t cell", mark2 = "helper t cell")
kamp_expectation_biv <- function(ppp_obj,
                             rvals = c(0, .05, .075, .1, .15, .2),
                             correction = "trans",
                             mark1 = "immune1",
                             mark2 = "immune2") {

  # gets number of points
  npts  <- npoints(ppp_obj)

  if (npts > 10000) {
    warning("Point pattern has more than 10,000 points. Switching to border correction")
    #correction = "border"
  }

  # Kest/Kcross need r to start at 0 and increase, so evaluate on sorted
  # radii (with 0 added) and map back to rvals
  rvals_sorted <- sort(unique(c(0, rvals)))

  # spatstat bins distances when r is unevenly spaced, which makes the border
  # estimate approximate (exact with evenly spaced r)
  if (correction == "border" &&
      !isTRUE(all.equal(diff(rvals_sorted), rep(rvals_sorted[2], length(rvals_sorted) - 1)))) {
    warning("rvals are not evenly spaced: border-corrected K from spatstat is approximate. ",
            "Use evenly spaced rvals for exact values.")
  }

  k_orig = Kcross(ppp_obj, i = mark1, j = mark2,
             r = rvals_sorted,
             correction = correction)

  kamp_df = Kest(ppp_obj,
              r = rvals_sorted,
              correction = correction) %>%
    as_tibble()


  if (correction == "trans") {
    kamp_df = kamp_df %>%
      mutate(kamp_csr = trans,
             k = k_orig$trans,
             theo_csr = k_orig$theo,
             kamp = k - kamp_csr) %>%
      select(r, k, theo_csr, kamp_csr, kamp)
  } else if (correction == "iso") {
    kamp_df = kamp_df %>%
      mutate(kamp_csr = iso,
             k = k_orig$iso,
             theo_csr = k_orig$theo,
             kamp = k - kamp_csr) %>%
      select(r, k, theo_csr, kamp_csr, kamp)
  } else if (correction == "border") {
    kamp_df = kamp_df %>%
      mutate(kamp_csr = border,
             k = k_orig$border,
             theo_csr = k_orig$theo,
             kamp = k - kamp_csr) %>%
      select(r, k, theo_csr, kamp_csr, kamp)
  } else if (correction == "none") {
    kamp_df = kamp_df %>%
      mutate(kamp_csr = un, # spatstat names the uncorrected estimate "un"
             k = k_orig$un,
             theo_csr = k_orig$theo,
             kamp = k - kamp_csr) %>%
      select(r, k, theo_csr, kamp_csr, kamp)
  }

  return(kamp_df[match(rvals, rvals_sorted), ])
}
