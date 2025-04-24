#' Title KAMP Expectation
#' @title get_kamp_expectation
#' @description Calculates the KAMP expectation for a given point pattern object.
#'
#' @param ppp_obj A point pattern object of class "ppp" from the spatstat package.
#' @param rvec A vector of radii at which to calculate the KAMP expectation.
#' @param correction Type of border correction
#' @param markvar The variable used to mark the points in the point pattern object. Default is "immune".
#'
#' @returns A tibble containing the sample_id, r values, K statistic, theoretical CSR, KAMP CSR, FUNDIFF
#' @importFrom spatstat.explore Kcross Kest
#' @importFrom spatstat.geom area.owin ppp as.owin
#' @importFrom dplyr mutate select
#' @importFrom tibble as_tibble
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
get_kamp_expectation <- function(ppp_obj,
                                 rvec = c(0, .05, .075,.1, .15, .2),
                                 correction = "trans",
                                 markvar = "immune") { # not sure if markvar is necessary?

  # Pre-existing code that uses spatstat
  # Gets original K using translational correction
  k = Kcross(ppp_obj, i = markvar, j = markvar,
             r = rvec,
             correction = correction)


  kamp = Kest(ppp_obj,
              r = rvec,
              correction = correction) %>%
    as_tibble() %>%
    mutate(kamp_csr = trans,
           k = k$trans, # takes calculated K from Kcross
           theo_csr = k$theo,
           kamp_fundiff = k - kamp_csr) %>% # takes calculated theoretical K from Kcross
    select(r, theo_csr, kamp_csr, k, kamp_fundiff) # gets r, kamp CSR, K statistic, method

  return(kamp)
}


#' Title
#' @title get_kamp_expectation_mat
#' @description Wrapper function for mapping across multiple radii to calculate the KAMP expectation.
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
