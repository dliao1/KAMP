#' Title KAMP Variance
#' @title get_kamp_variance
#' @description Wrapper function for mapping across multiple radii, believe furrr is parallelized version but not sure i should do that
#' @param ppp_obj A point pattern object of class "ppp" from the spatstat package.
#' @param rvec A vector of radii at which to calculate the KAMP variance.
#' @param markvar The variable used to mark the points in the point pattern object. Default is "immune".
#'
#' @returns A tibble containing the radius, K statistic, expectation, variance, Z statistic, p-value, and time taken to compute the variance.
#' @export
#'
#' @examples
get_kamp_variance = function(ppp_obj, rvec, markvar = "immune") {
  return(map_dfr(rvec, get_kamp_variance_helper, ppp_obj = ppp_obj, .progress = TRUE))
}
