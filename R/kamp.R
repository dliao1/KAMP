#' Title
#'
#' @param ppp_obj
#' @param rvals
#' @param univariate
#' @param marksvar1
#' @param marksvar2
#' @param correction
#' @param variance
#' @param thin
#' @param p_thin
#' @param background
#' @param ...
#'
#' @returns
#' @export
#'
#' @examples
kamp = function(ppp_obj, rvals, univariate = TRUE, marksvar1, marksvar2 = NULL, correction = "trans", variance = FALSE, thin = FALSE, p_thin = 0.5,
                background = NULL,...){

  input_checks <- check_inputs(ppp_obj,
                     rvals,
                     univariate,
                     correction,
                     marksvar1,
                     marksvar2,
                     variance,
                     thin,
                     p_thin,
                     background) # will this stop on its own if any inputs fail?

  if (thin == TRUE) {
    ppp_obj = rthin(ppp_obj, 1 - p_thin)
  }

  if (univariate == TRUE && variance == FALSE) {
    results <- kamp_expectation(ppp_obj = ppp_obj,
                                rvals = rvals,
                                correction = correction,
                                marksvar1 = marksvar1)
  } else if (univariate == TRUE && variance == TRUE) {
    results <- kamp_variance(ppp_obj = ppp_obj,
                             rvals = rvals,
                             correction = correction,
                             marksvar1 = marksvar1)
  }
  #else if (univariate == FALSE && variance == FALSE) {
  #  results <- kamp_expectation_biv(ppp_obj = ppp_obj,
  #                                  rvals = rvals,
  #                                  marksvar1 = marksvar1,
  #                                  marksvar2 = marksvar2,
  #                                  correction = correction)
  #} else if (univariate == FALSE && variance == TRUE) {
  #  results <- kamp_variance_biv(ppp_obj = ppp_obj,
  #                               rvals = rvals,
  #                               marksvar1 = marksvar1,
  #                               marksvar2 = marksvar2,
  #                               correction = correction)
  #}
  return(results)

}
