#' Computes KAMP expectation and variance
#'
#' @title KAMP
#'
#' @description
#' This function computes the KAMP expectation and variance for a given point pattern.
#' Calculates Ripley's K using both the traditional Ripley's K method
#' (based on `Kcross`) and the KAMP-adjusted CSR baseline (based on `Kest`).
#'
#' The KAMP-adjusted CSR represents a more robust baseline for K that accounts for spatial clustering or inhomogeneity
#' in a point pattern compared to the traditional CSR assumption, while
#' avoiding the computational burden of permuting the point pattern.
#'
#' For expectation, this function uses the `spatstat` package under the hood, which
#' automatically uses border correction when the number of
#' points in the point pattern is more than 3000.
#'
#' For variance, this function uses the Rcpp-accelerated implementation.
#'
#' See `?Kcross` and `?Kest` for more details on the K calculation methods.
#'
#' @param df Either a point pattern object of class `ppp` from the `spatstat` package, or a
#' data.frame with `x`/`y` columns and a marks column named by `mark_var` (used to build a `ppp`
#' internally via a convex-hull window). Should contain one point process at a time.
#' @param rvals A vector of distances at which to compute the KAMP expectation and variance.
#' @param univariate A logical value indicating whether to compute univariate KAMP (default is TRUE).
#' @param mark_var Column name in `df` containing the marks. Required when `df` is a data.frame
#' (an error is thrown if it is left as `NULL`); ignored when `df` is already a `ppp` object.
#' @param mark1 Variable used to mark the points in the point pattern object for the first type.
#' @param mark2 Variable used to mark the points in the point pattern object for the second type (optional, only used if `univariate` is FALSE).
#' @param correction Type of edge correction. Defaults to translational.
#' @param variance A logical value indicating whether to compute the variance of KAMP (default is FALSE).
#' @param thin A logical value indicating whether to thin the point pattern before computing KAMP (default is FALSE), called KAMP-lite.
#' @param p_thin Proportion of points to remove when `thin = TRUE`. Default is 0.5.
#' @param ... Additional arguments passed to the underlying functions.
#'
#' @returns
#' A dataframe with the following columns:
#' \describe{
#'  \item{r}{The radius at which K was calculated.}
#'  \item{k}{The observed K value}
#'  \item{theo_csr}{The theoretical K under CSR}
#'  \item{kamp_csr}{The adjusted CSR representing the KAMP expectation.}
#'  \item{kamp}{The difference between observed K and KAMP CSR}
#'  \item{var}{If variance = TRUE, variance of K under the KAMP null distribution}
#'  \item{pval}{If variance = TRUE, p-value that estimates the probability of observing a deviation
#'  from the expected KAMP-adjusted value as large or larger than the one
#'  observed, under the null hypothesis of CSR). Calculated using a normal approximation.}
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
#' @examples
#' # one image from the ovarian cancer data in the package
#' ids <- unique(ovarian_df$sample_id)
#' img <- ovarian_df[ovarian_df$sample_id == ids[4], ]
#' r_vals <- c(0, 25, 50, 75, 100)
#'
#' # univariate, marks column given by mark_var
#' kamp(img, rvals = r_vals, mark_var = "immune", mark1 = "immune")
#'
#' # with variance
#' kamp(img, rvals = r_vals, mark_var = "immune", mark1 = "immune", variance = TRUE)
#'
#' # bivariate
#' kamp(img, rvals = r_vals, univariate = FALSE, mark_var = "phenotype",
#'      mark1 = "cytotoxic t cell", mark2 = "helper t cell")
#'
#' # KAMP-lite: randomly remove 30% of points first
#' kamp(img, rvals = r_vals, mark_var = "immune", mark1 = "immune",
#'      thin = TRUE, p_thin = 0.3)
kamp = function(df, # expect this to just be one point process
                rvals,
                univariate = TRUE,
                mark_var = NULL,
                mark1,
                mark2 = NULL,
                variance = FALSE,
                correction = "trans",
                thin = FALSE,
                p_thin = 0.5,
                ...){

  ppp_obj  <- check_inputs(df,
                     rvals,
                     univariate,
                     correction,
                     mark_var,
                     mark1,
                     mark2,
                     variance,
                     thin,
                     p_thin)

  if (is.null(ppp_obj)) {
    stop("Input checks failed and point process object could not be created.")
  }

  # normalize full-name correction aliases to the short form the functions below expect
  correction <- switch(correction,
                        "translational" = "trans",
                        "isotropic" = "iso",
                        correction)

  if (thin == TRUE) {
    ppp_obj = rthin(ppp_obj, 1 - p_thin)
  }

  results <- NULL
  if (univariate == TRUE && variance == FALSE) {
    results <- kamp_expectation(ppp_obj = ppp_obj,
                                rvals = rvals,
                                correction = correction,
                                mark1 = mark1)
  } else if (univariate == TRUE && variance == TRUE) {
    results <- kamp_variance(ppp_obj = ppp_obj,
                                  rvals = rvals,
                                  correction = correction,
                                  mark1 = mark1)
  } else if (univariate == FALSE && variance == FALSE) {
    results <- kamp_expectation_biv(ppp_obj = ppp_obj,
                                    rvals = rvals,
                                    mark1 = mark1,
                                    mark2 = mark2,
                                    correction = correction)
  } else if (univariate == FALSE && variance == TRUE) {
    results <- kamp_variance_biv(ppp_obj = ppp_obj,
                                      rvals = rvals,
                                      mark1 = mark1,
                                      mark2 = mark2,
                                      correction = correction)
  }
  return(results)

}
