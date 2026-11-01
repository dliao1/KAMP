#' Computes GAMP expectation and variance
#'
#' @title GAMP
#'
#' @description
#' This function computes the GAMP expectation and variance for a given point pattern.
#' GAMP is the nearest-neighbor G analogue of KAMP: it gives the exact mean (and
#' optionally variance) of the uncorrected G function under random relabeling of
#' marks, i.e. the quantities a permutation test would estimate, without running
#' any permutations.
#'
#' Univariate G is the nearest-neighbor distance distribution of the `mark1` points
#' (matches `Gest(..., correction = "none")`). Bivariate G is the distribution of
#' the distance from each `mark1` point to its nearest `mark2` point (matches
#' `Gcross(..., correction = "none")`). Only no edge correction is supported.
#'
#' Neighbor counts are computed from `spatstat.geom::closepairs()`, so only pairs
#' within `max(rvals)` are touched. For variance, this function uses the
#' Rcpp-accelerated implementation.
#'
#' @param df Either a point pattern object of class `ppp` from the `spatstat` package, or a
#' data.frame with `x`/`y` columns and a marks column named by `mark_var` (used to build a `ppp`
#' internally via a convex-hull window). Should contain one point process at a time.
#' @param rvals A vector of distances at which to compute the GAMP expectation and variance.
#' @param univariate A logical value indicating whether to compute univariate GAMP (default is TRUE).
#' @param mark_var Column name in `df` containing the marks. Required when `df` is a data.frame
#' (an error is thrown if it is left as `NULL`); ignored when `df` is already a `ppp` object.
#' @param mark1 Variable used to mark the points in the point pattern object for the first type.
#' @param mark2 Variable used to mark the points in the point pattern object for the second type (optional, only used if `univariate` is FALSE).
#' @param variance A logical value indicating whether to compute the variance of GAMP (default is FALSE).
#' @param correction Type of edge correction. Only "none" (no edge correction) is currently supported.
#' @param thin A logical value indicating whether to thin the point pattern before computing GAMP (default is FALSE).
#' @param p_thin Proportion of points to remove when `thin = TRUE`. Default is 0.5.
#' @param ... Additional arguments passed to the underlying functions.
#'
#' @returns
#' A dataframe with the following columns:
#' \describe{
#'  \item{r}{The radius at which G was calculated.}
#'  \item{g}{The observed G value}
#'  \item{theo_csr}{The theoretical G under CSR, `1 - exp(-lambda * pi * r^2)`, where
#'  lambda is the intensity of the `mark1` (univariate) or `mark2` (bivariate) points}
#'  \item{gamp_csr}{The adjusted CSR representing the GAMP expectation.}
#'  \item{gamp}{The difference between observed G and GAMP CSR}
#'  \item{var}{If variance = TRUE, variance of G under the GAMP null distribution}
#'  \item{pvalue}{If variance = TRUE, p-value for G being larger than expected
#'  (clustering), calculated using a normal approximation: pnorm(-z).}
#' }
#'
#' @importFrom spatstat.random rthin
#' @export
#' @examples
#' set.seed(100)
#' win <- spatstat.geom::owin(c(0, 1), c(0, 1))
#' pp <- spatstat.random::rpoispp(lambda = 150, win = win)
#' marks <- sample(c("immune1", "immune2", "background"), pp$n, replace = TRUE)
#' marked_pp <- spatstat.geom::ppp(pp$x, pp$y, window = win, marks = factor(marks))
#'
#' # univariate
#' gamp_result <- gamp(marked_pp, rvals = seq(0, 0.1, 0.01),
#'                     univariate = TRUE, mark1 = "immune1", variance = TRUE)
#' head(gamp_result)
#'
#' # bivariate
#' gamp_biv <- gamp(marked_pp, rvals = seq(0, 0.1, 0.01),
#'                  univariate = FALSE, mark1 = "immune1", mark2 = "immune2")
#' head(gamp_biv)
gamp = function(df, # expect this to just be one point process
                rvals,
                univariate = TRUE,
                mark_var = NULL,
                mark1,
                mark2 = NULL,
                variance = FALSE,
                correction = "none",
                thin = FALSE,
                p_thin = 0.5,
                ...){

  ppp_obj <- check_inputs(df,
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

  if (thin == TRUE) {
    ppp_obj = rthin(ppp_obj, 1 - p_thin)
  }

  results <- NULL
  if (univariate == TRUE && variance == FALSE) {
    results <- gamp_expectation(ppp_obj = ppp_obj,
                                rvals = rvals,
                                correction = correction,
                                mark1 = mark1)
  } else if (univariate == TRUE && variance == TRUE) {
    results <- gamp_variance(ppp_obj = ppp_obj,
                             rvals = rvals,
                             correction = correction,
                             mark1 = mark1)
  } else if (univariate == FALSE && variance == FALSE) {
    results <- gamp_expectation_biv(ppp_obj = ppp_obj,
                                    rvals = rvals,
                                    correction = correction,
                                    mark1 = mark1,
                                    mark2 = mark2)
  } else if (univariate == FALSE && variance == TRUE) {
    results <- gamp_variance_biv(ppp_obj = ppp_obj,
                                 rvals = rvals,
                                 correction = correction,
                                 mark1 = mark1,
                                 mark2 = mark2)
  }
  return(results)

}
