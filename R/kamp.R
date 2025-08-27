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
#' For variance, this function utlizes a matrix-based implementation.
#'
#' See `?Kcross` and `?Kest` for more details on the K calculation methods.
#'
#' See `kamp_expectation_mat` for the matrix-based implementation (Note: currently
#' not recommended due to slow speed)
#'
#' @param ppp_obj A point pattern object of class `ppp` from the `spatstat` package.
#' @param rvals A vector of distances at which to compute the KAMP expectation and variance.
#' @param univariate A logical value indicating whether to compute univariate KAMP (default is TRUE).
#' @param mark1 Variable used to mark the points in the point pattern object for the first type.
#' @param mark2 Variable used to mark the points in the point pattern object for the second type (optional, only used if `univariate` is FALSE).
#' @param correction Type of edge correction. Defaults to translational.
#' @param variance A logical value indicating whether to compute the variance of KAMP (default is FALSE).
#' @param thin A logical value indicating whether to thin the point pattern before computing KAMP (default is FALSE), called KAMP-lite.
#' @param p_thin Percentage that determines how much to thin the amount of points in the point pattern object. Default is 0.
#' @param background Variable used to define the background for the point pattern object.
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
#' # Loads required packages
#' library(spatstat.geom)
#' library(spatstat.explore)
#'
#' # Simulates a simple marked point pattern
#' set.seed(100)
#' x_coords <- runif(100)
#' y_coords <- runif(100)
#' marks_vec <- sample(c("immune", "background"), 100, replace = TRUE)
#' win <- owin(c(0,1), c(0,1))
#' ppp_obj <- ppp(x_coords, y_coords, window = win, marks = marks_vec)
#'
#' # Defines radius values for K-function estimation
#' r_vals <- seq(0.01, 0.1, by = 0.01)
#'
#' # Computes univariate KAMP expectation
#' kamp_result <- kamp(ppp_obj = ppp_obj,
#'                     rvals = r_vals,
#'                     univariate = TRUE,
#'                     mark1 = "immune")
#' head(kamp_result)
#'
#' # Compute univariate KAMP expectation with thinning
#' kamp_thin <- kamp(ppp_obj = ppp_obj,
#'                   rvals = r_vals,
#'                   univariate = TRUE,
#'                   mark1 = "immune",
#'                   thin = TRUE,
#'                   p_thin = 0.3)
#' head(kamp_thin)
#'
#' # Use real data from VectraPolarisData in package
#' data(ovarian_df)
#' sample_id <- unique(ovarian_df$sample_id)[1]
#' ov_df <- subset(ovarian_df, sample_id == sample_id)
#' win <- convexhull.xy(ov_df$x, ov_df$y)
#' ppp_real <- ppp(ov_df$x, ov_df$y, window = win, marks = ov_df$immune)
#' kamp_real <- kamp(ppp_obj = ppp_real,
#'                   rvals = seq(0.01, 0.1, 0.01),
#'                   univariate = TRUE,
#'                   mark1 = "immune")
#' head(kamp_real)
kamp = function(df, # change to dataframe with x, y, mark_var, (factor) mark1, mark2
                         # expect this to just be one point process
                rvals,
                univariate = TRUE,
                mark_var,
                mark1,
                mark2 = NULL,
                variance = FALSE,
                correction = "trans",
                thin = FALSE,
                p_thin = 0.5,
                background = NULL,
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
                     p_thin,
                     background)

  if (is.null(ppp_obj)) {
    stop("Input checks failed and point process object could not be created.")
  }


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
  }
  else if (univariate == FALSE && variance == FALSE) {
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
