#' Helper function for KAMP Expectation (Bivariate matrix-based implementation)
#'
#' @title kamp_expectation_biv_mat_helper
#'
#' @description
#' Helper function to calculate the KAMP expectation for bivariate point patterns.
#'
#' @param ppp_obj A point pattern object "ppp" from the `spatstat` package.
#' @param rvalue A single radius
#' @param correction Type of edge correction. Defaults to translational.
#' @param marksvar1 Variable used to mark the points in the point pattern object for the first type. Default is "immune1".
#' @param marksvar2 Variable used to mark the points in the point pattern object for the second type. Default is "immune2".
#'
#' @returns
#' A single-row dataframe with the following columns:
#' \describe{
#'  \item{r}{The radius at which K was calculated.}
#'  \item{k}{The observed K value}
#'  \item{theo_csr}{The theoretical K under CSR}
#'  \item{kamp_csr}{The adjusted CSR representing the KAMP permuted expectation.}
#'  \item{kamp}{The difference between observed K and KAMP CSR}
#' }
#'
#' @importFrom spatstat.explore Kcross Kest edge.Trans edge.Ripley
#' @importFrom spatstat.geom area.owin ppp as.owin npoints Window
#' @importFrom dplyr mutate select
#' @importFrom tibble as_tibble
#' @importFrom magrittr %>%
#' @importFrom purrr map_dfr
#' @importFrom tictoc tic toc
#' @importFrom stats dist pnorm
#' @importFrom tibble tibble
#'
#' @keywords internal
#'
kamp_expectation_biv_mat_helper <- function(ppp_obj,
                                            rvalue,
                                            correction = "trans",
                                            marksvar1 = "immune1",
                                            marksvar2 = "immune2") {

  npts = npoints(ppp_obj)
  W = Window(ppp_obj)
  areaW = spatstat.geom::area(W)

  pp_df = as.data.frame(ppp_obj)
  W = as.matrix(dist(as.matrix(select(pp_df, x, y))))
  Wr <- NULL

  # TODO: implement border correction
  if (correction == "iso") {
    e = edge.Ripley(ppp_obj)
    W = ifelse(W <= rvalue, 1, 0)
    diag(W) = 0
    Wr = W * e
  } else if (correction == "trans") {
    e = edge.Trans(ppp_obj)
    W = ifelse(W <= rvalue, 1, 0)
    diag(W) = 0
    Wr = W * e
  }

  # Compute R0 term
  R0 = sum(Wr)

  # Compute expectation
  m1 = sum(ppp_obj$marks == marksvar1)
  m2 = sum(ppp_obj$marks == marksvar2)

  Kmat = Wr[which(ppp_obj$marks == marksvar1),which(ppp_obj$marks == marksvar2)]
  K = areaW * sum(Kmat) / m1 / m2 # Ripley's K based on translation correction
  mu_K = areaW * R0/npts / (npts-1) # expectation

  df <- tibble(
    r = rvalue,
    k = K,
    theo_csr = pi * rvalue^2, # theoretical CSR using area of circle
    kamp_csr = mu_K,
    kamp = k - kamp_csr # difference between K and CSR
  )

  return(df)

}

#' Helper function for bivariate KAMP Variance
#'
#' @title kamp_variance_biv_helper
#' @description
#' Helper function to calculate the KAMP variance for bivariate point patterns.
#'
#' @param ppp_obj A point pattern object "ppp" from the `spatstat` package.
#' @param rvalue A single radius
#' @param correction Type of edge correction. Defaults to translational.
#' @param marksvar1 Variable used to mark the points in the point pattern object for the first type. Default is "immune1".
#' @param marksvar2 Variable used to mark the points in the point pattern object for the second type. Default is "immune2".
#'
#' @returns
#' A single-row dataframe with the following columns:
#' \describe{
#'  \item{r}{The radius at which K was calculated.}
#'  \item{k}{The observed K value}
#'  \item{theo_csr}{The theoretical K under CSR}
#'  \item{kamp_csr}{The adjusted CSR representing the KAMP permuted expectation.}
#'  \item{var}{Variance of K under the permutation null distribution}
#'  \item{pval}{P-value, calculated using the formula: pnorm(-z)}
#'  }
#'
#' @importFrom spatstat.explore Kcross Kest edge.Trans edge.Ripley
#' @importFrom spatstat.geom area.owin ppp as.owin npoints Window
#' @importFrom dplyr mutate select
#' @importFrom tibble as_tibble
#' @importFrom magrittr %>%
#' @importFrom purrr map_dfr
#' @importFrom tictoc tic toc
#' @importFrom stats dist pnorm
#' @importFrom tibble tibble
#'
#' @keywords internal
#'
kamp_variance_biv_helper <- function(ppp_obj,
                                     rvalue,
                                     correction = "trans",
                                     marksvar1 = "immune1",
                                     marksvar2 = "immune2") {

  npts = npoints(ppp_obj)
  W = Window(ppp_obj)
  areaW = spatstat.geom::area(W)

  pp_df = as.data.frame(ppp_obj)
  W = as.matrix(dist(as.matrix(select(pp_df, x, y))))
  Wr <- NULL

  # TODO: implement border correction
  if (correction == "iso") {
    e = edge.Ripley(ppp_obj)
    W = ifelse(W <= rvalue, 1, 0)
    diag(W) = 0
    Wr = W * e
  } else if (correction == "trans") {
    e = edge.Trans(ppp_obj)
    W = ifelse(W <= rvalue, 1, 0)
    diag(W) = 0
    Wr = W * e
  }

  R0 = sum(Wr)
  R1 = sum(Wr^2)
  R2 = sum(rowSums(Wr)^2) - R1
  R3 = R0^2 - 2*R1 - 4*R2

  m1 = sum(ppp_obj$marks == marksvar1)
  m2 = sum(ppp_obj$marks == marksvar2)
  f1 = m1*m2/npts/(npts-1)
  f2 = f1*(m1+m2-2)/(npts-2)
  f3 = f1*(m1-1)*(m2-1)/(npts-2)/(npts-3)


  Kmat = Wr[which(ppp_obj$marks == marksvar1),which(ppp_obj$marks == marksvar2)]
  K = areaW * sum(Kmat) / m1 / m2 # Ripley's K based on translation correction
  mu_K = areaW * R0/npts / (npts-1) # expectation
  var_K = areaW^2*(R1*f1 + R2*f2 + R3*f3)/m1/m1/m2/m2 - mu_K^2

  Z_k = (K-mu_K) / sqrt(var_K) # Test statistic
  pval_appx = pnorm(-Z_k) # approximated p-value based on normal distribution

  result = tibble(
    r = rvalue,
    k = K,
    theo_csr = pi * rvalue^2, # theoretical CSR using area of circle
    kamp_csr = mu_K, # K expectation under permutation distributions
    var = var_K,
    pvalue = min(1, pval_appx)
  )

  return(result)

}


#' Checks the validity of inputs for bivariate KAMP functions
#'
#' @param ppp_obj A point pattern object "ppp" from the `spatstat` package.
#' @param rvec Vector of radii
#' @param correction Type of edge correction.
#' @param marksvar1 Variable used to mark the points in the point pattern object for the first type
#' @param marksvar2 Variable used to mark the points in the point pattern object for the second type
#' @param thin_pct Thinning percentage
#'
#' @returns TRUE only if all the parameter checks pass
#' @keywords internal
#'
check_valid_inputs_biv <- function(ppp_obj,
                                   rvec,
                                   correction,
                                   marksvar1,
                                   marksvar2,
                                   thin_pct) {
  # Makes sure ppp_obj is not NULL
  if (is.null(ppp_obj)) {
    stop("The point pattern object cannot be NULL.")
  }

  # Check if rvec is numeric
  if (!is.numeric(rvec)) {
    stop("rvec must be numeric.")
  }

  # Check if correction is translational or isotropic
  if (correction %in% c("trans", "iso") == FALSE) {
    stop("Currently only isotropic and translational edge correction are supported.")
  }

  # Ensure there are marks
  if (is.null(ppp_obj$marks)) {
    stop("The point pattern object must have marks.")
  }

  # Check if the specified marks are present
  if (marksvar1 %in% levels(ppp_obj$marks) == FALSE) {
    stop("marksvar1 is not a mark in the point pattern object.")
  }

  if (marksvar2 %in% levels(ppp_obj$marks) == FALSE) {
    stop("marksvar2 is not a mark in the point pattern object.")
  }

  # Check if the specified marks are different
  if (marksvar1 == marksvar2) {
    stop("marksvar1 and marksvar2 must be different.")
  }

  # Check if thin_pct is numeric
  if (!is.numeric(thin_pct)) {
    stop("thin_pct must be numeric.")
  }

  # Check if thin_pct is between 0 and 1
  if (thin_pct < 0 || thin_pct > 1) {
    stop("thin_pct must be between 0 and 1.")
  }

  return(TRUE)

}

