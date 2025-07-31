#' KAMP Expectation Helper
#' @title kamp_expectation_mat_helper
#'
#' @description Helper function to calculate the KAMP expectation for a given point pattern object and radius.
#'
#' @param ppp_obj A point pattern object "ppp" from the spatstat package.
#' @param rvalue A single radius
#' @param correction Type of edge correction. Defaults to translational
#' @param markvar The variable used to mark the points in the point pattern object. Defaults to "immune".
#'
#' @returns
#' A dataframe with the following columns:
#' \describe{
#'   \item{r}{The radius at which K was calculated.}
#'   \item{k}{The observed K value from `Kcross`}
#'   \item{theo_csr}{The theoretical K under CSR from `Kcross`}
#'   \item{kamp_csr}{The adjusted CSR representing the permuted expectation.}
#'   \item{kamp_fundiff}{The difference between observed K and KAMP CSR}
#' }
#'
#' @importFrom spatstat.explore Kcross Kest edge.Trans
#' @importFrom spatstat.geom area.owin ppp as.owin npoints Window
#' @importFrom dplyr mutate select
#' @importFrom tibble as_tibble
#' @importFrom magrittr %>%
#' @importFrom purrr map_dfr
#' @importFrom stats dist pnorm
#' @importFrom tibble tibble
#' @keywords internal
kamp_expectation_mat_helper = function(ppp_obj,
                                       rvalue,
                                       correction = "trans",
                                       markvar = "immune") {

  npts = npoints(ppp_obj)
  W = Window(ppp_obj)
  areaW = spatstat.geom::area(W)

  pp_df = as.data.frame(ppp_obj)
  W = as.matrix(dist(as.matrix(select(pp_df, x, y))))
  Wr <- NULL

  # TODO: implement border correction
  # Bottleneck here for edge correction when lots of points
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
  m = sum(ppp_obj$marks == markvar)

  Kmat = Wr[which(ppp_obj$marks == markvar),which(ppp_obj$marks == markvar)]
  K = areaW * sum(Kmat)/ m / (m - 1) # Ripley's K
  mu_K = areaW * R0 / npts / (npts - 1) # expectation


  df <- tibble(
    r = rvalue,
    k = K,
    theo_csr = pi * rvalue^2, # theoretical CSR using area of circle
    kamp_csr = mu_K,
    kamp_fundiff = k - kamp_csr # difference between K and CSR
  )

  return(df)
}

#' Helper function for KAMP Variance
#' @title kamp_variance_helper
#'
#' @description Helper function to calculate the KAMP variance for a point pattern object and single radius.
#'
#' @param ppp_obj A point pattern object "ppp" from the spatstat package.
#' @param rvalue A single radius
#' @param correction Type of edge correction. Defaults to translational.
#' @param markvar The variable used to mark the points in the point pattern object. Defaults to "immune".
#'
#' @returns
#' A single-row dataframe with the following columns:
#' \describe{
#'   \item{r}{The current radius at which K was calculated.}
#'   \item{k}{The observed K value}
#'   \item{theo_csr}{The theoretical K under CSR}
#'   \item{kamp_csr}{The adjusted CSR representing the KAMP permuted expectation.}
#'   \item{kamp_fundiff}{The difference between observed K and KAMP CSR}
#'   \item{var}{Variance of K under the permutation null distribution}
#'   \item{z}{Z statistic, calculated by normalizing K using the formula: (K - KAMP)/sqrt(var)}
#'   \item{pval}{P-value, calculated using the formula: pnorm(-z)}
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
#' @keywords internal
kamp_variance_helper = function(ppp_obj,
                                rvalue,
                                correction = "trans",
                                markvar = "immune") {
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
  m = sum(ppp_obj$marks == markvar)
  R1 = sum(Wr^2)
  R2 = sum(rowSums(Wr)^2) - R1
  R3 = R0^2 - 2*R1 - 4*R2

  npairs =  npts * (npts - 1)
  f1 = m*(m-1)/npts/(npts-1)
  f2 = f1*(m-2)/(npts-2)
  f3 = f2*(m-3)/(npts-3)

  Kmat = Wr[which(ppp_obj$marks == markvar),which(ppp_obj$marks == markvar)]
  K = areaW * sum(Kmat)/ m / (m - 1) # Ripley's K
  mu_K = areaW * R0 / npts / (npts - 1) # expectation
  var_K = areaW^2 * (2 * R1 * f1 + 4 * R2 * f2 + R3 * f3) / m / m / (m - 1) / (m - 1) - mu_K^2   # variance


  Z_k = (K - mu_K) / sqrt(var_K) # Test statistic
  pval_appx = pnorm(-Z_k) # approximated p-value based on normal distribution


  result = tibble(
    r = rvalue,
    k = K,
    theo_csr = pi * rvalue^2, # theoretical CSR using area of circle
    kamp_csr = mu_K, # K expectation under permutation distributions
    var = var_K,
    z = Z_k, # test statistic
    pvalue = min(1, pval_appx)
  )

  return(result)
}

#' Checks the validity of inputs for univariate KAMP functions
#'
#' @param ppp_obj A point pattern object "ppp" from the `spatstat` package.
#' @param rvec Vector of radii
#' @param correction Type of edge correction.
#' @param markvar The variable used to mark the points in the point pattern object.
#' @param thin_pct Thinning percentage
#'
#' @returns TRUE only if all the parameter checks pass
#' @keywords internal
#'
check_valid_inputs_univ <- function(ppp_obj,
                                    rvec,
                                    correction,
                                    markvar,
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
  if (markvar %in% levels(ppp_obj$marks) == FALSE) {
    stop("markvar is not a mark in the point pattern object.")
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


#' Title
#'
#' @param ppp_obj
#' @param rvals
#' @param univariate
#' @param marksvar1
#' @param marksvar2
#' @param variance
#' @param thin
#' @param p_thin
#' @param background
#' @param ...
#'
#' @returns
#' @keywords internal
#'
check_inputs <- function(ppp_obj,
                         rvals,
                         univariate,
                         correction,
                         marksvar1,
                         marksvar2,
                         variance,
                         thin,
                         p_thin,
                         background,...) {

  # Makes sure ppp_obj is not NULL
  if (is.null(ppp_obj)) {
    stop("The point pattern object cannot be NULL.")
  }

  # Check if rvec is numeric
  if (!is.numeric(rvals) || any(rvals < 0)) {
    stop("rvec must be numeric and 0 or more")
  }

  # Check if correction is translational or isotropic - default to trans maybe?
  if (correction %in% c("trans", "iso") == FALSE) {
    stop("Currently only isotropic and translational edge correction are supported.")
  }

  # Ensure there are marks
  if (is.null(ppp_obj$marks)) {
    stop("The point pattern object must have marks.")
  }

  if (!is.logical(thin)) {
    stop("Argument 'thin' must be TRUE or FALSE.")
  }

  if (thin == TRUE) {
    # Check if p_thin is numeric
    if (!is.numeric(p_thin)) {
      stop("p_thin must be numeric.")
    }

    # Check if p_thin is between 0 and 1
    if (p_thin < 0 || p_thin > 1) {
      stop("p_thin must be between 0 and 1.")
    }

    # Lets user know that variance = TRUE with KAMP lite is not supported
    if (variance == TRUE) {
      message("Variance calculation is not supported with KAMP lite")
    }
  }


  if (univariate == TRUE) {
    all_marks <- unique(marks(ppp_obj))

    if (!is.null(marksvar1) && !(marksvar1 %in% all_marks)) {
      stop(paste0("marksvar1 ('", marksvar1, "') not found in point pattern marks."))
    }
  } else { #must be bivariate
    all_marks <- unique(marks(ppp_obj))

    if (!is.null(marksvar1) && !(marksvar1 %in% all_marks)) {
      stop(paste0("marksvar1 ('", marksvar1, "') not found in point pattern marks."))
    }

    if (!is.null(marksvar2) && !(marksvar2 %in% all_marks)) {
      stop(paste0("marksvar2 ('", marksvar2, "') not found in point pattern marks."))
    }

    if (is.null(marksvar1) || is.null(marksvar2)) {
      stop("Both marksvar1 and marksvar2 must be specified for bivariate KAMP.")
    }
  }

  return(TRUE)
}
