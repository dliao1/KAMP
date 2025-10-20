#' Helper function for KAMP Variance
#' @title kamp_variance_helper
#'
#' @description Helper function to calculate the KAMP variance for a point pattern object and single radius.
#'
#' @param ppp_obj A point pattern object "ppp" from the `spatstat` package.
#' @param rvalue A single radius
#' @param correction Type of edge correction. Defaults to translational.
#' @param mark1 Value used to mark the points in the point pattern object. Default is "immune".
#'
#' @returns
#' A single-row dataframe with the following columns:
#' \describe{
#'   \item{r}{The current radius at which K was calculated.}
#'   \item{k}{The observed K value}
#'   \item{theo_csr}{The theoretical K under CSR}
#'   \item{kamp_csr}{The adjusted CSR representing the KAMP permuted expectation.}
#'   \item{kamp}{The difference between observed K and KAMP CSR}
#'   \item{var}{Variance of K under the permutation null distribution}
#'   \item{pval}{P-value}
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
                                mark1 = "immune") {
  npts = npoints(ppp_obj)
  W = Window(ppp_obj)
  areaW = spatstat.geom::area(W)

  pp_df = as.data.frame(ppp_obj)
  W = as.matrix(dist(as.matrix(select(pp_df, x, y))))
  Wr <- NULL


  # TODO: implement border correction

 if (correction == "trans") {
    e = edge.Trans(ppp_obj)
    W = ifelse(W <= rvalue, 1, 0)
    diag(W) = 0
    Wr = W * e
  }


  # Compute R0 term
  R0 = sum(Wr)

  # Compute expectation
  m = sum(ppp_obj$marks == mark1)
  R1 = sum(Wr^2)
  R2 = sum(rowSums(Wr)^2) - R1
  R3 = R0^2 - 2*R1 - 4*R2

  npairs =  npts * (npts - 1)
  f1 = m*(m-1)/npts/(npts-1)
  f2 = f1*(m-2)/(npts-2)
  f3 = f2*(m-3)/(npts-3)

  Kmat = Wr[which(ppp_obj$marks == mark1),which(ppp_obj$marks == mark1)]
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
    kamp = k - kamp_csr, # difference between K and KAMP CSR
    var = var_K,
    pvalue = min(1, pval_appx)
  )

  return(result)
}


#' Checks inputs for KAMP functions
#'
#' @param df A dataframe containing the point pattern data. Will be converted into a `ppp` object.
#' @param rvals Vector of radius values at which to compute the KAMP expectation.
#' @param univariate Logical indicating whether to compute univariate KAMP expectation. Defaults to TRUE.
#' @param mark_var Column name in `df` that contains the marks for the point pattern object.
#' @param mark1 Value used to mark the points in the point pattern object.
#' @param mark2 Value used to mark the points in the point pattern object for the second type (optional, only used if `univariate` is FALSE).
#' @param variance Logical indicating whether to compute the variance of KAMP (default is FALSE).
#' @param thin Logical indicating whether to thin the point pattern before computing KAMP (default is FALSE), called KAMP-lite.
#' @param p_thin Percentage that determines how much to thin
#' @param background Value used to define the background for the point pattern object.
#' @param ...
#'
#' @returns TRUE if all inputs are valid, otherwise throws an error with a descriptive message.
#' @keywords internal
#'
check_inputs <- function(df,
                         rvals,
                         univariate,
                         correction,
                         mark_var,
                         mark1,
                         mark2,
                         variance,
                         thin,
                         p_thin,
                         background,...) {
  # Initialize ppp_obj
  ppp_obj <- NULL

  # Check if df is a dataframe
  if (!is.data.frame(df)) {
    stop("Input df must be a dataframe.")
  }

  # Checks if mark_var is supplied
  if (is.null(mark_var) || mark_var == "") {
    stop("mark_var must be supplied and cannot be NULL or empty.")
  }

  # Check if df has x and y columns
  if (!all(c("x", "y") %in% colnames(df))) {
    stop("Input dataframe must contain 'x' and 'y' columns.")
  }

  # Factorize mark_var
  if (!mark_var %in% colnames(df)) {
    stop(paste0("mark_var '", mark_var, "' not found in dataframe columns."))
  }

  df$mark_var <- as.factor(df[[mark_var]])


  # Check if mark_var has at least two unique values
  if (nlevels(df$mark_var) < 2) {
    stop("The mark_var column must have at least two unique values.")
  }

  # Convert df to ppp object
  win <- convexhull.xy(df$x, df$y)
  ppp_obj <- ppp(df$x, df$y, window = win, marks = df$mark_var)
  all_marks <- unique(marks(ppp_obj))

  message("We expect the dataframe to be a single point process. If you have multiple point processes, subset the dataframe by ID and please run KAMP separately for each process.")

  # Makes sure ppp_obj is not NULL
  if (is.null(ppp_obj)) {
    stop("The point pattern object cannot be NULL. Conversion of dataframe to ppp object failed.")
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
  if (is.null(all_marks) || length(all_marks) == 0) {
    stop("The point pattern object does not have marks.")
  }

  if (mark1 %in% all_marks == FALSE) {
    stop("mark1 is not a mark in the point pattern object.")
  }

  if (is.null(mark2) == FALSE && mark2 %in% all_marks == FALSE) {
    stop("mark2 is not a mark in the point pattern object.")
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

    if (!is.null(mark1) && !(mark1 %in% all_marks)) {
      stop(paste0("mark1 ('", mark1, "') not found in point pattern marks."))
    }

    if(is.null(mark2) == FALSE) {
      message("mark2 is not used in univariate KAMP. It will be ignored.")
    }

  } else { #must be bivariate

    if (!is.null(mark1) && !(mark1 %in% all_marks)) {
      stop(paste0("mark1 ('", mark1, "') not found in point pattern marks."))
    }

    if (!is.null(mark2) && !(mark2 %in% all_marks)) {
      stop(paste0("mark2 ('", mark2, "') not found in point pattern marks."))
    }

    if (is.null(mark1) || is.null(mark2)) {
      stop("Both mark1 and mark2 must be specified for bivariate KAMP.")
    }

    if (mark1 == mark2) {
      stop("mark1 and mark2 cannot be the same for bivariate KAMP.")
    }
  }

  if (thin == TRUE && variance == TRUE) {
    message("Variance calculation with KAMP lite is not recommended. Variance will still be computed, but interpret with caution.")
  }

  # End of most input sanitization checks

  if (sum(marks(ppp_obj) == mark1) < 5) {
    message(paste0("Less than 5 target cells marked as '", mark1, "'. This may lead to unreliable results."))
  }

  if (npoints(ppp_obj) <= 10000) {
    message("The point pattern object has more than 10000 points. Switching to border correction")
  }



  return(ppp_obj)
}
