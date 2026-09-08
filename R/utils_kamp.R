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
#' @importFrom stats dist pnorm fft
#' @importFrom tibble tibble
#' @keywords internal
kamp_variance_helper = function(ppp_obj,
                                rvalue,
                                correction = "trans",
                                mark1 = "immune") {
  npts = npoints(ppp_obj)
  ppp_window = Window(ppp_obj)
  areaW = spatstat.geom::area(ppp_window)
  marks_vec <- as.character(ppp_obj$marks)
  m <- sum(marks_vec == mark1)

  if (correction %in% c("trans", "translational", "border", "none")) {
    # closepairs-based -- only touches pairs within rvalue, never builds
    # the full n x n distance/weight matrix
    cp <- spatstat.geom::closepairs(ppp_obj, rmax = rvalue, what = "all", twice = TRUE)
    i_idx <- cp$i
    j_idx <- cp$j

    if (correction == "none") {
      # no edge correction at all: every close pair counts with weight 1
      e_vals <- rep(1, length(i_idx))
      denom_m <- m - 1
      denom_npts <- npts - 1

    } else if (correction %in% c("trans", "translational")) {
      if (ppp_window$type == "rectangle") {
        W_width  <- diff(ppp_window$xrange)
        W_height <- diff(ppp_window$yrange)
        area_int <- pmax(0, W_width - abs(cp$dx)) * pmax(0, W_height - abs(cp$dy))
        e_vals   <- ifelse(area_int > 0, areaW / area_int, 0)
      } else {
        # polygonal/mask window -- same pixel approximation edge.Trans uses
        # internally, via a single 2D FFT autocorrelation of the mask, then
        # a per-pair lookup by displacement instead of an n x n matrix
        W_mask <- spatstat.geom::as.mask(ppp_window)
        m_num  <- as.numeric(W_mask$m)
        ny <- nrow(W_mask$m); nx <- ncol(W_mask$m)
        xstep <- W_mask$xstep; ystep <- W_mask$ystep

        F_m <- apply(matrix(m_num, ny, nx), 2, fft)
        F_m <- t(apply(F_m, 1, fft))

        ac_step <- apply(Mod(F_m)^2, 2, fft, inverse = TRUE) / ny
        ac <- Re(t(apply(ac_step, 1, fft, inverse = TRUE))) / nx * xstep * ystep

        di <- round(cp$dy / ystep)
        dj <- round(cp$dx / xstep)
        ridx <- ((-di) %% ny) + 1L
        cidx <- ((-dj) %% nx) + 1L

        overlap <- ac[cbind(ridx, cidx)]
        e_vals  <- ifelse(overlap > 0, areaW / overlap, 0)
      }

      denom_m <- m - 1
      denom_npts <- npts - 1

    } else {
      # border correction: only pairs whose *center* point i is far enough
      # from the boundary count, each with weight 1 so we basically
      # to zeroing out ineligible rows of the 0/1 adjacency matrix
      dist_to_boundary <- spatstat.geom::bdist.points(ppp_obj)
      eligible <- dist_to_boundary >= rvalue

      keep_elig <- eligible[i_idx]
      i_idx <- i_idx[keep_elig]
      j_idx <- j_idx[keep_elig]
      e_vals <- rep(1, length(i_idx))

      denom_m <- sum(eligible & marks_vec == mark1) # m_elig
      denom_npts <- sum(eligible)                  # n_elig
    }

    R0 <- sum(e_vals)
    R1 <- sum(e_vals^2)

    row_sums <- numeric(npts)
    if (length(i_idx) > 0) {
      agg <- tapply(e_vals, i_idx, sum)
      row_sums[as.integer(names(agg))] <- agg
    }
    R2 <- sum(row_sums^2) - R1
    R3 <- R0^2 - 2*R1 - 4*R2

    keep <- marks_vec[i_idx] == mark1 & marks_vec[j_idx] == mark1
    Ksum <- sum(e_vals[keep])

    f1 <- m*denom_m/npts/denom_npts
    f2 <- f1*(denom_m-1)/(denom_npts-1)
    f3 <- f2*(denom_m-2)/(denom_npts-2)

    K <- areaW * Ksum / m / denom_m
    mu_K <- areaW * R0 / npts / denom_npts
    var_K <- areaW^2 * (2*R1*f1 + 4*R2*f2 + R3*f3) / m / m / denom_m / denom_m - mu_K^2

    Z_k <- (K - mu_K) / sqrt(var_K)
    pval_appx <- pnorm(-Z_k)

    return(tibble(
      r = rvalue,
      k = K,
      theo_csr = pi * rvalue^2,
      kamp_csr = mu_K,
      kamp = K - mu_K,
      var = var_K,
      pvalue = min(1, pval_appx)
    ))

 } else if (correction %in% c("iso", "isotropic")) {
   # matrix-based: edge.Ripley needs the full n x n distance matrix, so this
   # path doesn't benefit from the closepairs shortcut used above
   pp_df <- as.data.frame(ppp_obj)
   W <- as.matrix(dist(as.matrix(select(pp_df, x, y))))

   e <- spatstat.explore::edge.Ripley(ppp_obj, r = W)
   W <- (W <= rvalue) * 1.0
   diag(W) <- 0
   Wr <- W * e

   R0 <- sum(Wr)
   R1 <- sum(Wr^2)
   R2 <- sum(rowSums(Wr)^2) - R1
   R3 <- R0^2 - 2*R1 - 4*R2

   Kmat <- Wr[which(marks_vec == mark1), which(marks_vec == mark1)]

   denom_m <- m - 1
   denom_npts <- npts - 1

   f1 <- m*denom_m/npts/denom_npts
   f2 <- f1*(denom_m-1)/(denom_npts-1)
   f3 <- f2*(denom_m-2)/(denom_npts-2)

   K <- areaW * sum(Kmat) / m / denom_m
   mu_K <- areaW * R0 / npts / denom_npts
   var_K <- areaW^2 * (2*R1*f1 + 4*R2*f2 + R3*f3) / m / m / denom_m / denom_m - mu_K^2

   Z_k <- (K - mu_K) / sqrt(var_K)
   pval_appx <- pnorm(-Z_k)

   return(tibble(
     r = rvalue,
     k = K,
     theo_csr = pi * rvalue^2,
     kamp_csr = mu_K,
     kamp = K - mu_K,
     var = var_K,
     pvalue = min(1, pval_appx)
   ))

 } else {
   stop("Only translational and isotropic edge corrections are implemented.")
 }
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
#' @param ... Additional arguments (currently unused).
#'
#' @returns The `ppp` point pattern object (built from `df` if it was a data.frame) if all
#' inputs are valid, otherwise throws an error with a descriptive message.
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
                         p_thin,...) {
  ppp_obj <- NULL
  # If it's already a ppp, use it and DO NOT run dataframe checks
  if (inherits(df, "ppp")) {
    ppp_obj <- df
    all_marks <- unique(spatstat.geom::marks(ppp_obj))

  } else if (is.data.frame(df)) {

    if (is.null(mark_var) || mark_var == "") {
      stop("mark_var must be supplied and cannot be NULL or empty.")
    }
    if (!all(c("x", "y") %in% colnames(df))) {
      stop("Input dataframe must contain 'x' and 'y' columns.")
    }
    if (!mark_var %in% colnames(df)) {
      stop(paste0("mark_var '", mark_var, "' not found in dataframe columns."))
    }

    df$mark_var <- as.factor(df[[mark_var]])
    if (nlevels(df$mark_var) < 2) {
      stop("The mark_var column must have at least two unique values.")
    }

    win <- spatstat.geom::convexhull.xy(df$x, df$y)
    ppp_obj <- spatstat.geom::ppp(df$x, df$y, window = win, marks = df$mark_var)
    all_marks <- unique(spatstat.geom::marks(ppp_obj))

  } else {
    stop("Input must be either a data.frame or a spatstat ppp object.")
  }

  message("We expect the dataframe to be a single point process. If you have multiple point processes, subset the dataframe by ID and please run KAMP separately for each process.")

  # Makes sure ppp_obj is not NULL
  if (is.null(ppp_obj)) {
    stop("The point pattern object cannot be NULL. Conversion of dataframe to ppp object failed.")
  }

  # Check if rvec is numeric
  if (!is.numeric(rvals) || any(rvals < 0)) {
    stop("rvals must be numeric and 0 or more")
  }

  # "translational"/"isotropic" are accepted as full-name aliases for "trans"/"iso"
  # and get normalized to the short form by kamp() before being dispatched onward
  if (correction %in% c("trans", "translational", "iso", "isotropic", "none") == FALSE) {
    stop("correction must be one of 'trans', 'translational', 'iso', 'isotropic', or 'none'.")
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

  if (sum(spatstat.geom::marks(ppp_obj) == mark1) < 5) {
    message(paste0("Less than 5 target cells marked as '", mark1, "'. This may lead to unreliable results."))
  }


  if (npoints(ppp_obj) > 100000) {
    message("Point pattern has more than 100,000 points. At this time, it is not recommended to use edge correction with KAMP on datasets of this size. Results may be unreliable. Consider using KAMP-lite (thin = TRUE) or subsetting your data.")
  }

  return(ppp_obj)
}
