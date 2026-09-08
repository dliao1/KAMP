#' Helper function for bivariate KAMP Variance
#'
#' @title kamp_variance_biv_helper
#' @description
#' Helper function to calculate the KAMP variance for bivariate point patterns.
#'
#' @param ppp_obj A point pattern object "ppp" from the `spatstat` package.
#' @param rval A single radius
#' @param correction Type of edge correction. Defaults to translational.
#' @param mark1 Variable used to mark the points in the point pattern object for the first type. Default is "immune1".
#' @param mark2 Variable used to mark the points in the point pattern object for the second type. Default is "immune2".
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
#' @importFrom stats dist pnorm fft
#' @importFrom tibble tibble
#'
#' @keywords internal
#'
kamp_variance_biv_helper <- function(ppp_obj,
                                     rval,
                                     correction = "trans",
                                     mark1 = "immune1",
                                     mark2 = "immune2") {

  npts = npoints(ppp_obj)
  ppp_window = Window(ppp_obj)
  areaW = spatstat.geom::area(ppp_window)

  if (correction %in% c("trans", "translational", "none")) {
    # closepairs-based -- only touches pairs within rval, never builds
    # the full n x n distance/weight matrix
    cp <- spatstat.geom::closepairs(ppp_obj, rmax = rval, what = "all", twice = TRUE)
    i_idx <- cp$i
    j_idx <- cp$j

    if (correction == "none") {
      # no edge correction at all: every close pair counts with weight 1
      e_vals <- rep(1, length(i_idx))

    } else if (ppp_window$type == "rectangle") {
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

    R0 <- sum(e_vals)
    R1 <- sum(e_vals^2)

    row_sums <- numeric(npts)
    if (length(i_idx) > 0) {
      agg <- tapply(e_vals, i_idx, sum)
      row_sums[as.integer(names(agg))] <- agg
    }
    R2 <- sum(row_sums^2) - R1
    R3 <- R0^2 - 2*R1 - 4*R2

    marks_vec <- as.character(ppp_obj$marks)
    m1 <- sum(marks_vec == mark1)
    m2 <- sum(marks_vec == mark2)
    keep <- marks_vec[i_idx] == mark1 & marks_vec[j_idx] == mark2
    Ksum <- sum(e_vals[keep])

    f1 <- m1*m2/npts/(npts-1)
    f2 <- f1*(m1+m2-2)/(npts-2)
    f3 <- f1*(m1-1)*(m2-1)/(npts-2)/(npts-3)

    K <- areaW * Ksum / m1 / m2
    mu_K <- areaW * R0 / npts / (npts-1)
    var_K <- areaW^2*(R1*f1 + R2*f2 + R3*f3)/m1/m1/m2/m2 - mu_K^2

    Z_k <- (K-mu_K) / sqrt(var_K)
    pval_appx <- pnorm(-Z_k)

    return(tibble(
      r = rval,
      k = K,
      theo_csr = pi * rval^2,
      kamp_csr = mu_K,
      kamp = K - mu_K,
      var = var_K,
      pvalue = min(1, pval_appx)
    ))

  } else if (correction %in% c("iso", "isotropic")) {
    pp_df = as.data.frame(ppp_obj)
    W = as.matrix(dist(as.matrix(select(pp_df, x, y))))

    e = edge.Ripley(ppp_obj, r = W)
    W = ifelse(W <= rval, 1, 0)
    diag(W) = 0
    Wr = W * e

    R0 = sum(Wr)
    R1 = sum(Wr^2)
    R2 = sum(rowSums(Wr)^2) - R1
    R3 = R0^2 - 2*R1 - 4*R2

    m1 = sum(ppp_obj$marks == mark1)
    m2 = sum(ppp_obj$marks == mark2)
    f1 = m1*m2/npts/(npts-1)
    f2 = f1*(m1+m2-2)/(npts-2)
    f3 = f1*(m1-1)*(m2-1)/(npts-2)/(npts-3)

    Kmat = Wr[which(ppp_obj$marks == mark1),which(ppp_obj$marks == mark2)]
    K = areaW * sum(Kmat) / m1 / m2 # Ripley's K based on translation correction
    mu_K = areaW * R0/npts / (npts-1) # expectation
    var_K = areaW^2*(R1*f1 + R2*f2 + R3*f3)/m1/m1/m2/m2 - mu_K^2

    Z_k = (K-mu_K) / sqrt(var_K) # Test statistic
    pval_appx = pnorm(-Z_k) # approximated p-value based on normal distribution

    return(tibble(
      r = rval,
      k = K,
      theo_csr = pi * rval^2, # theoretical CSR using area of circle
      kamp_csr = mu_K, # K expectation under permutation distributions
      kamp = K - mu_K,
      var = var_K,
      pvalue = min(1, pval_appx)
    ))

  } else {
    stop("Only translational, and isotropic edge correction are currently supported.")
  }

}




