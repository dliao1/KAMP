#' KAMP Variance
#' @title KAMP univariate variance - Rcpp Version
#'
#' @description
#' Computes the KAMP (K-function Adjusted for Marked Permutations) variance
#' for a given spatial point pattern. Also returns the KAMP expectation,
#' z-statistic, and p-value. Calculations are performed across a vector of
#' radii in a single call.
#'
#' Note: this a matrix-based implementation of the KAMP variance. It relies on
#' Rcpp functions (`kamp_pair_sums_trans_iso`, `kamp_pair_sums_bord`) to do the
#' pairwise sums, and `spatstat` for point pattern distances/edge corrections.
#' Translational, isotropic, border, and no ("none") edge corrections are all supported.
#'
#'
#' @param ppp_obj A point pattern object of class "ppp" from the spatstat package.
#' @param rvals A vector of radii at which to calculate the KAMP variance.
#' @param correction Type of edge correction. Defaults to translational.
#' @param mark1 The variable used to mark the points in the point pattern object. Default is "immune".
#'
#' @returns
#' A dataframe with the following columns:
#' \describe{
#'   \item{r}{The radius at which K was calculated.}
#'   \item{k}{The observed K value}
#'   \item{theo_csr}{The theoretical K under CSR}
#'   \item{kamp_csr}{The adjusted CSR representing the KAMP permuted expectation.}
#'   \item{kamp}{The difference between observed K and KAMP CSR}
#'   \item{var}{Variance of K under the permutation null distribution}
#'   \item{pval}{P-value, calculated using the formula: pnorm(-z)}
#' }
#'
#' @export
#' @importFrom purrr map_dfr
#' @importFrom dplyr %>%
#' @importFrom tibble as_tibble
#' @useDynLib KAMP, .registration = TRUE
#' @importFrom Rcpp sourceCpp evalCpp
#' @examples
#' win <- spatstat.geom::owin(c(0, 1), c(0, 1))
#' pp <- spatstat.random::rpoispp(lambda = 150, win = win)
#' marks <- sample(c("immune", "background"), pp$n, replace = TRUE, prob = c(0.4, 0.6))
#' marked_pp <- spatstat.geom::ppp(pp$x, pp$y, window = win, marks = factor(marks))
#'
#' result <- kamp_variance_Rcpp(marked_pp, rvals = c(0.05, 0.1, 0.15), mark1 = "immune")
#' print(result)
kamp_variance_Rcpp = function(ppp_obj,
                              rvals = c(0, .05, .075, .1, .15, .2),
                              correction = "trans",
                              mark1 = "immune") {
  npts  <- npoints(ppp_obj)
  Wwin  <- Window(ppp_obj)
  areaW <- spatstat.geom::area(Wwin)

  marks_vec <- spatstat.geom::marks(ppp_obj) == mark1
  m <- sum(marks_vec) # number of marked points

  rmax <- max(rvals) # maximum radius
  cp <- spatstat.geom::closepairs(ppp_obj, rmax = rmax, what = "all") # get all close pairs up to rmax

  # Sanity check for indices
  if (any(cp$i < 1 | cp$i > npts)) {
    stop("Invalid i indices from closepairs")
  }
  if (any(cp$j < 1 | cp$j > npts)) {
    stop("Invalid j indices from closepairs")
  }


  bdist <- NULL
  if (correction == "trans") {
    # paired = TRUE so this only weights the pairs we actually have,
    # instead of building the whole n x n matrix and subsetting it
    w <- edge.Trans(ppp_obj[cp$i], ppp_obj[cp$j], paired = TRUE)

    if (any(is.na(w)) || any(!is.finite(w))) {
      stop("Invalid edge correction weights")
    }
  } else if (correction == "border") {
    bdist <- spatstat.geom::bdist.points(ppp_obj)
    w <- rep(1, length(cp$i))

  } else if (correction == "iso") {
    # same idea -- per-pair, not the full n x n distance/weight matrices
    bdist_all <- spatstat.geom::bdist.points(ppp_obj)
    w <- as.numeric(spatstat.explore::edge.Ripley(ppp_obj[cp$i],
                                                   matrix(cp$d, ncol = 1),
                                                   bdistX = bdist_all[cp$i]))

    if (any(is.na(w)) || any(!is.finite(w))) {
      stop("Invalid isotropic edge correction weights")
    }

  } else if (correction == "none") {
    # no edge correction at all: every close pair counts with weight 1
    w <- rep(1, length(cp$i))

  } else {
    stop("correction must be one of 'trans', 'iso', 'border', or 'none'.")
  }

  ord <- order(cp$d) # orders close pairs by distance
  d_sorted <- cp$d[ord] # sorted distances
  i_sorted <- cp$i[ord] # rearranges i and j indices according to sorted distances
  j_sorted <- cp$j[ord]
  w_sorted <- w[ord] # sorted weights


  # Returns if no close pairs
  if (length(d_sorted) == 0) {
    return(tibble::tibble(
      r = rvals,
      k = 0,
      theo_csr = pi * rvals^2,
      kamp_csr = 0,
      kamp_diff = 0,
      var = 0,
      pvalue = 1
    ))
  }

  # Calls Rcpp helper for sums -- trans/iso/none all reduce to the same
  # generic weighted pair-sum kernel, differing only in how w was built above
  rcpp_sums <- if (correction == "border") {
    kamp_pair_sums_bord(i = i_sorted,
                         j = j_sorted,
                         d = d_sorted,
                         rvals = rvals,
                         is_mark1 = marks_vec,
                         bdist = bdist,
                         npts = npts)
  } else {
    kamp_pair_sums_trans_iso(i = i_sorted,
                              j = j_sorted,
                              d = d_sorted,
                              w = w_sorted,
                              rvals = rvals,
                              is_mark1 = marks_vec,
                              npts = npts)
  }

  R0 <- rcpp_sums$R0
  R1 <- rcpp_sums$R1
  R2 <- rcpp_sums$R2
  Ksum <- rcpp_sums$Ksum

  R3 <- R0^2 - 2 * R1 - 4 * R2

  if (correction == "border") {
    # border correction restricts eligible reference points per radius, so
    # the denominators depend on r via n_elig/m_elig instead of npts/m
    n_elig <- vapply(rvals, function(rv) sum(bdist >= rv), numeric(1))
    m_elig <- vapply(rvals, function(rv) sum(bdist >= rv & marks_vec), numeric(1))

    f1 <- m_elig * m / (n_elig * npts)
    f2 <- f1 * (m_elig - 1) / (n_elig - 1)
    f3 <- f2 * (m_elig - 2) / (n_elig - 2)

    K <- (areaW * Ksum) / (m * m_elig)
    mu_K <- (areaW * R0) / (n_elig * npts)
    var_K <- (areaW^2 * (2 * R1 * f1 + 4 * R2 * f2 + R3 * f3) /
                (m^2 * m_elig^2)) - mu_K^2
  } else {
    f1 <- m * (m - 1) / (npts * (npts - 1))
    f2 <- f1 * (m - 2) / (npts - 2)
    f3 <- f2 * (m - 3) / (npts - 3)

    K <- (areaW * Ksum) / (m * (m - 1))
    mu_K <- (areaW * R0) / (npts * (npts - 1))
    var_K <- (areaW^2 * (2 * R1 * f1 + 4 * R2 * f2 + R3 * f3) /
                (m^2 * (m - 1)^2)) - mu_K^2
  }

  var_K <- pmax(var_K, 0)

  Z_k = (K - mu_K) / sqrt(var_K)
  pval_appx = pnorm(-Z_k)


  tibble::tibble(
    r = rvals,
    k = K,
    theo_csr = pi * rvals^2,
    kamp_csr = mu_K,
    kamp = k - kamp_csr,
    var = var_K,
    pvalue = pmin(1, pval_appx)
  )
}


