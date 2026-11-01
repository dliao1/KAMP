#' GAMP expectation helper
#'
#' @description Exact permutation expectation of G, as in PermFree-NearestNeighborG:
#' `E(G) = 1 - (1/N) * sum_p choose(N - 1 - n_p, t) / choose(N - 1, t)`, where n_p(r)
#' is the number of neighbors of point p within r and t is the number of mark2 points
#' among the other N - 1 points (m1 - 1 univariate, m2 bivariate).
#'
#' @param ppp_obj A point pattern object "ppp" from the `spatstat` package.
#' @param rvals Vector of radii.
#' @param mark1 Mark of the points G is computed for (e.g. immune).
#' @param mark2 Mark of the nearest-neighbor points for bivariate G, or `NULL` for univariate G.
#' @param correction Type of edge correction. Only "none" is supported.
#'
#' @returns A dataframe with columns `r`, `g`, `theo_csr`, `gamp_csr`, `gamp`.
#' @keywords internal
gamp_expectation_helper <- function(ppp_obj, rvals, mark1, mark2 = NULL, correction = "none") {
  if (correction != "none") {
    stop("GAMP currently only supports correction = 'none'.")
  }

  npts <- npoints(ppp_obj)
  areaW <- spatstat.geom::area(Window(ppp_obj))

  univariate <- is.null(mark2)
  marks_vec <- spatstat.geom::marks(ppp_obj)
  is1 <- marks_vec == mark1 # vector for TRUE/FALSE, TRUE if point is of type mark1 and FALSE otherwise
  m1 <- sum(is1) # number of points of type mark1

  # m_nn = number of points G looks for as nearest neighbors (mark1 if univariate)
  # t1 = how many of those are among the other N - 1 points, given p is mark1
  if (univariate == TRUE) {
    m_nn <- m1
    t1 <- m1 - 1
  } else {
    m_nn <- sum(marks_vec == mark2)
    t1 <- m_nn
  }

  cp <- spatstat.geom::closepairs(ppp_obj, rmax = max(rvals), what = "ijd", twice = TRUE)

  # number of neighbors of each point within each radius

  # for each radius curr_rval, keep only the pairs with distance <= curr_rval, and
  # then count how many neighbors each point has within that radius

  # Ex: One radius: tabulate(cp$i[cp$d <= curr_rval], nbins = npts)

  # For curr_rval = 0.5:
  #  1. cp$d <= 0.5 -> TRUE TRUE FALSE FALSE FALSE FALSE
  #2. cp$i[...] -> c(1, 2), the points that have a neighbor within 0.5
  #3. tabulate(c(1, 2), nbins = 4) → 1 1 0 0

  counts <- vapply(rvals, function(curr_rval) tabulate(cp$i[cp$d <= curr_rval], nbins = npts), integer(npts))
  counts <- matrix(counts, nrow = npts)

  # P(none of a point's neighbors is marked) = choose(N - 1 - n_p, t1) / choose(N - 1, t1)
  binom_vals <- exp(lchoose(npts - 1 - counts, t1) - lchoose(npts - 1, t1))
  exact_csr_values <- 1 - (1 / npts) * colSums(binom_vals) # fraction of marked points whose nearest other marked point is within r


  rvals_sorted <- sort(unique(c(0, rvals)))
  g_fun <- if (univariate == TRUE) {
    spatstat.explore::Gest(ppp_obj[is1], r = rvals_sorted, correction = "none")
  } else {
    spatstat.explore::Gcross(ppp_obj, i = mark1, j = mark2, r = rvals_sorted, correction = "none")
  }

  tibble::tibble(
    r = rvals,
    g = g_fun$raw[match(rvals, rvals_sorted)],
    theo_csr = 1 - exp(-(m_nn / areaW) * pi * rvals^2),
    gamp_csr = exact_csr_values,
    gamp = g - gamp_csr
  )
}


#' GAMP variance helper
#'
#' @description Exact permutation expectation and variance of G, with the pair
#' sums done in Rcpp (`gamp_sums`) from `spatstat.geom::closepairs()`.
#'
#' @inheritParams gamp_expectation_helper
#'
#' @returns A dataframe with columns `r`, `g`, `theo_csr`, `gamp_csr`, `gamp`, `var`, `pvalue`.
#' @keywords internal
gamp_variance_helper <- function(ppp_obj, rvals, mark1, mark2 = NULL, correction = "none") {
  if (correction != "none") {
    stop("GAMP currently only supports correction = 'none'.")
  }

  npts <- npoints(ppp_obj)
  areaW <- spatstat.geom::area(Window(ppp_obj))

  univariate <- is.null(mark2)
  marks_vec <- spatstat.geom::marks(ppp_obj)
  is1 <- marks_vec == mark1
  m1 <- sum(is1)

  # m_nn = number of points G looks for as nearest neighbors (mark1 if univariate);
  # gamp_sums takes m2 = 0 for univariate
  if (univariate == TRUE) {
    m_nn <- m1
    m2 <- 0L
  } else {
    m_nn <- sum(marks_vec == mark2)
    m2 <- m_nn
  }

  cp <- spatstat.geom::closepairs(ppp_obj, rmax = max(rvals), what = "ijd", twice = TRUE)

  ord <- order(cp$i, cp$d)
  start <- c(0L, cumsum(tabulate(cp$i, nbins = npts)))

  sums <- gamp_sums(nbr = cp$j[ord],
                    d = cp$d[ord],
                    start = as.integer(start),
                    rvals = rvals,
                    npts = npts,
                    m1 = m1,
                    m2 = m2,
                    univariate = univariate)


  rvals_sorted <- sort(unique(c(0, rvals)))
  g_fun <- if (univariate == TRUE) {
    spatstat.explore::Gest(ppp_obj[is1], r = rvals_sorted, correction = "none")
  } else {
    spatstat.explore::Gcross(ppp_obj, i = mark1, j = mark2, r = rvals_sorted, correction = "none")
  }
  g_obs <- g_fun$raw[match(rvals, rvals_sorted)]
  var_G <- sums$var_G
  Z_g <- (g_obs - sums$mu_G) / sqrt(var_G)

  tibble::tibble(
    r = rvals,
    g = g_obs,
    theo_csr = 1 - exp(-(m_nn / areaW) * pi * rvals^2),
    gamp_csr = sums$mu_G,
    gamp = g - gamp_csr,
    var = var_G,
    pvalue = ifelse(var_G > 0, pnorm(-Z_g), 1)
  )
}
