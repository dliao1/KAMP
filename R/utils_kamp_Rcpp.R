#' Close pairs and edge correction weights for the Rcpp KAMP variance
#'
#' @description Finds all close pairs up to `max(rvals)`, builds their edge
#'   correction weights, and sorts them by distance, as expected by the
#'   `kamp_pair_sums_*` Rcpp functions.
#'
#' @param ppp_obj A point pattern object of class "ppp".
#' @param rvals A vector of radii.
#' @param correction One of "trans", "iso", "border", or "none". Border pairs
#'   get weight 1; eligibility is handled in `kamp_pair_sums_bord`.
#'
#' @returns A list with sorted pair indices `i`, `j`, distances `d`, and weights `w`.
#' @keywords internal
kamp_pair_prep <- function(ppp_obj, rvals, correction) {
  npts <- npoints(ppp_obj)
  cp <- spatstat.geom::closepairs(ppp_obj, rmax = max(rvals), what = "all") # get all close pairs up to rmax

  # Sanity check for indices
  if (any(cp$i < 1 | cp$i > npts) || any(cp$j < 1 | cp$j > npts)) {
    stop("Invalid indices from closepairs")
  }

  if (length(cp$i) == 0L) {
    return(list(i = integer(0), j = integer(0), d = numeric(0), w = numeric(0)))
  }

  if (correction == "trans") {
    # paired = TRUE so this only weights the pairs we actually have,
    # instead of building the whole n x n matrix and subsetting it
    w <- spatstat.explore::edge.Trans(ppp_obj[cp$i], ppp_obj[cp$j], paired = TRUE)

    if (any(is.na(w)) || any(!is.finite(w))) {
      stop("Invalid edge correction weights")
    }
  } else if (correction == "iso") {
    # same idea -- per-pair, not the full n x n distance/weight matrices
    bdist_all <- spatstat.geom::bdist.points(ppp_obj)
    w <- as.numeric(spatstat.explore::edge.Ripley(ppp_obj[cp$i],
                                                   matrix(cp$d, ncol = 1),
                                                   bdistX = bdist_all[cp$i]))

    if (any(is.na(w)) || any(!is.finite(w))) {
      stop("Invalid isotropic edge correction weights")
    }
  } else {
    # border and none: every close pair counts with weight 1
    w <- rep(1, length(cp$i))
  }

  ord <- order(cp$d) # orders close pairs by distance
  list(i = cp$i[ord], j = cp$j[ord], d = cp$d[ord], w = w[ord])
}

#' KAMP variance result when no pairs fall within `max(rvals)`
#' @noRd
kamp_empty_result <- function(rvals) {
  tibble::tibble(
    r = rvals,
    k = 0,
    theo_csr = pi * rvals^2,
    kamp_csr = 0,
    kamp = 0,
    var = 0,
    pvalue = 1
  )
}
