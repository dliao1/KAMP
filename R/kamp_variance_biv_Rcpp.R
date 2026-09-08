#' Computes KAMP Variance for Bivariate Point Patterns
#'
#' @title KAMP bivariate variance
#' @description
#' Computes the KAMP (K-function Adjusted for Marked Permutations) variance
#' for bivariate point patterns. This function calculates Ripley's K
#' using both the traditional Ripley's K method and the KAMP-adjusted CSR
#' baseline by using a matrix-based implementation.
#'
#' The KAMP-adjusted CSR represents a more realistic baseline for K (compared
#' to traditional CSR) that accounts for spatial clustering or inhomogeneity
#' in a point pattern compared to the traditional CSR assumption, while
#' avoiding the computational burden of permuting the point pattern.
#'
#' Note: This function implements a slower, matrix-based implementation
#' of the KAMP variance. It is a wrapper around the `kamp_variance_biv_helper`
#' function that calculates the KAMP variance at one radius and maps it
#' over a vector of radii.
#'
#' @param ppp_obj A point pattern object from the `spatstat.geom` package.
#' @param rvals A vector of radii at which to calculate the KAMP expectation. Defaults to c(0, 0.05, 0.075, 0.1, 0.15, 0.2).
#' @param correction Type of edge correction. Defaults to translational.
#' @param mark1 Variable used to mark the points in the point pattern object for the first type. Default is "immune1".
#' @param mark2 Variable used to mark the points in the point pattern object for the second type. Default is "immune2".
#'
#' @importFrom spatstat.explore Kcross Kest edge.Trans edge.Ripley
#' @importFrom spatstat.geom area.owin ppp as.owin npoints Window
#' @importFrom dplyr mutate select
#' @importFrom tibble as_tibble
#' @importFrom magrittr %>%
#' @importFrom purrr map_dfr
#' @importFrom stats dist pnorm
#' @importFrom tibble tibble
#'
#' @returns
#' A dataframe with the following columns:
#' \describe{
#' \item{r}{The radius at which K was calculated.}
#' \item{k}{The observed K value}
#' \item{theo_csr}{The theoretical K under CSR}
#' \item{kamp_csr}{The adjusted CSR representing the KAMP permuted expectation.}
#' \item{var}{Variance of K under the permutation null distribution}
#' \item{z}{Z statistic, calculated by normalizing K using the formula: (K - KAMP)/sqrt(var)}
#' \item{pval}{P-value, calculated using the formula: pnorm(-z)}
#' }
#'
#' @export
#' @examples
#' win <- spatstat.geom::owin(c(0, 1), c(0, 1))
#' pp <- spatstat.random::rpoispp(lambda = 150, win = win)
#' mark_labels <- c("immune1", "immune2", "background")
#' marks <- sample(mark_labels, pp$n, replace = TRUE, prob = c(0.3, 0.3, 0.4))
#' marked_pp <- spatstat.geom::ppp(pp$x, pp$y, window = win, marks = factor(marks))
#'
#' result <- kamp_variance_biv_Rcpp(marked_pp, rvals = c(0.05, 0.1),
#'                                  mark1 = "immune1", mark2 = "immune2")
#' print(result)
kamp_variance_biv_Rcpp <- function(ppp_obj,
                                   rvals = c(0, .05, .075, .1, .15, .2),
                                   correction = "trans",
                                   mark1 = "immune1",
                                   mark2 = "immune2") {

  npts = npoints(ppp_obj)
  ppp_window = Window(ppp_obj)
  areaW = spatstat.geom::area(ppp_window)

  m1 = sum(ppp_obj$marks == mark1)
  m2 = sum(ppp_obj$marks == mark2)


  is_m1 <- spatstat.geom::marks(ppp_obj) == mark1
  is_m2 <- spatstat.geom::marks(ppp_obj) == mark2

  rmax <- max(rvals)
  cp <- spatstat.geom::closepairs(ppp_obj, rmax = rmax, what = "all")
  w <- NULL

  # no pairs within rmax
  if (length(cp$i) == 0L) {
    R0 <- R1 <- R2 <- Ksum12 <- rep(0, length(rvals))
  } else {

    if (correction == "trans") {
      # paired = TRUE so this only weights the pairs we actually have,
      # instead of building the whole n x n matrix and subsetting it
      w <- spatstat.explore::edge.Trans(ppp_obj[cp$i], ppp_obj[cp$j], paired = TRUE)
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
      stop("correction must be one of 'trans', 'iso', or 'none'.")
    }

  }

  ord <- order(cp$d) # orders close pairs by distance
  d_sorted <- cp$d[ord] # sorted distances
  i_sorted <- cp$i[ord] # rearranges i and j indices according to sorted distances
  j_sorted <- cp$j[ord]
  w_sorted <- w[ord] # sorted weights

  sums <- kamp_pair_sums_trans_iso_biv(
    i = i_sorted,
    j = j_sorted,
    d = d_sorted,
    w = w_sorted,
    rvals = rvals,
    is_mark1 = is_m1,
    is_mark2 = is_m2,
    npts = npts
  )

  R0 <- sums$R0
  R1 <- sums$R1
  R2 <- sums$R2
  Ksum <- sums$Ksum


  R3 <- R0^2 - 2 * R1 - 4 * R2

  f1 <- m1 * m2 / (npts * (npts - 1))
  f2 <- f1 * (m1 + m2 - 2) / (npts - 2)
  f3 <- f1 * (m1 - 1) * (m2 - 1) / ((npts - 2) * (npts - 3))

  K    <- areaW * Ksum / (m1 * m2)
  mu_K <- areaW * R0 / (npts * (npts - 1))
  var_K <- areaW^2 * (R1 * f1 + R2 * f2 + R3 * f3) / (m1^2 * m2^2) - mu_K^2
  var_K <- pmax(var_K, 0)

  Z <- (K - mu_K) / sqrt(var_K)
  p <- pnorm(-Z)

  tibble::tibble(
    r = rvals,
    k = K,
    theo_csr = pi * rvals^2,
    kamp_csr = mu_K,
    kamp = k - kamp_csr,
    var = var_K,
    pvalue = pmin(1, p)
  )
}
