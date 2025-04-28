#' Title KAMP Expectation Helper
#' @title get_kamp_expectation_mat_helper
#' @description Helper function to calculate the KAMP expectation for a given point pattern object and radius.
#' @param ppp_obj A point pattern object of class "ppp" from the spatstat package.
#' @param rvalue A single radius at which to calculate the KAMP expectation.
#' @param markvar The variable used to mark the points in the point pattern object. Default is "immune".
#'
#' @returns A tibble containing the radius, expectation, and time taken to compute the expectation.
#' @importFrom spatstat.explore Kcross Kest edge.Trans
#' @importFrom spatstat.geom area.owin ppp as.owin npoints Window
#' @importFrom dplyr mutate select
#' @importFrom tibble as_tibble
#' @importFrom magrittr %>%
#' @importFrom purrr map_dfr
#' @importFrom tictoc tic toc
#' @importFrom stats dist pnorm
#' @importFrom tibble tibble
#' @examples
get_kamp_expectation_mat_helper = function(rvalue,
                                           ppp_obj,
                                           markvar = "immune") {

  tic()
  npts = npoints(ppp_obj)
  W = Window(ppp_obj)
  areaW = spatstat.geom::area(W)

  # Compute translational edge correction weights, believe exact = FALSE by default
  e = edge.Trans(ppp_obj)  # Believe this is where the bottleneck is? might switch to border correction if too many cells

  pp_df = as.data.frame(ppp_obj)
  W = as.matrix(dist(as.matrix(select(pp_df, x, y))))
  W = ifelse(W <= rvalue, 1, 0)
  diag(W) = 0
  Wr = W * e  # Apply edge correction

  # Compute R0 term
  R0 = sum(Wr)

  # Compute expectation
  m = sum(ppp_obj$marks == markvar)
  mu_K = areaW * R0 / (m * (m - 1))

  time_k = toc()
  time <- time_k$toc - time_k$tic

  return(tibble(r = rvalue,
                expectation = mu_K,
                time = time))
}

#' Title
#'
#' @param ppp_obj A point pattern object of class "ppp" from the spatstat package.
#' @param rvalue A radius at which to calculate the KAMP variance.
#' @param correction Type of border correction
#'
#' @returns A tibble containing the radius, K statistic, expectation, variance, Z statistic, p-value, and time taken to compute the variance.
#' @importFrom spatstat.explore Kcross Kest edge.Trans
#' @importFrom spatstat.geom area.owin ppp as.owin npoints Window
#' @importFrom dplyr mutate select
#' @importFrom tibble as_tibble
#' @importFrom magrittr %>%
#' @importFrom purrr map_dfr
#' @importFrom tictoc tic toc
#' @importFrom stats dist pnorm
#' @importFrom tibble tibble
#' @examples
get_kamp_variance_helper = function(rvalue,
                                    ppp_obj,
                                    correction = "trans") {
  npts = npoints(ppp_obj)
  W = Window(ppp_obj)
  areaW = spatstat.geom::area(W)

  # Compute translational edge correction weights, believe exact = FALSE by default
  e = edge.Trans(ppp_obj)  # Believe this is where the bottleneck is?

  pp_df = as.data.frame(ppp_obj)
  W = as.matrix(dist(as.matrix(select(pp_df, x, y))))
  W = ifelse(W <= rvalue, 1, 0)
  diag(W) = 0
  Wr = W * e  # Apply edge correction

  # Compute R0 term
  R0 = sum(Wr)

  # Compute expectation
  m = sum(ppp_obj$marks == "immune")
  R1 = sum(Wr^2)
  R2 = sum(rowSums(Wr)^2) - R1
  R3 = R0^2 - 2*R1 - 4*R2

  npairs =  npts * (npts - 1)
  f1 = m*(m-1)/npts/(npts-1)
  f2 = f1*(m-2)/(npts-2)
  f3 = f2*(m-3)/(npts-3)

  Kmat = Wr[which(ppp_obj$marks == "immune"),which(ppp_obj$marks == "immune")]
  K = areaW * sum(Kmat)/ m / (m - 1) # Ripley's K based on translation correction
  mu_K = areaW * R0 / npts / (npts - 1) # expectation
  var_K = areaW^2 * (2 * R1 * f1 + 4 * R2 * f2 + R3 * f3) / m / m / (m - 1) / (m - 1) - mu_K^2   # variance


  Z_k = (K - mu_K) / sqrt(var_K) # Test statistic
  pval_appx = pnorm(-Z_k) # approximated p-value based on normal distribution


  result = tibble(
    r = rvalue,
    khat = K,
    expectation = mu_K, # K expectation under permutation distributions
    var = var_K,
    Z = Z_k, # test statistic
    pvalue = min(1, pval_appx)
  )

  return(result)
}

