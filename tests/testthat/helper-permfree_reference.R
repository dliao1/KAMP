# Reference implementation of the exact permutation variance of univariate G,
# copied from PermFree-NearestNeighborG (functions.R, exact_G_variance)
# loops over every pair of points in R, so it is slow; only used to check gamp_sums()
# neighbors are D < r here, while gamp uses d <= r (same unless a distance equals r exactly)

log_comb <- function(n, k) {
  lgamma(n + 1) - (lgamma(k + 1) + lgamma(n - k + 1))
}

binom_prob <- function(n, k, x) {
  exp(log_comb(n - x, k) - log_comb(n, k))
}

exact_G_variance <- function(df, mark_col, mark, radii) {
  N  = nrow(df)
  n  = sum(df[[mark_col]] == mark)
  if (n < 2) stop("Need at least two marked points.")
  ell = n - 2
  D   = as.matrix(dist(df[, c("x", "y")]))

  res = numeric(length(radii))

  for (i in seq_along(radii)) {
    r = radii[i]
    neigh  = (D < r) & (D > 0)
    n_p    = rowSums(neigh)
    log_num = log_comb(N - 1 - n_p, n - 1)
    log_den = log_comb(N - 1, n - 1)
    mu_p    = (n / N) * (1 - exp(log_num - log_den))

    logC_Nm2_ell = log_comb(N - 2, ell)
    logC_N_n     = log_comb(N,     n)
    log_common   = logC_Nm2_ell - logC_N_n
    common_fac   = exp(log_common)

    pi_mat = matrix(0, N, N)
    for (p in 1:(N - 1)) {
      Np_all = which(neigh[p, ])
      for (q in (p + 1):N) {
        if (neigh[p, q]) {
          pi = common_fac
        } else {
          Nq_all  = which(neigh[q, ])
          both    = intersect(Np_all, Nq_all)
          only_p  = setdiff(Np_all, c(q, both))
          only_q  = setdiff(Nq_all, c(p, both))

          a = length(only_p)
          b = length(only_q)
          w = length(both)

          fail_p    = binom_prob(N - 2, ell, a + w)
          fail_q    = binom_prob(N - 2, ell, b + w)
          fail_both = binom_prob(N - 2, ell, a + b + w)

          prob_valid = 1 - fail_p - fail_q + fail_both

          pi = common_fac * prob_valid
        }
        pi_mat[p, q] = pi_mat[q, p] = pi
      }
    }
    var_G = (1 / n^2) * (
      sum(mu_p * (1 - mu_p)) +
        2 * sum(pi_mat[upper.tri(pi_mat)] -
                  outer(mu_p, mu_p, "*")[upper.tri(pi_mat)])
    )
    res[i] = pmax(var_G, 0)
  }
  res
}
