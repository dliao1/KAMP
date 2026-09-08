make_univ_pp <- function(seed = 1, lambda = 150, p_immune = 0.4) {
  set.seed(seed)
  win <- spatstat.geom::owin(c(0, 1), c(0, 1))
  pp <- spatstat.random::rpoispp(lambda = lambda, win = win)
  marks <- sample(c("immune", "background"), pp$n, replace = TRUE,
                   prob = c(p_immune, 1 - p_immune))
  spatstat.geom::ppp(pp$x, pp$y, window = win, marks = factor(marks))
}

test_that("kamp_variance returns expected columns and defaults to translational correction", {
  marked_pp <- make_univ_pp()
  rvals <- c(0.05, 0.1)

  result <- kamp_variance(marked_pp, rvals = rvals)
  expected_default <- kamp_variance(marked_pp, rvals = rvals, correction = "trans")

  expect_equal(names(result), c("r", "k", "theo_csr", "kamp_csr", "kamp", "var", "pvalue"))
  expect_equal(nrow(result), length(rvals))
  expect_true(is.numeric(result$var))
  expect_true(all(result$pvalue >= 0 & result$pvalue <= 1))
  expect_equal(result, expected_default)
})

test_that("kamp_variance respects a custom rvals vector", {
  marked_pp <- make_univ_pp()
  rvals <- c(0.02, 0.05, 0.08)

  result <- kamp_variance(marked_pp, rvals = rvals)
  expect_equal(result$r, rvals)
})

test_that("kamp_variance_Rcpp matches kamp_variance for translational correction", {
  marked_pp <- make_univ_pp(seed = 7)
  rvals <- c(0.05, 0.1, 0.15)

  base_result <- kamp_variance(marked_pp, rvals = rvals)
  rcpp_result <- kamp_variance_Rcpp(marked_pp, rvals = rvals)

  expect_equal(names(rcpp_result), names(base_result))
  expect_equal(rcpp_result$k, base_result$k, tolerance = 1e-8)
  expect_equal(rcpp_result$kamp_csr, base_result$kamp_csr, tolerance = 1e-8)
  expect_equal(rcpp_result$var, base_result$var, tolerance = 1e-8)
})

test_that("kamp with variance = TRUE validates edge correction argument", {
  set.seed(8)
  df <- data.frame(x = runif(50), y = runif(50),
                    immune = sample(c("immune", "background"), 50, replace = TRUE))

  expect_error(
    kamp(df, rvals = c(0.05, 0.1), mark_var = "immune", mark1 = "immune",
         variance = TRUE, correction = "border"),
    "Currently only isotropic and translational edge correction are supported"
  )
})
