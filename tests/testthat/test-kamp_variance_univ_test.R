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

test_that("kamp_variance_Rcpp matches kamp_variance_helper for border correction", {
  marked_pp <- make_univ_pp(seed = 9, lambda = 200)
  rvals <- c(0.02, 0.05, 0.1)

  rcpp_result <- kamp_variance_Rcpp(marked_pp, rvals = rvals, correction = "border")
  helper_result <- purrr::map_dfr(
    rvals,
    ~kamp_variance_helper(marked_pp, rvalue = .x, correction = "border")
  )

  expect_equal(rcpp_result$k, helper_result$k, tolerance = 1e-8)
  expect_equal(rcpp_result$kamp_csr, helper_result$kamp_csr, tolerance = 1e-8)
})

test_that("kamp with variance = TRUE validates edge correction argument", {
  set.seed(8)
  df <- data.frame(x = runif(50), y = runif(50),
                    immune = sample(c("immune", "background"), 50, replace = TRUE))

  expect_error(
    kamp(df, rvals = c(0.05, 0.1), mark_var = "immune", mark1 = "immune",
         variance = TRUE, correction = "border"),
    "correction must be one of 'trans', 'translational', 'iso', 'isotropic', or 'none'"
  )
})

test_that("kamp accepts full-name correction aliases and 'none' for variance", {
  marked_pp <- make_univ_pp()
  rvals <- c(0.05, 0.1)

  trans_result <- kamp(marked_pp, rvals = rvals, mark1 = "immune",
                        variance = TRUE, correction = "trans")
  translational_result <- kamp(marked_pp, rvals = rvals, mark1 = "immune",
                                variance = TRUE, correction = "translational")
  expect_equal(trans_result, translational_result)

  iso_result <- kamp(marked_pp, rvals = rvals, mark1 = "immune",
                      variance = TRUE, correction = "iso")
  isotropic_result <- kamp(marked_pp, rvals = rvals, mark1 = "immune",
                            variance = TRUE, correction = "isotropic")
  expect_equal(iso_result, isotropic_result)

  none_result <- kamp(marked_pp, rvals = rvals, mark1 = "immune",
                       variance = TRUE, correction = "none")
  expect_equal(names(none_result), c("r", "k", "theo_csr", "kamp_csr", "kamp", "var", "pvalue"))
  expect_true(all(none_result$pvalue >= 0 & none_result$pvalue <= 1))
})
