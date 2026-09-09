make_biv_pp <- function(seed = 1, lambda = 150) {
  set.seed(seed)
  win <- spatstat.geom::owin(c(0, 1), c(0, 1))
  pp <- spatstat.random::rpoispp(lambda = lambda, win = win)
  marks <- sample(c("immune1", "immune2", "background"), pp$n, replace = TRUE)
  spatstat.geom::ppp(pp$x, pp$y, window = win, marks = factor(marks))
}

test_that("kamp_variance_biv returns expected columns and defaults to translational correction", {
  marked_pp <- make_biv_pp()
  rvals <- c(0.05, 0.1)

  result <- kamp_variance_biv(marked_pp, rvals = rvals)
  expected_default <- kamp_variance_biv(marked_pp, rvals = rvals, correction = "trans")

  expect_equal(names(result), c("r", "k", "theo_csr", "kamp_csr", "kamp", "var", "pvalue"))
  expect_equal(nrow(result), length(rvals))
  expect_equal(result, expected_default)
})

test_that("kamp_variance_biv respects custom rvals and mark1/mark2", {
  marked_pp <- make_biv_pp()
  rvals <- c(0.02, 0.05, 0.08)

  result <- kamp_variance_biv(marked_pp, rvals = rvals,
                               mark1 = "immune1", mark2 = "immune2")
  expect_equal(result$r, rvals)
  expect_true(all(result$pvalue >= 0 & result$pvalue <= 1))
})

test_that("kamp_variance_biv_Rcpp matches kamp_variance_biv for translational correction", {
  marked_pp <- make_biv_pp(seed = 9)
  rvals <- c(0.05, 0.1, 0.15)

  base_result <- kamp_variance_biv(marked_pp, rvals = rvals)
  rcpp_result <- kamp_variance_biv_Rcpp(marked_pp, rvals = rvals)

  expect_equal(names(rcpp_result), names(base_result))
  expect_equal(rcpp_result$k, base_result$k, tolerance = 1e-8)
  expect_equal(rcpp_result$kamp_csr, base_result$kamp_csr, tolerance = 1e-8)
  expect_equal(rcpp_result$var, base_result$var, tolerance = 1e-8)
})

test_that("kamp dispatches to kamp_variance_biv for bivariate variance calls", {
  set.seed(10)
  df <- data.frame(x = runif(100), y = runif(100),
                    phenotype = sample(c("immune1", "immune2", "background"), 100, replace = TRUE))

  result <- kamp(df, rvals = c(0.1, 0.2), univariate = FALSE, mark_var = "phenotype",
                  mark1 = "immune1", mark2 = "immune2", variance = TRUE)

  expect_true(all(c("r", "k", "kamp_csr", "var", "pvalue") %in% names(result)))
  expect_equal(nrow(result), 2)
})

test_that("kamp_variance_biv matches kamp_variance_biv_helper and accepts the isotropic alias", {
  marked_pp <- make_biv_pp(seed = 11)
  rvals <- c(0.05, 0.1)

  iso_result <- kamp_variance_biv(marked_pp, rvals = rvals, correction = "iso")
  isotropic_result <- purrr::map_dfr(
    rvals,
    ~kamp_variance_biv_helper(marked_pp, rval = .x, correction = "isotropic")
  )

  expect_equal(iso_result, isotropic_result)
  expect_true(all(iso_result$var >= 0))
  expect_true(all(iso_result$pvalue >= 0 & iso_result$pvalue <= 1))
})

test_that("kamp_variance_biv_Rcpp matches kamp_variance_biv for isotropic correction", {
  marked_pp <- make_biv_pp(seed = 12, lambda = 120)
  rvals <- c(0.02, 0.05, 0.1)

  base_result <- kamp_variance_biv(marked_pp, rvals = rvals, correction = "iso")
  rcpp_result <- kamp_variance_biv_Rcpp(marked_pp, rvals = rvals, correction = "iso")

  expect_equal(rcpp_result$k, base_result$k, tolerance = 1e-6)
  expect_equal(rcpp_result$kamp_csr, base_result$kamp_csr, tolerance = 1e-6)
  expect_equal(rcpp_result$var, base_result$var, tolerance = 1e-6)
})

test_that("kamp_variance_biv_Rcpp matches kamp_variance_biv for no edge correction", {
  marked_pp <- make_biv_pp(seed = 13, lambda = 120)
  rvals <- c(0.02, 0.05, 0.1)

  base_result <- kamp_variance_biv(marked_pp, rvals = rvals, correction = "none")
  rcpp_result <- kamp_variance_biv_Rcpp(marked_pp, rvals = rvals, correction = "none")

  expect_equal(rcpp_result$k, base_result$k, tolerance = 1e-8)
  expect_equal(rcpp_result$kamp_csr, base_result$kamp_csr, tolerance = 1e-8)
  expect_equal(rcpp_result$var, base_result$var, tolerance = 1e-8)
})

test_that("kamp() bivariate variance accepts iso/none correction", {
  set.seed(14)
  df <- data.frame(x = runif(100), y = runif(100),
                    phenotype = sample(c("immune1", "immune2", "background"), 100, replace = TRUE))

  iso_result <- kamp(df, rvals = c(0.1, 0.2), univariate = FALSE, mark_var = "phenotype",
                      mark1 = "immune1", mark2 = "immune2", variance = TRUE, correction = "iso")
  expect_true(all(c("r", "k", "kamp_csr", "var", "pvalue") %in% names(iso_result)))
  expect_true(all(iso_result$pvalue >= 0 & iso_result$pvalue <= 1))

  none_result <- kamp(df, rvals = c(0.1, 0.2), univariate = FALSE, mark_var = "phenotype",
                       mark1 = "immune1", mark2 = "immune2", variance = TRUE, correction = "none")
  expect_true(all(c("r", "k", "kamp_csr", "var", "pvalue") %in% names(none_result)))
  expect_true(all(none_result$pvalue >= 0 & none_result$pvalue <= 1))
})
