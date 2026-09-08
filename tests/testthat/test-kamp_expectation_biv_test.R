make_biv_pp <- function(seed = 1, lambda = 150) {
  set.seed(seed)
  win <- spatstat.geom::owin(c(0, 1), c(0, 1))
  pp <- spatstat.random::rpoispp(lambda = lambda, win = win)
  marks <- sample(c("immune1", "immune2", "background"), pp$n, replace = TRUE)
  spatstat.geom::ppp(pp$x, pp$y, window = win, marks = factor(marks))
}

test_that("kamp_expectation_biv returns expected columns and defaults to translational correction", {
  marked_pp <- make_biv_pp()
  rvals <- c(0, 0.05, 0.1)

  result <- kamp_expectation_biv(marked_pp, rvals = rvals)
  expected_default <- kamp_expectation_biv(marked_pp, rvals = rvals, correction = "trans")

  expect_equal(names(result), c("r", "k", "theo_csr", "kamp_csr", "kamp"))
  expect_equal(nrow(result), length(rvals))
  expect_equal(result, expected_default)
})

test_that("kamp_expectation_biv respects custom rvals and mark1/mark2", {
  marked_pp <- make_biv_pp()
  rvals <- c(0, 0.02, 0.05)

  result <- kamp_expectation_biv(marked_pp, rvals = rvals,
                                  mark1 = "immune1", mark2 = "immune2")
  expect_equal(nrow(result), length(rvals))
  expect_true(is.numeric(result$k))
})

test_that("kamp dispatches to kamp_expectation_biv for bivariate, non-variance calls", {
  set.seed(6)
  df <- data.frame(x = runif(100), y = runif(100),
                    phenotype = sample(c("immune1", "immune2", "background"), 100, replace = TRUE))

  result <- kamp(df, rvals = c(0, 0.1, 0.2), univariate = FALSE, mark_var = "phenotype",
                  mark1 = "immune1", mark2 = "immune2")

  expect_equal(names(result), c("r", "k", "theo_csr", "kamp_csr", "kamp"))
  expect_equal(nrow(result), 3)
})
