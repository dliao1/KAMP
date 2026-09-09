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

test_that("kamp_expectation_biv matches spatstat Kcross/Kest for iso correction", {
  marked_pp <- make_biv_pp()
  rvals <- c(0, 0.05, 0.1)

  iso_result <- kamp_expectation_biv(marked_pp, rvals = rvals, correction = "iso",
                                      mark1 = "immune1", mark2 = "immune2")
  iso_k <- spatstat.explore::Kcross(marked_pp, i = "immune1", j = "immune2",
                                     r = rvals, correction = "iso")
  iso_csr <- spatstat.explore::Kest(marked_pp, r = rvals, correction = "iso")

  expect_equal(iso_result$k, iso_k$iso)
  expect_equal(iso_result$kamp_csr, iso_csr$iso)
  expect_equal(iso_result$kamp, iso_result$k - iso_result$kamp_csr)
})

test_that("kamp_expectation_biv matches spatstat Kcross/Kest for no edge correction", {
  marked_pp <- make_biv_pp()
  rvals <- c(0, 0.05, 0.1)

  none_result <- kamp_expectation_biv(marked_pp, rvals = rvals, correction = "none",
                                       mark1 = "immune1", mark2 = "immune2")
  none_k <- spatstat.explore::Kcross(marked_pp, i = "immune1", j = "immune2",
                                      r = rvals, correction = "none")
  none_csr <- spatstat.explore::Kest(marked_pp, r = rvals, correction = "none")

  expect_equal(none_result$k, none_k$un)
  expect_equal(none_result$kamp_csr, none_csr$un)
  expect_equal(none_result$kamp, none_result$k - none_result$kamp_csr)
})

test_that("kamp() bivariate accepts iso/none correction and full-name aliases", {
  set.seed(6)
  df <- data.frame(x = runif(100), y = runif(100),
                    phenotype = sample(c("immune1", "immune2", "background"), 100, replace = TRUE))

  iso_result <- kamp(df, rvals = c(0, 0.1), univariate = FALSE, mark_var = "phenotype",
                      mark1 = "immune1", mark2 = "immune2", correction = "iso")
  isotropic_result <- kamp(df, rvals = c(0, 0.1), univariate = FALSE, mark_var = "phenotype",
                            mark1 = "immune1", mark2 = "immune2", correction = "isotropic")
  expect_equal(iso_result, isotropic_result)

  none_result <- kamp(df, rvals = c(0, 0.1), univariate = FALSE, mark_var = "phenotype",
                       mark1 = "immune1", mark2 = "immune2", correction = "none")
  expect_equal(names(none_result), c("r", "k", "theo_csr", "kamp_csr", "kamp"))
  expect_true(all(!is.na(none_result$kamp_csr)))
})

test_that("check_inputs notes suggested correction for very large bivariate point patterns", {
  win <- spatstat.geom::owin(c(0, 1), c(0, 1))
  n <- 100001
  big_pp <- spatstat.geom::ppp(runif(n), runif(n), window = win,
                                marks = factor(sample(c("immune1", "immune2", "background"),
                                                       n, replace = TRUE)))

  expect_message(
    suppressWarnings(
      check_inputs(big_pp, rvals = c(0, 0.01), univariate = FALSE, correction = "trans",
                   mark_var = NULL, mark1 = "immune1", mark2 = "immune2", variance = FALSE,
                   thin = FALSE, p_thin = 0)
    ),
    "more than 100,000 points"
  )
})

test_that("kamp() surfaces the >100,000 points message for bivariate calls", {
  set.seed(13)
  n <- 100001
  win <- spatstat.geom::owin(c(0, 1), c(0, 1))
  big_pp <- spatstat.geom::ppp(runif(n), runif(n), window = win,
                                marks = factor(sample(c("immune1", "immune2", "background"),
                                                       n, replace = TRUE)))

  expect_message(
    result <- kamp(big_pp, rvals = c(0, 0.05), univariate = FALSE,
                    mark1 = "immune1", mark2 = "immune2", correction = "none",
                    thin = TRUE, p_thin = 0.995),
    "more than 100,000 points"
  )

  expect_equal(names(result), c("r", "k", "theo_csr", "kamp_csr", "kamp"))
  expect_equal(nrow(result), 2)
})
