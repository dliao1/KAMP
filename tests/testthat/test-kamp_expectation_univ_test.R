make_univ_pp <- function(seed = 1, lambda = 150, p_immune = 0.4) {
  set.seed(seed)
  win <- spatstat.geom::owin(c(0, 1), c(0, 1))
  pp <- spatstat.random::rpoispp(lambda = lambda, win = win)
  marks <- sample(c("immune", "background"), pp$n, replace = TRUE,
                   prob = c(p_immune, 1 - p_immune))
  spatstat.geom::ppp(pp$x, pp$y, window = win, marks = factor(marks))
}

test_that("kamp_expectation returns expected columns and defaults to translational correction", {
  marked_pp <- make_univ_pp()
  rvals <- c(0, 0.05, 0.1)

  result <- kamp_expectation(marked_pp, rvals = rvals)
  expected_default <- kamp_expectation(marked_pp, rvals = rvals, correction = "trans")

  expect_equal(names(result), c("r", "k", "theo_csr", "kamp_csr", "kamp"))
  expect_equal(nrow(result), length(rvals))
  expect_true(is.numeric(result$k))
  expect_equal(result, expected_default)
})

test_that("kamp_expectation respects a custom rvals vector and iso correction", {
  marked_pp <- make_univ_pp()
  rvals <- c(0, 0.02, 0.05, 0.08)

  result <- kamp_expectation(marked_pp, rvals = rvals, correction = "iso")
  expect_equal(nrow(result), length(rvals))
  expect_equal(result$r, rvals)
  expect_true(all(!is.na(result$kamp_csr)))
})

test_that("kamp_expectation supports 'none' (uncorrected) and kamp() accepts full-name aliases", {
  marked_pp <- make_univ_pp()
  rvals <- c(0, 0.05, 0.1)

  none_result <- kamp_expectation(marked_pp, rvals = rvals, correction = "none")
  expect_equal(names(none_result), c("r", "k", "theo_csr", "kamp_csr", "kamp"))
  expect_equal(nrow(none_result), length(rvals))
  expect_true(all(!is.na(none_result$kamp_csr)))

  translational_result <- kamp(marked_pp, rvals = rvals, mark1 = "immune",
                                correction = "translational")
  trans_result <- kamp(marked_pp, rvals = rvals, mark1 = "immune", correction = "trans")
  expect_equal(translational_result, trans_result)

  isotropic_result <- kamp(marked_pp, rvals = rvals, mark1 = "immune",
                            correction = "isotropic")
  iso_result <- kamp(marked_pp, rvals = rvals, mark1 = "immune", correction = "iso")
  expect_equal(isotropic_result, iso_result)
})

test_that("check_inputs notes suggested correction for very large point patterns", {
  win <- spatstat.geom::owin(c(0, 1), c(0, 1))
  n <- 100001
  big_pp <- spatstat.geom::ppp(runif(n), runif(n), window = win,
                                marks = factor(sample(c("immune", "background"), n, replace = TRUE)))

  expect_message(
    suppressWarnings(
      check_inputs(big_pp, rvals = c(0, 0.01), univariate = TRUE, correction = "trans",
                   mark_var = NULL, mark1 = "immune", mark2 = NULL, variance = FALSE,
                   thin = FALSE, p_thin = 0)
    ),
    "more than 100,000 points"
  )
})

test_that("check_inputs stays silent about large-N correction below the 100,000 threshold", {
  win <- spatstat.geom::owin(c(0, 1), c(0, 1))
  n <- 100000
  pp <- spatstat.geom::ppp(runif(n), runif(n), window = win,
                            marks = factor(sample(c("immune", "background"), n, replace = TRUE)))

  msgs <- character(0)
  withCallingHandlers(
    suppressWarnings(
      check_inputs(pp, rvals = c(0, 0.01), univariate = TRUE, correction = "trans",
                   mark_var = NULL, mark1 = "immune", mark2 = NULL, variance = FALSE,
                   thin = FALSE, p_thin = 0)
    ),
    message = function(m) { msgs[[length(msgs) + 1]] <<- conditionMessage(m); invokeRestart("muffleMessage") }
  )

  expect_false(any(grepl("more than 100,000 points", msgs)))
})

test_that("kamp dispatches to kamp_expectation for univariate, non-variance calls", {
  set.seed(2)
  df <- data.frame(x = runif(80), y = runif(80),
                    immune = sample(c("immune", "background"), 80, replace = TRUE))

  result <- kamp(df, rvals = c(0, 0.1, 0.2), mark_var = "immune", mark1 = "immune")

  expect_equal(names(result), c("r", "k", "theo_csr", "kamp_csr", "kamp"))
  expect_equal(nrow(result), 3)
})

test_that("kamp dispatches to kamp_variance when variance = TRUE", {
  set.seed(3)
  df <- data.frame(x = runif(80), y = runif(80),
                    immune = sample(c("immune", "background"), 80, replace = TRUE))

  result <- kamp(df, rvals = c(0.1, 0.2), mark_var = "immune", mark1 = "immune",
                  variance = TRUE)

  expect_true(all(c("r", "k", "kamp_csr", "var", "pvalue") %in% names(result)))
  expect_equal(nrow(result), 2)
})

test_that("kamp validates the input dataframe", {
  expect_error(
    kamp(data.frame(x = 1, y = 1), rvals = c(0, 0.1), mark_var = NULL, mark1 = "immune"),
    "mark_var must be supplied"
  )

  expect_error(
    kamp(data.frame(a = 1, b = 1), rvals = c(0, 0.1), mark_var = "immune", mark1 = "immune"),
    "must contain 'x' and 'y' columns"
  )

  df <- data.frame(x = runif(10), y = runif(10), immune = "immune")
  expect_error(
    kamp(df, rvals = c(0, 0.1), mark_var = "immune", mark1 = "immune"),
    "at least two unique values"
  )
})

test_that("kamp validates mark1/mark2 and thinning arguments", {
  set.seed(4)
  df <- data.frame(x = runif(50), y = runif(50),
                    immune = sample(c("immune", "background"), 50, replace = TRUE))

  expect_error(
    kamp(df, rvals = c(0, 0.1), mark_var = "immune", mark1 = "alien"),
    "mark1 is not a mark in the point pattern object"
  )

  expect_error(
    kamp(df, rvals = c(0, 0.1), mark_var = "immune", mark1 = "immune",
         thin = TRUE, p_thin = "a lot"),
    "p_thin must be numeric"
  )

  expect_error(
    kamp(df, rvals = c(0, 0.1), mark_var = "immune", mark1 = "immune",
         thin = TRUE, p_thin = 1.5),
    "p_thin must be between 0 and 1"
  )

  expect_error(
    kamp(df, rvals = c(0, 0.1), mark_var = "immune", mark1 = "immune", thin = "yes"),
    "Argument 'thin' must be TRUE or FALSE"
  )
})

test_that("kamp errors when bivariate marks are missing or identical", {
  set.seed(5)
  df <- data.frame(x = runif(50), y = runif(50),
                    phenotype = sample(c("immune1", "immune2", "background"), 50, replace = TRUE))

  expect_error(
    kamp(df, rvals = c(0, 0.1), univariate = FALSE, mark_var = "phenotype",
         mark1 = "immune1", mark2 = "immune1"),
    "cannot be the same"
  )

  expect_error(
    kamp(df, rvals = c(0, 0.1), univariate = FALSE, mark_var = "phenotype",
         mark1 = "immune1", mark2 = NULL),
    "Both mark1 and mark2 must be specified"
  )
})
