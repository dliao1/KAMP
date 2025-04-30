test_that("kamp_variance_univ_default", {
  win <- spatstat.geom::owin(c(0, 1), c(0, 1))

  pp <- spatstat.random::rpoispp(lambda = 100, win = win)
  marks <- sample(c("immune","background"), pp$n, replace = TRUE)
  marked_pp <- spatstat.geom::ppp(pp$x, pp$y, window = win, marks = factor(marks))

  result <- kamp_variance(marked_pp, markvar = "immune")

  expect_true(all(c("r", "k", "theo_csr", "kamp_csr", "var", "z", "pvalue") %in% names(result)))
  expect_equal(nrow(result), length(c(0, .05, .075, .1, .15, .2)))
  expect_true(is.numeric(result$k))
  expect_true(any(!is.na(result$k)))

})


test_that("kamp_variance_univ_custom_r_values", {
  win <- spatstat.geom::owin(c(0, 1), c(0, 1))
  pp <- spatstat.random::rpoispp(lambda = 100, win = win)
  marks <- sample(c("immune", "background"), pp$n, replace = TRUE)
  marked_pp <- spatstat.geom::ppp(pp$x, pp$y, window = win, marks = factor(marks))

  rvec_custom <- c(0.01, 0.02, 0.05)
  result <- kamp_variance(marked_pp, rvec = rvec_custom)

  expect_equal(nrow(result), length(rvec_custom))
})

test_that("kamp_variance_univ_thinning", {
  win <- spatstat.geom::owin(c(0, 1), c(0, 1))
  pp <- spatstat.random::rpoispp(lambda = 100, win = win)
  marks <- sample(c("immune", "background"), pp$n, replace = TRUE)
  marked_pp <- spatstat.geom::ppp(pp$x, pp$y, window = win, marks = factor(marks))

  thin_result <- kamp_variance(marked_pp, thin_pct = 0.5)
  full_result <- kamp_variance(marked_pp, thin_pct = 0)

  expect_false(all(thin_result$k == full_result$k))
})

test_that("kamp_variance_univ_errors", {
  win <- spatstat.geom::owin(c(0, 1), c(0, 1))
  pp <- spatstat.random::rpoispp(lambda = 100, win = win)
  marks <- sample(c("immune", "background"), pp$n, replace = TRUE)
  marked_pp <- spatstat.geom::ppp(pp$x, pp$y, window = win, marks = factor(marks))

  thin_result <- kamp_variance(marked_pp, thin_pct = 0.5)
  full_result <- kamp_variance(marked_pp, thin_pct = 0)

  expect_false(all(thin_result$k == full_result$k))
})

test_that("kamp_variance_univ_errors", {
  win <- spatstat.geom::owin(c(0, 1), c(0, 1))
  pp <- spatstat.random::rpoispp(lambda = 100, win = win)

  expect_error(kamp_variance(pp),
               "The point pattern object must have marks")

  marks <- sample(c("immune", "background"), pp$n, replace = TRUE)
  marked_pp <- spatstat.geom::ppp(pp$x, pp$y, window = win, marks = factor(marks))
  expect_error(kamp_variance(marked_pp, correction = "border"),
               "Currently only isotropic and translational edge correction are supported")

  expect_error(kamp_variance(marked_pp, thin_pct = "a lot"),
               "thin_pct must be numeric")

  expect_error(kamp_variance(marked_pp, thin_pct = -0.1),
               "thin_pct must be between 0 and 1")
  expect_error(kamp_variance(marked_pp, thin_pct = 1.5),
               "thin_pct must be between 0 and 1")

  expect_error(kamp_variance(marked_pp, markvar = "alien"),
               "markvar is not a mark in the point pattern object")
})



