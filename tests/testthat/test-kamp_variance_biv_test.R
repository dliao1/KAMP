test_that("kamp_variance_biv_default", {
  win <- spatstat.geom::owin(c(0, 1), c(0, 1))

  pp <- spatstat.random::rpoispp(lambda = 100, win = win)
  marks <- sample(c("immune1", "immune2", "background"), pp$n, replace = TRUE)
  marked_pp <- spatstat.geom::ppp(pp$x, pp$y, window = win, marks = factor(marks))

  result <- kamp_variance_biv(marked_pp, markvar1 = "immune1", markvar2 = "immune2")

  expect_true(all(c("r", "k", "theo_csr", "kamp_csr", "var", "z", "pvalue") %in% names(result)))
  expect_equal(nrow(result), length(c(0, .05, .075, .1, .15, .2)))
  expect_true(is.numeric(result$k))
  expect_true(any(!is.na(result$k)))

})


test_that("kamp_variance_biv_custom_r_values", {
  win <- spatstat.geom::owin(c(0, 1), c(0, 1))
  pp <- spatstat.random::rpoispp(lambda = 100, win = win)
  marks <- sample(c("immune1", "immune2", "background"), pp$n, replace = TRUE)
  marked_pp <- spatstat.geom::ppp(pp$x, pp$y, window = win, marks = factor(marks))

  rvec_custom <- c(0.01, 0.02, 0.05)
  result <- kamp_variance_biv(marked_pp, rvec = rvec_custom)

  expect_equal(nrow(result), length(rvec_custom))
})

test_that("kamp_variance_biv_thinning", {
  win <- spatstat.geom::owin(c(0, 1), c(0, 1))
  pp <- spatstat.random::rpoispp(lambda = 100, win = win)
  marks <- sample(c("immune1", "immune2", "background"), pp$n, replace = TRUE)
  marked_pp <- spatstat.geom::ppp(pp$x, pp$y, window = win, marks = factor(marks))

  thin_result <- kamp_variance_biv(marked_pp, thin_pct = 0.5)
  full_result <- kamp_variance_biv(marked_pp, thin_pct = 0)

  expect_false(all(thin_result$k == full_result$k))
})


