# toy pattern: 1500 uniform points, seed chosen so blocking leaves some residual points
set.seed(2)
bk_n <- 1500
bk_x <- runif(bk_n, 0, 10)
bk_y <- runif(bk_n, 0, 10)
bk_mk  <- sample(c("immune", "background"), bk_n, replace = TRUE, prob = c(0.2, 0.8))
bk_mkb <- sample(c("a", "b", "c"), bk_n, replace = TRUE, prob = c(0.2, 0.2, 0.6))
bk_win <- spatstat.geom::owin(c(0, 10), c(0, 10))
bk_pp  <- spatstat.geom::ppp(bk_x, bk_y, window = bk_win, marks = factor(bk_mk))
bk_ppb <- spatstat.geom::ppp(bk_x, bk_y, window = bk_win, marks = factor(bk_mkb))
bk_rel <- c(0.05, 0.1, 0.15)

bk_blocks <- suppressWarnings(bkamp_block(bk_pp, mark1 = "immune"))

test_that("bkamp_block returns non-overlapping rectangular ppp blocks inside the window", {
  expect_type(bk_blocks, "list")
  expect_true(all(vapply(bk_blocks, inherits, logical(1), "ppp")))
  expect_true(all(vapply(bk_blocks, function(p) p$window$type == "rectangle", logical(1))))

  inside <- vapply(bk_blocks, function(p) {
    all(p$window$xrange >= bk_win$xrange[1] - 1e-8 & p$window$xrange <= bk_win$xrange[2] + 1e-8) &&
      all(p$window$yrange >= bk_win$yrange[1] - 1e-8 & p$window$yrange <= bk_win$yrange[2] + 1e-8)
  }, logical(1))
  expect_true(all(inside))

  overlap <- 0
  for (i in seq_along(bk_blocks)) for (j in seq_along(bk_blocks)) if (i < j) {
    a <- spatstat.geom::intersect.owin(bk_blocks[[i]]$window, bk_blocks[[j]]$window, fatal = FALSE)
    if (!is.null(a)) overlap <- max(overlap, spatstat.geom::area(a))
  }
  expect_lt(overlap, 1e-8)
})

test_that("bkamp_block assigns each point to at most one block and keeps its mark", {
  key <- paste(bk_pp$x, bk_pp$y)
  block_key <- unlist(lapply(bk_blocks, function(p) paste(p$x, p$y)))
  expect_false(anyDuplicated(block_key) > 0)
  expect_true(all(block_key %in% key))

  block_marks <- unlist(lapply(bk_blocks, function(p) as.character(p$marks)))
  expect_equal(block_marks, as.character(bk_pp$marks)[match(block_key, key)])
})

test_that("bkamp_block warns about residual points with the right count", {
  n_blocked <- sum(vapply(bk_blocks, spatstat.geom::npoints, numeric(1)))
  expect_warning(bkamp_block(bk_pp, mark1 = "immune"),
                 sprintf("%d points remain in residual cells", bk_pp$n - n_blocked))
})

test_that("every block meets the default minimum counts", {
  m1 <- sum(bk_pp$marks == "immune")
  n_min <- ceiling(c(sqrt(m1), min(sqrt(m1), sqrt(bk_pp$n - m1))))
  counts <- t(vapply(bk_blocks, function(p) c(sum(p$marks == "immune"), sum(p$marks != "immune")), numeric(2)))
  expect_true(all(counts[, 1] >= n_min[1]))
  expect_true(all(counts[, 2] >= n_min[2]))

  blocks_n <- suppressWarnings(bkamp_block(bk_pp, mark1 = "immune", n_min = c(60, 60)))
  counts_n <- t(vapply(blocks_n, function(p) c(sum(p$marks == "immune"), sum(p$marks != "immune")), numeric(2)))
  expect_true(all(counts_n >= 60))
})

test_that("bkamp_block is deterministic", {
  expect_identical(suppressWarnings(bkamp_block(bk_pp, mark1 = "immune")), bk_blocks)
})

test_that("bkamp and bkamp_Rcpp agree", {
  a <- suppressWarnings(bkamp(bk_pp, rvals = bk_rel, mark1 = "immune"))
  b <- suppressWarnings(bkamp_Rcpp(bk_pp, rvals = bk_rel, mark1 = "immune"))
  expect_equal(a, b, tolerance = 1e-8)

  a2 <- suppressWarnings(bkamp(bk_ppb, rvals = bk_rel, univariate = FALSE,
                               mark1 = "a", mark2 = "b"))
  b2 <- suppressWarnings(bkamp_Rcpp(bk_ppb, rvals = bk_rel, univariate = FALSE,
                                    mark1 = "a", mark2 = "b"))
  expect_equal(a2, b2, tolerance = 1e-8)
})

test_that("bkamp_Rcpp combines block z-statistics with weights n_i / p_i", {
  res <- suppressWarnings(bkamp_Rcpp(bk_pp, rvals = bk_rel, mark1 = "immune"))

  z <- sapply(bk_blocks, function(p) {
    radii <- bk_rel * min(diff(p$window$xrange), diff(p$window$yrange))
    k <- kamp_variance_Rcpp(p, rvals = radii, mark1 = "immune")
    k$kamp / sqrt(k$var)
  })
  w <- vapply(bk_blocks, function(p) p$n / (sum(p$marks == "immune") / p$n), numeric(1))
  z_comb <- as.vector(z %*% w) / sqrt(sum(w^2))

  expect_equal(res$r, bk_rel)
  expect_equal(res$z, z_comb, tolerance = 1e-10)
  expect_equal(res$pvalue, pnorm(-z_comb), tolerance = 1e-10)
})

test_that("bkamp_Rcpp keeps the order of unsorted rvals and gives NA at r = 0", {
  sorted <- suppressWarnings(bkamp_Rcpp(bk_pp, rvals = c(0, bk_rel), mark1 = "immune"))
  unsorted <- suppressWarnings(bkamp_Rcpp(bk_pp, rvals = c(0.15, 0, 0.05, 0.1), mark1 = "immune"))
  expect_equal(unsorted$r, c(0.15, 0, 0.05, 0.1))
  expect_equal(unsorted, sorted[c(4, 1, 2, 3), ], ignore_attr = TRUE)
  expect_true(is.na(sorted$z[1]))
})

test_that("bkamp input checks", {
  hull <- spatstat.geom::ppp(bk_x, bk_y,
                             window = spatstat.geom::convexhull.xy(bk_x, bk_y),
                             marks = factor(bk_mk))
  expect_error(bkamp_block(hull, mark1 = "immune"), "rectangular observation window")
  expect_error(bkamp_block(bk_pp, mark1 = "nope"), "mark1 is not a mark")
  expect_error(bkamp_block(bk_ppb, univariate = FALSE, mark1 = "a"),
               "Both mark1 and mark2 must be specified")
  expect_error(bkamp_block(bk_ppb, univariate = FALSE, mark1 = "a", mark2 = "a"),
               "cannot be the same")
  expect_error(bkamp_block(bk_pp, mark1 = "immune", rho1 = 3, rho2 = 2), "rho2 must be at least rho1")
  expect_error(bkamp_block(bk_pp, mark1 = "immune", split_max = 0), "split_max must be finite")
  expect_error(bkamp_block(bk_pp, mark1 = "immune", n_min = c(1, 2, 3)), "n_min must contain")
  expect_error(bkamp_block(bk_pp, mark1 = "immune", n_min = c(1e6, 1)), "n_min exceeds")
  expect_error(bkamp(bk_pp, rvals = c(0.1, 0.1), mark1 = "immune"), "distinct nonnegative")
  expect_error(bkamp(bk_pp, rvals = -0.1, mark1 = "immune"), "distinct nonnegative")
  expect_error(bkamp_Rcpp(bk_ppb, univariate = FALSE, correction = "border",
                          mark1 = "a", mark2 = "b"),
               "Unsupported edge correction")
})
