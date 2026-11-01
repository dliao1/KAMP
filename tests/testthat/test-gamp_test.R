# tiny pattern so every relabeling can be enumerated
# univariate: choose(10, 4) = 210; bivariate: choose(10, 3) * choose(7, 3) = 4200
set.seed(333)
N <- 10
gamp_df <- data.frame(x = runif(N, 0, 10),
                      y = runif(N, 0, 10),
                      marks = sample(c(rep("pos", 4), rep("neg", N - 4))))
gamp_df$marks3 <- sample(c(rep("a", 3), rep("b", 3), rep("c", 4)))
gamp_rvals <- seq(0, 8, 0.5)

gamp_D <- as.matrix(dist(gamp_df[, c("x", "y")]))
diag(gamp_D) <- Inf

# G at each radius for every relabeling: distance from each mark1 point to its nearest mark2 point
enum_G <- function(ref_sets) {
  sapply(ref_sets, function(s) {
    nn <- apply(gamp_D[s$ref, s$tgt, drop = FALSE], 1, min)
    vapply(gamp_rvals, function(r) mean(nn <= r), numeric(1))
  })
}
pop_var <- function(G) apply(G, 1, function(v) mean((v - mean(v))^2))

test_that("univariate gamp matches exact enumeration of all relabelings", {
  res <- suppressMessages(gamp(gamp_df, gamp_rvals, univariate = TRUE, mark_var = "marks",
                               mark1 = "pos", variance = TRUE))

  sets <- lapply(utils::combn(N, 4, simplify = FALSE), function(s) list(ref = s, tgt = s))
  G_all <- enum_G(sets)

  expect_equal(res$gamp_csr, rowMeans(G_all), tolerance = 1e-10)
  expect_equal(res$var, pop_var(G_all), tolerance = 1e-10)
})

test_that("bivariate gamp matches exact enumeration of all relabelings", {
  res <- suppressMessages(gamp(gamp_df, gamp_rvals, univariate = FALSE, mark_var = "marks3",
                               mark1 = "a", mark2 = "b", variance = TRUE))

  sets <- list()
  for (s1 in utils::combn(N, 3, simplify = FALSE)) {
    rest <- setdiff(seq_len(N), s1)
    for (s2 in utils::combn(rest, 3, simplify = FALSE)) {
      sets[[length(sets) + 1]] <- list(ref = s1, tgt = s2)
    }
  }
  G_all <- enum_G(sets)

  expect_equal(res$gamp_csr, rowMeans(G_all), tolerance = 1e-10)
  expect_equal(res$var, pop_var(G_all), tolerance = 1e-10)
})

test_that("gamp observed G matches Gest / Gcross with no correction", {
  pp <- spatstat.geom::ppp(gamp_df$x, gamp_df$y,
                           window = spatstat.geom::convexhull.xy(gamp_df$x, gamp_df$y),
                           marks = factor(gamp_df$marks3))

  res_u <- suppressMessages(gamp(pp, gamp_rvals, mark1 = "a"))
  ge <- spatstat.explore::Gest(pp[pp$marks == "a"], r = gamp_rvals, correction = "none")
  expect_equal(res_u$g, ge$raw)
  expect_equal(res_u$theo_csr, ge$theo)
  expect_false("var" %in% names(res_u))

  res_b <- suppressMessages(gamp(pp, gamp_rvals, univariate = FALSE, mark1 = "a", mark2 = "b"))
  gc <- spatstat.explore::Gcross(pp, i = "a", j = "b", r = gamp_rvals, correction = "none")
  expect_equal(res_b$g, gc$raw)
  expect_equal(res_b$theo_csr, gc$theo)
})

test_that("gamp dispatches to the matching expectation / variance function", {
  pp <- spatstat.geom::ppp(gamp_df$x, gamp_df$y,
                           window = spatstat.geom::convexhull.xy(gamp_df$x, gamp_df$y),
                           marks = factor(gamp_df$marks3))

  expect_equal(suppressMessages(gamp(pp, gamp_rvals, mark1 = "a")),
               gamp_expectation(pp, gamp_rvals, mark1 = "a"))
  expect_equal(suppressMessages(gamp(pp, gamp_rvals, mark1 = "a", variance = TRUE)),
               gamp_variance(pp, gamp_rvals, mark1 = "a"))
  expect_equal(suppressMessages(gamp(pp, gamp_rvals, univariate = FALSE, mark1 = "a", mark2 = "b")),
               gamp_expectation_biv(pp, gamp_rvals, mark1 = "a", mark2 = "b"))
  expect_equal(suppressMessages(gamp(pp, gamp_rvals, univariate = FALSE, mark1 = "a", mark2 = "b",
                                     variance = TRUE)),
               gamp_variance_biv(pp, gamp_rvals, mark1 = "a", mark2 = "b"))
})

test_that("gamp expectation and variance functions agree on the mean", {
  pp <- spatstat.geom::ppp(gamp_df$x, gamp_df$y,
                           window = spatstat.geom::convexhull.xy(gamp_df$x, gamp_df$y),
                           marks = factor(gamp_df$marks3))

  expect_equal(gamp_variance(pp, gamp_rvals, mark1 = "a")$gamp_csr,
               gamp_expectation(pp, gamp_rvals, mark1 = "a")$gamp_csr)
  expect_equal(gamp_variance_biv(pp, gamp_rvals, mark1 = "a", mark2 = "b")$gamp_csr,
               gamp_expectation_biv(pp, gamp_rvals, mark1 = "a", mark2 = "b")$gamp_csr)
})

test_that("gamp errors with too few marked points", {
  df1 <- gamp_df
  df1$marks <- c("pos", rep("neg", N - 1))
  expect_error(suppressMessages(gamp(df1, gamp_rvals, mark_var = "marks", mark1 = "pos")),
               "at least two")
  expect_error(suppressMessages(gamp(gamp_df, gamp_rvals, univariate = FALSE, mark_var = "marks3",
                                     mark1 = "a", mark2 = "a")),
               "cannot be the same")
})

test_that("gamp only accepts correction = 'none'", {
  expect_error(gamp(gamp_df, gamp_rvals, mark_var = "marks", mark1 = "pos", correction = "trans"),
               "only supports correction = 'none'")
  pp <- spatstat.geom::ppp(gamp_df$x, gamp_df$y,
                           window = spatstat.geom::convexhull.xy(gamp_df$x, gamp_df$y),
                           marks = factor(gamp_df$marks))
  expect_error(gamp_variance(pp, gamp_rvals, mark1 = "pos", correction = "iso"),
               "only supports correction = 'none'")
})

test_that("gamp errors when a data.frame is passed without mark_var", {
  expect_error(suppressMessages(gamp(gamp_df, gamp_rvals, mark1 = "pos")),
               "mark_var must be supplied")
})

test_that("univariate gamp variance matches PermFree-NearestNeighborG exact_G_variance", {
  set.seed(2024)
  df <- data.frame(x = runif(80, 0, 10), y = runif(80, 0, 10))
  df$marks <- sample(c("pos", "neg"), 80, replace = TRUE, prob = c(0.3, 0.7))
  radii <- c(0.5, 1, 1.5, 2, 3)

  res <- suppressMessages(gamp(df, radii, mark_var = "marks", mark1 = "pos", variance = TRUE))

  expect_equal(res$var, exact_G_variance(df, "marks", "pos", radii), tolerance = 1e-10)
})
