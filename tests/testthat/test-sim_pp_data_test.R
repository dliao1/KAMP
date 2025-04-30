test_that("sim_pp_data works with hom and no clust", {
  pp <- sim_pp_data(lambda_n = 100, abundance = 0.3, clust = FALSE, distribution = "hom")
  expect_s3_class(pp, "ppp")
  expect_true(all(levels(pp$marks) %in% c("immune", "background")))
})

test_that("sim_pp_data works with inhom and no clust", {
  pp <- sim_pp_data(lambda_n = 100, abundance = 0.3, clust = FALSE, distribution = "inhom")
  expect_s3_class(pp, "ppp")
  expect_true(all(levels(pp$marks) %in% c("immune", "background")))

})

test_that("sim_pp_data works with hom and clustering", {
  pp <- sim_pp_data(lambda_n = 100, abundance = 0.3, clust = TRUE, distribution = "hom")
  expect_s3_class(pp, "ppp")
  expect_true(all(levels(pp$marks) %in% c("immune", "background")))
})

test_that("sim_pp_data works with inhom and no clustering", {
  pp <- sim_pp_data(lambda_n = 100, abundance = 0.3, clust = TRUE, distribution = "inhom")
  expect_s3_class(pp, "ppp")
  expect_true(all(levels(pp$marks) %in% c("immune", "background")))
})

test_that("sim_pp_data_biv works with hom and no clust", {
  pp <- sim_pp_data_biv(lambda_n = 100, abundance = 0.3, clust = FALSE, distribution = "hom")
  expect_s3_class(pp, "ppp")
  expect_true(all(levels(pp$marks) %in% c("immune1", "immune2", "background")))
})

test_that("sim_pp_data_biv works with inhom and no clust", {
  pp <- sim_pp_data_biv(lambda_n = 100, abundance = 0.3, clust = FALSE, distribution = "inhom")
  expect_s3_class(pp, "ppp")
  expect_true(all(levels(pp$marks) %in% c("immune1", "immune2", "background")))

})

test_that("sim_pp_data_biv works with hom and clustering", {
  pp <- sim_pp_data_biv(lambda_n = 100, abundance = 0.3, clust = TRUE, distribution = "hom")
  expect_s3_class(pp, "ppp")
  expect_true(all(levels(pp$marks) %in% c("immune1", "immune2", "background")))
})

test_that("sim_pp_data_biv works with inhom and no clustering", {
  pp <- sim_pp_data_biv(lambda_n = 100, abundance = 0.3, clust = TRUE, distribution = "inhom")
  expect_s3_class(pp, "ppp")
  expect_true(all(levels(pp$marks) %in% c("immune1", "immune2", "background")))
})
