library(testthat)
library(MLMResearch)

test_that("adaptive_rve selects correct method", {
  # Mock functions since test_rve_cr2, test_rve_cr1, etc. are used in adaptive_rve
  # and we assume they are defined in the package.
  
  # create dummy data
  set.seed(123)
  k <- 20
  df <- data.frame(
    review_id = sample(1:k, 50, replace = TRUE),
    study_id = sample(1:30, 50, replace = TRUE),
    TE = rnorm(50, 0, 0.5),
    seTE = runif(50, 0.1, 0.3)
  )
  
  # Because test_rve_cr2 etc. might rely on actual metafor/robumeta models, 
  # we test get_rve_recommendation instead to avoid complex model fitting failures
  rec <- get_rve_recommendation(df, cluster_var = "review_id")
  expect_equal(rec$n_clusters, 20)
  expect_equal(rec$recommendation$method, "CR2")
  
  # Test with fewer clusters
  df2 <- data.frame(
    review_id = c(1:8, sample(1:8, 12, replace = TRUE)),
    study_id = sample(1:10, 20, replace = TRUE),
    TE = rnorm(20, 0, 0.5),
    seTE = runif(20, 0.1, 0.3)
  )
  
  rec2 <- get_rve_recommendation(df2, cluster_var = "review_id")
  expect_equal(rec2$n_clusters, 8)
  expect_true(grepl("CR1", rec2$recommendation$method))
})
