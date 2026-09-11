test_that("SeroMulti initializes joint models for two through four biomarkers", {
  titre <- cbind(IgG = 1:8, IgA = 11:18)
  infected <- rep(c(0L, 1L), 4)

  model <- SeroMulti$new(titre, infected)

  expect_s3_class(model, "SeroMulti")
  expect_equal(model$biomarker_names, c("IgG", "IgA"))
  expect_equal(model$titre, titre)
  expect_equal(model$priors$slope_scale, 1)
  expect_error(model$predict(), "has not been fitted")
  expect_error(model$surface_data(), "has not been fitted")
})

test_that("SeroMulti validates its joint-model dimensions and inputs", {
  infected <- c(0L, 1L, 0L, 1L)

  expect_error(SeroMulti$new(1:4, infected), "matrix")
  expect_error(SeroMulti$new(matrix(1:4, ncol = 1), infected), "between 2 and 4")
  expect_error(SeroMulti$new(matrix(1:20, ncol = 5), infected), "between 2 and 4")
  expect_error(SeroMulti$new(matrix(1:8, ncol = 2), c(0L, 1L)), "equal")
  expect_error(SeroMulti$new(matrix(1:8, ncol = 2), c(0L, 1L, 0L, 2L)), "binary")
  expect_error(SeroMulti$new(matrix(1:8, ncol = 2), infected, c("A", "A")), "unique")
})

test_that("SeroMulti prior updates validate positive scale parameters", {
  model <- SeroMulti$new(matrix(1:8, ncol = 2), c(0L, 1L, 0L, 1L))

  model$definePrior(ec50_mean = 3, ec50_sd = 2, slope_scale = 1.5)
  expect_equal(model$priors$ec50_mean, 3)
  expect_equal(model$priors$slope_scale, 1.5)
  expect_error(model$definePrior(slope_scale = 0), "must be positive")
})