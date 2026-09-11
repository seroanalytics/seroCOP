test_that("overview vignette workflow initializes single, weighted, and grouped models", {
  simulated <- simulate_serocop_data(n = 30, seed = 123)
  groups <- rep(c("Younger", "Older"), each = 15)
  weights <- rep(c(1, 2), length.out = 30)

  model <- NULL
  expect_message(
    model <- SeroCOP$new(simulated$titre, simulated$infected),
    "SeroCOP initialized"
  )
  weighted_model <- SeroCOP$new(simulated$titre, simulated$infected, weights = weights)
  grouped_model <- SeroCOP$new(simulated$titre, simulated$infected, group = groups)

  expect_s3_class(model, "SeroCOP")
  expect_equal(weighted_model$weights, weights)
  expect_equal(levels(grouped_model$group), c("Older", "Younger"))
  expect_equal(model$priors$ec50_mean, mean(range(simulated$titre)))

  model$definePrior(ec50_mean = 2, ec50_sd = 1, slope_sd = 1)
  expect_equal(model$priors$ec50_mean, 2)
  expect_equal(model$priors$ec50_sd, 1)
  expect_equal(model$priors$slope_sd, 1)
})

test_that("overview vignette inputs reject malformed values", {
  expect_error(SeroCOP$new("titre", c(0, 1)), "titre must be numeric")
  expect_error(SeroCOP$new(c(1, 2), c(0, 2)), "infected must be binary")
  expect_error(SeroCOP$new(c(1, 2), c(0, 1), weights = c(-1, 1)), "non-negative")
  expect_error(SeroCOP$new(c(1, 2), c(0, 1), group = "A"), "same length")
})