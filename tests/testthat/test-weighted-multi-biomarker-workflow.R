test_that("weighted multi-biomarker vignette accepts matrix and broadcast weights", {
  titre <- matrix(seq_len(24), ncol = 3)
  infected <- rep(c(0, 1), length.out = nrow(titre))
  observation_weights <- seq_len(nrow(titre)) / 10
  weight_matrix <- cbind(observation_weights, observation_weights * 2, observation_weights * 3)

  broadcast_model <- SeroCOPCompare$new(titre, infected, weights = observation_weights)
  weighted_model <- SeroCOPCompare$new(titre, infected, weights = weight_matrix)

  expect_equal(broadcast_model$weights, matrix(observation_weights, nrow(titre), ncol(titre)))
  expect_equal(weighted_model$weights, weight_matrix)
})

test_that("weighted multi-biomarker vignette validates weight dimensions and values", {
  titre <- matrix(seq_len(12), ncol = 2)
  infected <- rep(c(0, 1), 3)

  expect_error(SeroCOPCompare$new(titre, infected, weights = c(1, 2)), "length")
  expect_error(SeroCOPCompare$new(titre, infected, weights = matrix(1, 3, 2)), "same dimensions")
  expect_error(SeroCOPCompare$new(titre, infected, weights = rep(-1, nrow(titre))), "non-negative")
  expect_error(SeroCOPCompare$new(titre, infected, weights = rep(NA_real_, nrow(titre))), "Missing values")
})