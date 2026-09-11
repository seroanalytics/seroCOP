test_that("simulation-recovery vignette data are reproducible and internally consistent", {
  true_params <- list(floor = 0.05, ceiling = 0.70, ec50 = 1.5, slope = 2.0)

  simulated <- simulate_serocop_data(
    n = 300,
    floor = true_params$floor,
    ceiling = true_params$ceiling,
    ec50 = true_params$ec50,
    slope = true_params$slope,
    titre_mean = 2.0,
    titre_sd = 1.5,
    seed = 2025
  )

  repeated <- simulate_serocop_data(
    n = 300,
    floor = true_params$floor,
    ceiling = true_params$ceiling,
    ec50 = true_params$ec50,
    slope = true_params$slope,
    titre_mean = 2.0,
    titre_sd = 1.5,
    seed = 2025
  )

  expected_probability <- true_params$floor +
    (true_params$ceiling - true_params$floor) /
      (1 + exp(true_params$slope * (simulated$titre - true_params$ec50)))

  expect_equal(simulated, repeated)
  expect_length(simulated$titre, 300)
  expect_true(all(simulated$infected %in% c(0, 1)))
  expect_true(all(simulated$prob_true >= 0 & simulated$prob_true <= 1))
  expect_equal(simulated$prob_true, expected_probability)
  expect_equal(simulated$params[names(true_params)], true_params)
})

test_that("simulation parameters and Brier-score inputs are validated", {
  expect_error(simulate_serocop_data(floor = 0.8, ceiling = 0.2), "floor must be less")
  expect_error(simulate_serocop_data(slope = 0), "slope must be positive")
  expect_error(simulate_serocop_data(n = 0), "n must be positive")

  expect_equal(brier_score(c(0, 1), c(0.2, 0.8)), 0.04)
  expect_error(brier_score(c(0, 1), 0.5), "same length")
  expect_error(brier_score(c(0, 2), c(0.2, 0.8)), "binary")
  expect_error(brier_score(c(0, 1), c(-0.1, 0.8)), "probabilities")
})