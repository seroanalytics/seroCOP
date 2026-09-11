test_that("multi-biomarker vignette data initialize with stable biomarker names", {
  set.seed(2025)
  titre <- cbind(
    IgG = rnorm(40, 2, 1),
    IgA = rnorm(40, 1.5, 1),
    Nonspecific = rnorm(40, 2, 1)
  )
  infected <- rbinom(40, 1, plogis(1 - titre[, "IgG"]))

  model <- SeroCOPCompare$new(titre, infected)

  expect_s3_class(model, "SeroCOPCompare")
  expect_equal(model$biomarker_names, colnames(titre))
  expect_equal(model$titre, titre)
  expect_equal(model$infected, infected)
  expect_error(model$compare_biomarkers(), "have not been fitted")
  expect_error(model$plot_all_curves(), "have not been fitted")
})

test_that("multi-biomarker CoP transformation preserves names and values", {
  model <- SeroCOPCompare$new(matrix(c(1, 2, 3, 4), ncol = 2), c(0, 1))
  risks <- list(IgG = c(0, 0.35, 0.7), IgA = c(0.14, 0.28, 0.42))

  protection <- model$extract_cop_multi(risks, upper_bound = 0.7)

  expect_named(protection, names(risks))
  expect_equal(protection$IgG, c(1, 0.5, 0))
  expect_equal(protection$IgA, c(0.8, 0.6, 0.4))
  expect_error(model$extract_cop_multi(c(0.1, 0.2)), "must be a list")
  expect_error(model$extract_cop_multi(list(IgG = 1.1)), "between 0 and 1")
})