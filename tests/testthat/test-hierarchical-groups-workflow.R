test_that("hierarchical-groups vignette data retain group labels and sample sizes", {
  set.seed(2025)
  group_sizes <- c(Young = 15, Middle = 9, Old = 6)
  group <- rep(names(group_sizes), group_sizes)
  titre <- rnorm(sum(group_sizes), mean = rep(c(2.2, 1.8, 1.4), group_sizes), sd = 0.8)
  infected <- rbinom(length(titre), 1, plogis(1 - titre))

  model <- SeroCOP$new(titre, infected, group = group)
  binned <- aggregate(infected ~ group + bin,
    data = transform(data.frame(group, titre, infected), bin = cut(titre, breaks = 3)),
    FUN = mean
  )

  expect_equal(as.integer(table(model$group)), unname(group_sizes[levels(model$group)]))
  expect_equal(nlevels(model$group), 3)
  expect_true(all(binned$infected >= 0 & binned$infected <= 1))
  expect_error(model$extract_group_parameters(), "has not been fitted")
})