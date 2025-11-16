# Example: Using seroCOP for Correlates of Protection Analysis
# This script demonstrates the basic workflow

library(seroCOP)

# 1. Simulate Data --------------------------------------------------------

# Set seed for reproducibility
set.seed(2025)

# Simulate data with known parameters
sim_data <- simulate_serocop_data(
  n = 300,
  floor = 0.05,      # 5% baseline infection risk
  ceiling = 0.90,    # 90% maximum infection risk
  ec50 = 1.5,        # Titre at 50% infection probability
  slope = 2.0,       # Steepness of the curve
  titre_mean = 2.0,
  titre_sd = 1.5,
  seed = 2025
)

# Examine the simulated data
cat("Simulated data summary:\n")
cat(sprintf("  Sample size: %d\n", length(sim_data$titre)))
cat(sprintf("  Infection rate: %.1f%%\n", mean(sim_data$infected) * 100))
cat(sprintf("  Titre range: [%.2f, %.2f]\n", 
            min(sim_data$titre), max(sim_data$titre)))


# 2. Initialize SeroCOP Object --------------------------------------------

model <- SeroCOP$new(
  titre = sim_data$titre,
  infected = sim_data$infected
)


# 3. Fit the Bayesian Model -----------------------------------------------

# Fit the model with Stan
# Note: Adjust cores based on your system
model$fit_model(
  chains = 4,
  iter = 2000,
  warmup = 1000,
  cores = 4
)


# 4. Examine Parameter Estimates ------------------------------------------

# Get parameter summary from the model
cat("\nParameter estimates:\n")
param_summary <- model$summary()
print(param_summary)

# Extract parameters with credible intervals
params <- extract_parameters(model, prob = 0.95)
print(params)

# Compare to true values
cat("\nTrue parameters:\n")
print(sim_data$params[c("floor", "ceiling", "ec50", "slope")])


# 5. Calculate Performance Metrics ----------------------------------------

# Get all performance metrics
metrics <- model$get_metrics()

# Individual metrics are also accessible
cat(sprintf("\nROC AUC: %.3f\n", metrics$roc_auc))
cat(sprintf("Brier Score: %.3f\n", metrics$brier_score))


# 6. Generate Predictions -------------------------------------------------

# Predictions for original data
predictions <- model$predict()
pred_mean <- colMeans(predictions)

# Predictions for new data
new_titres <- seq(min(sim_data$titre), max(sim_data$titre), length.out = 100)
new_predictions <- model$predict(newdata = new_titres)
new_pred_mean <- colMeans(new_predictions)


# 7. Visualize Results ----------------------------------------------------

# Plot the fitted curve with uncertainty
p1 <- model$plot_curve()
print(p1)

# Plot ROC curve
p2 <- model$plot_roc()
print(p2)

# Plot parameter posterior distributions
# (requires tidyr package)
if (requireNamespace("tidyr", quietly = TRUE)) {
  p3 <- model$plot_parameters()
  print(p3)
}


# 8. Model Diagnostics ----------------------------------------------------

# Print LOO-CV results
print(model$loo)

# Access the Stan fit object directly for more detailed diagnostics
# Check convergence (Rhat should be < 1.01)
stan_summary <- rstan::summary(model$fit)
rhat_values <- stan_summary$summary[, "Rhat"]
cat(sprintf("\nMax Rhat: %.3f\n", max(rhat_values, na.rm = TRUE)))

# Check effective sample size
neff_values <- stan_summary$summary[, "n_eff"]
cat(sprintf("Min n_eff: %.0f\n", min(neff_values, na.rm = TRUE)))


# 9. Save Results (Optional) ----------------------------------------------

# Save the fitted model
# saveRDS(model, "serocop_model.rds")

# Save plots
# ggsave("curve_plot.png", p1, width = 8, height = 6)
# ggsave("roc_plot.png", p2, width = 6, height = 6)

# Save parameter estimates
# write.csv(params, "parameter_estimates.csv", row.names = FALSE)

cat("\nAnalysis complete!\n")
