## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----setup--------------------------------------------------------------------
library(seroCOP)
library(ggplot2)

## ----simulate-----------------------------------------------------------------
# Set parameters for simulation
set.seed(2025)

true_params <- list(
  floor = 0.05,
  ceiling = 0.70,
  ec50 = 1.5,
  slope = 2.0
)

# Simulate data
sim_data <- simulate_serocop_data(
  n = 300,
  floor = true_params$floor,
  ceiling = true_params$ceiling,
  ec50 = true_params$ec50,
  slope = true_params$slope,
  titre_mean = 2.0,
  titre_sd = 1.5,
  seed = 2025
)

# Examine the simulated data
cat(sprintf("Sample size: %d\n", length(sim_data$titre)))
cat(sprintf("Infection rate: %.1f%%\n", mean(sim_data$infected) * 100))
cat(sprintf("Titre range: [%.2f, %.2f]\n", 
            min(sim_data$titre), max(sim_data$titre)))

## ----plot-sim-data------------------------------------------------------------
# Create a plot showing the true relationship
plot_df <- data.frame(
  titre = sim_data$titre,
  infected = sim_data$infected,
  prob_true = sim_data$prob_true
)

ggplot(plot_df, aes(x = titre, y = prob_true)) +
  geom_line(color = "red", linewidth = 1, alpha = 0.7) +
  geom_point(aes(y = infected), alpha = 0.3) +
  labs(
    title = "Simulated Data: True Infection Probability",
    x = "Antibody Titre (log scale)",
    y = "Probability of Infection"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

## ----init-model---------------------------------------------------------------
# Initialize SeroCOP object
model <- SeroCOP$new(
  titre = sim_data$titre,
  infected = sim_data$infected
)

# View default priors (automatically set based on data)
cat("Default priors:\n")
print(model$priors)

## ----custom-priors, eval=FALSE------------------------------------------------
# # Example: Set custom priors (not run in this vignette)
# model$definePrior(
#   floor_alpha = 2, floor_beta = 18,      # Stronger prior: E[floor] ≈ 0.1 (10% residual risk at high titre)
#   ceiling_alpha = 18, ceiling_beta = 2,  # Stronger prior: E[ceiling] ≈ 0.9 (90% risk at low titre)
#   ec50_mean = 2.0, ec50_sd = 1.0,        # Informative prior on ec50
#   slope_mean = 0, slope_sd = 1           # More conservative slope
# )

## ----fit-model, results='hide', message=FALSE---------------------------------
# Fit the model (using default data-driven priors)
# Note: Using fewer iterations for vignette speed; 
# use more (e.g., iter=2000) for real analysis
model$fit_model(
  chains = 4,
  iter = 1000,
  warmup = 500,
  cores = 1  # Adjust based on your system
)

## ----param-summary------------------------------------------------------------
# Get parameter summary
param_est <- extract_parameters(model, prob = 0.95)
print(param_est)

# Compare to true values
cat("\nTrue Parameters:\n")
cat(sprintf("  floor:   %.3f\n", true_params$floor))
cat(sprintf("  ceiling: %.3f\n", true_params$ceiling))
cat(sprintf("  ec50:    %.3f\n", true_params$ec50))
cat(sprintf("  slope:   %.3f\n", true_params$slope))

## ----recovery-plot------------------------------------------------------------
# Create comparison plot
recovery_df <- data.frame(
  parameter = param_est$parameter,
  estimated = param_est$mean,
  lower = param_est$lower,
  upper = param_est$upper,
  true = c(true_params$floor, true_params$ceiling, 
           true_params$ec50, true_params$slope)
)

ggplot(recovery_df, aes(x = parameter)) +
  geom_pointrange(aes(y = estimated, ymin = lower, ymax = upper),
                  color = "steelblue", size = 1) +
  geom_point(aes(y = true), color = "red", size = 3, shape = 4) +
  labs(
    title = "Parameter Recovery",
    subtitle = "Blue: Estimated (95% CI) | Red X: True Value",
    x = "Parameter",
    y = "Value"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

## ----metrics------------------------------------------------------------------
# Get performance metrics
metrics <- model$get_metrics()

## ----plot-curve---------------------------------------------------------------
model$plot_curve()

## ----cop-calculation----------------------------------------------------------
# Extract protection probabilities directly from Stan model
protection_samples <- model$predict_protection()
correlate_of_protection <- colMeans(protection_samples)

# Also extract risk for comparison
risk_samples <- model$predict()
correlate_of_risk <- colMeans(risk_samples)

# Display summary
summary_df <- data.frame(
  titre = sim_data$titre,
  correlate_of_risk = correlate_of_risk,
  correlate_of_protection = correlate_of_protection
)

head(summary_df, 10)

## ----plot-cop-curve-----------------------------------------------------------
# Create prediction grid
titre_grid <- seq(min(sim_data$titre), max(sim_data$titre), length.out = 100)

# Extract protection probabilities from Stan model
cop_matrix <- model$predict_protection(newdata = titre_grid)

# Calculate summary statistics
cop_mean <- colMeans(cop_matrix)
cop_lower <- apply(cop_matrix, 2, quantile, probs = 0.025)
cop_upper <- apply(cop_matrix, 2, quantile, probs = 0.975)

# Create plot data
cop_plot_df <- data.frame(
  titre = titre_grid,
  cop = cop_mean,
  lower = cop_lower,
  upper = cop_upper
)

# Plot
ggplot(cop_plot_df, aes(x = titre, y = cop)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "steelblue", alpha = 0.3) +
  geom_line(color = "steelblue", linewidth = 1) +
  labs(
    title = "Correlate of Protection Curve",
    x = "Antibody Titre (log scale)",
    y = "Correlate of Protection"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

## ----plot-roc-----------------------------------------------------------------
model$plot_roc()

## ----posterior-predictive-----------------------------------------------------
# Get predictions
predictions <- model$predict()

# Calculate prediction intervals
pred_mean <- colMeans(predictions)
pred_lower <- apply(predictions, 2, quantile, probs = 0.025)
pred_upper <- apply(predictions, 2, quantile, probs = 0.975)

# Plot calibration
calib_df <- data.frame(
  observed = sim_data$infected,
  predicted = pred_mean,
  titre = sim_data$titre
)

# Binned calibration plot
breaks <- quantile(calib_df$titre, probs = seq(0, 1, 0.1))
calib_df$bin <- cut(calib_df$titre, breaks = breaks, include.lowest = TRUE)

calib_summary <- aggregate(
  cbind(observed, predicted) ~ bin,
  data = calib_df,
  FUN = mean
)

ggplot(calib_summary, aes(x = predicted, y = observed)) +
  geom_point(size = 3, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    title = "Calibration Plot",
    x = "Predicted Probability",
    y = "Observed Proportion"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1))

## ----loo-details--------------------------------------------------------------
# Print LOO details
print(model$loo)

## ----session-info-------------------------------------------------------------
sessionInfo()

