## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 8,
  fig.height = 6
)

## ----setup--------------------------------------------------------------------
library(seroCOP)
library(ggplot2)
set.seed(2025)

## ----simulate-data------------------------------------------------------------
n <- 250  # Sample size

# Biomarker 1: Strong CoP (IgG)
# High ceiling-floor difference, clear dose-response
titre_IgG <- rnorm(n, mean = 2.5, sd = 1.2)
prob_IgG <- 0.02 + 0.68 / (1 + exp(2.5 * (titre_IgG - 2.0)))

# Biomarker 2: Weak CoP (IgA)  
# Moderate ceiling-floor difference, weaker slope
titre_IgA <- rnorm(n, mean = 1.8, sd = 1.5)
prob_IgA <- 0.15 + 0.55 / (1 + exp(1.0 * (titre_IgA - 1.5)))

# Biomarker 3: No CoP (Non-specific)
# No relationship with infection - flat line
titre_Nonspec <- rnorm(n, mean = 3.0, sd = 1.0)
prob_Nonspec <- rep(0.35, n)  # Constant probability

# Generate infection outcomes
# Use weighted average with noise
prob_combined <- 0.5 * prob_IgG + 0.3 * prob_IgA + 0.2 * prob_Nonspec
infected <- rbinom(n, 1, prob_combined)

# Combine into matrix
titre_matrix <- cbind(
  IgG = titre_IgG,
  IgA = titre_IgA,
  Nonspecific = titre_Nonspec
)

# Store true parameters for later comparison
true_params <- list(
  IgG = list(floor = 0.02, ceiling = 0.70, ec50 = 2.0, slope = 2.5),
  IgA = list(floor = 0.15, ceiling = 0.70, ec50 = 1.5, slope = 1.0),
  Nonspecific = list(floor = 0.35, ceiling = 0.70, ec50 = 3.0, slope = 0.01)
)

cat(sprintf("Simulated %d samples with 3 biomarkers\n", n))
cat(sprintf("Overall infection rate: %.1f%%\n", mean(infected) * 100))

## ----plot-true-relationships--------------------------------------------------
# Create visualization of true relationships
plot_data <- data.frame(
  titre = c(titre_IgG, titre_IgA, titre_Nonspec),
  prob = c(prob_IgG, prob_IgA, prob_Nonspec),
  infected = rep(infected, 3),
  biomarker = rep(c("IgG (Strong CoP)", "IgA (Weak CoP)", 
                    "Nonspecific (No CoP)"), each = n)
)

ggplot(plot_data, aes(x = titre, y = prob)) +
  geom_line(color = "red", linewidth = 1, alpha = 0.7) +
  geom_point(aes(y = infected), alpha = 0.3, size = 0.8) +
  facet_wrap(~biomarker, scales = "free_x") +
  labs(
    title = "True Infection Probability by Biomarker",
    x = "Antibody Titre (log scale)",
    y = "Probability of Infection"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

## ----init-multi-model---------------------------------------------------------
# Initialize multi-biomarker model
multi_model <- SeroCOPMulti$new(
  titre = titre_matrix,
  infected = infected,
  biomarker_names = c("IgG", "IgA", "Nonspecific")
)

## ----fit-multi-model, results='hide', message=FALSE---------------------------
# Fit all models
# Note: Using reduced iterations for vignette speed
multi_model$fit_all(
  chains = 4,
  iter = 1000,
  warmup = 500,
  cores = 1
)

## ----plot-all-curves----------------------------------------------------------
multi_model$plot_all_curves()

## ----cop-multi-calculation----------------------------------------------------
# Extract protection probabilities for each biomarker
cop_results <- list()

for (biomarker in multi_model$biomarker_names) {
  model <- multi_model$models[[biomarker]]
  
  # Get protection probabilities directly from Stan model
  protection_samples <- model$predict_protection()
  correlate_of_protection <- colMeans(protection_samples)
  
  # Also get risk for comparison
  risk_samples <- model$predict()
  correlate_of_risk <- colMeans(risk_samples)
  
  cop_results[[biomarker]] <- data.frame(
    titre = model$titre,
    correlate_of_risk = correlate_of_risk,
    correlate_of_protection = correlate_of_protection
  )
}

# Display summary for each biomarker
for (biomarker in multi_model$biomarker_names) {
  cat(sprintf("\n=== %s ===\n", biomarker))
  print(head(cop_results[[biomarker]], 5))
}

## ----plot-cop-curves----------------------------------------------------------
# Combine CoP data for all biomarkers
cop_plot_data <- list()

for (i in seq_along(multi_model$biomarker_names)) {
  biomarker <- multi_model$biomarker_names[i]
  model <- multi_model$models[[biomarker]]
  
  # Create prediction grid
  titre_grid <- seq(min(model$titre), max(model$titre), length.out = 100)
  
  # Extract protection probabilities from Stan model
  cop_matrix <- model$predict_protection(newdata = titre_grid)
  
  # Calculate summary statistics
  cop_mean <- colMeans(cop_matrix)
  cop_lower <- apply(cop_matrix, 2, quantile, probs = 0.025)
  cop_upper <- apply(cop_matrix, 2, quantile, probs = 0.975)
  
  cop_plot_data[[i]] <- data.frame(
    biomarker = biomarker,
    titre = titre_grid,
    cop = cop_mean,
    lower = cop_lower,
    upper = cop_upper
  )
}

cop_plot_df <- do.call(rbind, cop_plot_data)

# Plot all CoP curves
ggplot(cop_plot_df, aes(x = titre, y = cop, color = biomarker, fill = biomarker)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  facet_wrap(~biomarker, scales = "free_x") +
  labs(
    title = "Correlate of Protection Curves by Biomarker",
    x = "Antibody Titre (log scale)",
    y = "Correlate of Protection"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

## ----compare-biomarkers-------------------------------------------------------
comparison <- multi_model$compare_biomarkers()

## ----plot-comparison----------------------------------------------------------
multi_model$plot_comparison()

## ----parameter-recovery-------------------------------------------------------
# Extract parameters for each biomarker
for (biomarker in c("IgG", "IgA", "Nonspecific")) {
  cat(sprintf("\n=== %s ===\n", biomarker))
  
  model <- multi_model$models[[biomarker]]
  params <- extract_parameters(model)
  
  cat("\nEstimated parameters:\n")
  print(params[, c("parameter", "mean", "lower", "upper")])
  
  cat("\nTrue parameters:\n")
  true <- true_params[[biomarker]]
  cat(sprintf("  floor:   %.3f\n", true$floor))
  cat(sprintf("  ceiling: %.3f\n", true$ceiling))
  cat(sprintf("  ec50:    %.3f\n", true$ec50))
  cat(sprintf("  slope:   %.3f\n", true$slope))
}

## ----recovery-visualization, fig.height=8-------------------------------------
# Create combined recovery plot
recovery_list <- list()

for (i in seq_along(multi_model$biomarker_names)) {
  biomarker <- multi_model$biomarker_names[i]
  model <- multi_model$models[[biomarker]]
  params <- extract_parameters(model)
  
  recovery_list[[i]] <- data.frame(
    biomarker = biomarker,
    parameter = params$parameter,
    estimated = params$mean,
    lower = params$lower,
    upper = params$upper,
    true = c(
      true_params[[biomarker]]$floor,
      true_params[[biomarker]]$ceiling,
      true_params[[biomarker]]$ec50,
      true_params[[biomarker]]$slope
    )
  )
}

recovery_df <- do.call(rbind, recovery_list)

ggplot(recovery_df, aes(x = parameter, color = biomarker)) +
  geom_pointrange(
    aes(y = estimated, ymin = lower, ymax = upper),
    position = position_dodge(width = 0.5),
    size = 0.8
  ) +
  geom_point(
    aes(y = true),
    shape = 4,
    size = 3,
    stroke = 1.5,
    position = position_dodge(width = 0.5)
  ) +
  facet_wrap(~parameter, scales = "free_y", ncol = 2) +
  labs(
    title = "Parameter Recovery Across Biomarkers",
    subtitle = "Points with bars: Estimated (95% CI) | X: True value",
    y = "Value",
    color = "Biomarker"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  )

## ----individual-roc, fig.height=4---------------------------------------------
# Plot ROC curves for each biomarker
par(mfrow = c(1, 3))

for (biomarker in multi_model$biomarker_names) {
  model <- multi_model$models[[biomarker]]
  
  pred <- colMeans(model$predict())
  roc_obj <- pROC::roc(model$infected, pred, quiet = TRUE)
  
  plot(roc_obj, 
       main = sprintf("%s\nAUC = %.3f", biomarker, pROC::auc(roc_obj)),
       col = "steelblue",
       lwd = 2)
  abline(a = 0, b = 1, lty = 2, col = "gray")
}

par(mfrow = c(1, 1))

## ----loo-comparison-----------------------------------------------------------
# Compare LOO-CV across biomarkers
cat("\n=== Leave-One-Out Cross-Validation Comparison ===\n\n")

for (biomarker in multi_model$biomarker_names) {
  cat(sprintf("--- %s ---\n", biomarker))
  model <- multi_model$models[[biomarker]]
  print(model$loo)
  cat("\n")
}

## ----session-info-------------------------------------------------------------
sessionInfo()

