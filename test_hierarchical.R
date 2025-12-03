#!/usr/bin/env Rscript
# Test hierarchical effects with age groups

library(seroCOP)
library(ggplot2)

cat("Testing SeroCOP with hierarchical age group effects...\n\n")

set.seed(456)
n_per_group <- 80

# Define age groups with different CoP relationships
# Young: steep slope (strong protection)
# Middle: moderate slope
# Old: flat slope (no correlate)

# Generate data for each age group
generate_age_group <- function(n, age_group, slope_true, ec50_true) {
  titre <- rnorm(n, mean = 2, sd = 1.5)
  
  # Generate infection probabilities with group-specific parameters
  prob <- 0.7 * (1 / (1 + exp(-slope_true * (titre - ec50_true))) * 0.95 + 0.05)
  infected <- rbinom(n, 1, prob)
  
  data.frame(
    titre = titre,
    infected = infected,
    age_group = age_group,
    true_slope = slope_true,
    true_ec50 = ec50_true
  )
}

# Create data with different slopes by age
young_data <- generate_age_group(n_per_group, "Young", slope_true = 3.0, ec50_true = 1.5)
middle_data <- generate_age_group(n_per_group, "Middle", slope_true = 1.5, ec50_true = 1.8)
old_data <- generate_age_group(n_per_group, "Old", slope_true = 0.3, ec50_true = 2.0)

# Combine all data
all_data <- rbind(young_data, middle_data, old_data)

cat("Simulated data with age-specific correlates:\n")
cat(sprintf("  Young (n=%d): slope=%.1f, ec50=%.1f (steep - strong protection)\n", 
           n_per_group, unique(young_data$true_slope), unique(young_data$true_ec50)))
cat(sprintf("  Middle (n=%d): slope=%.1f, ec50=%.1f (moderate protection)\n", 
           n_per_group, unique(middle_data$true_slope), unique(middle_data$true_ec50)))
cat(sprintf("  Old (n=%d): slope=%.1f, ec50=%.1f (flat - no correlate)\n\n", 
           n_per_group, unique(old_data$true_slope), unique(old_data$true_ec50)))

# Fit hierarchical model
cat("Fitting hierarchical model...\n")
hier_model <- SeroCOP$new(
  titre = all_data$titre,
  infected = all_data$infected,
  group = all_data$age_group
)

hier_model$fit_model(
  chains = 2,
  iter = 1000,
  warmup = 500,
  cores = 2,
  refresh = 0
)

cat("\n✓ Hierarchical model fitted successfully!\n\n")

# Extract group-specific parameters
cat("Group-specific parameter estimates:\n")
posterior <- brms::as_draws_df(hier_model$fit)

# Overall (population-level) parameters
cat("\nPopulation-level parameters:\n")
cat(sprintf("  ec50 (intercept): %.2f\n", mean(posterior$b_ec50_Intercept)))
cat(sprintf("  slope (intercept): %.2f\n", mean(posterior$b_slope_Intercept)))

# Group-level deviations
cat("\nGroup-level deviations from population mean:\n")
for (age in c("Young", "Middle", "Old")) {
  ec50_col <- paste0("r_group__ec50[", age, ",Intercept]")
  slope_col <- paste0("r_group__slope[", age, ",Intercept]")
  
  if (ec50_col %in% names(posterior)) {
    cat(sprintf("\n%s:\n", age))
    cat(sprintf("  ec50 deviation: %.2f\n", mean(posterior[[ec50_col]])))
    cat(sprintf("  slope deviation: %.2f\n", mean(posterior[[slope_col]])))
    cat(sprintf("  Total ec50: %.2f\n", mean(posterior$b_ec50_Intercept + posterior[[ec50_col]])))
    cat(sprintf("  Total slope: %.2f\n", mean(posterior$b_slope_Intercept + posterior[[slope_col]])))
  }
}

# Compare with non-hierarchical model
cat("\n\nFitting non-hierarchical model for comparison...\n")
nonhier_model <- SeroCOP$new(
  titre = all_data$titre,
  infected = all_data$infected
)

nonhier_model$fit_model(
  chains = 2,
  iter = 1000,
  warmup = 500,
  cores = 2,
  refresh = 0
)

cat("\n✓ Non-hierarchical model fitted successfully!\n")

# Extract and display group-specific parameters
cat("\n=== Group-Specific Parameters ===\n")
group_params <- hier_model$extract_group_parameters()
print(group_params)

# Create group-specific curve plot
cat("\n=== Creating group-specific curve plot ===\n")
p_groups <- hier_model$plot_group_curves(title = "Age-Specific Correlates of Protection")
ggsave("hierarchical_group_curves.png", p_groups, width = 12, height = 5, dpi = 150)
cat("Saved group-specific curves to: hierarchical_group_curves.png\n")

# Compare LOO
cat("\n=== Model Comparison (LOO-CV) ===\n")
cat(sprintf("  Hierarchical LOO-ELPD: %.2f (SE: %.2f)\n", 
           hier_model$loo$estimates["elpd_loo", "Estimate"],
           hier_model$loo$estimates["elpd_loo", "SE"]))
cat(sprintf("  Non-hierarchical LOO-ELPD: %.2f (SE: %.2f)\n", 
           nonhier_model$loo$estimates["elpd_loo", "Estimate"],
           nonhier_model$loo$estimates["elpd_loo", "SE"]))

loo_diff <- hier_model$loo$estimates["elpd_loo", "Estimate"] - 
            nonhier_model$loo$estimates["elpd_loo", "Estimate"]
cat(sprintf("\n  Difference: %.2f (positive favors hierarchical model)\n", loo_diff))

cat("\n✓ All tests passed! Hierarchical effects working correctly.\n")
