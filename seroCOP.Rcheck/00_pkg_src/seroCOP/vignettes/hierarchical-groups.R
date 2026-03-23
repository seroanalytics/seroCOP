## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 10,
  fig.height = 6,
  warning = FALSE,
  message = FALSE
)

## ----setup--------------------------------------------------------------------
library(seroCOP)
library(ggplot2)
set.seed(456)

## ----simulate-data------------------------------------------------------------
n_per_group <- 80  # Sample size per age group

# Function to simulate data for one age group
simulate_group <- function(n, group_name, slope, ec50) {
  # Generate antibody titres
  titre <- rnorm(n, mean = 2.5, sd = 2.0)
  
  # Calculate infection probability using the logistic model
  # floor = 0.05, ceiling = 0.70
  floor <- 0.05
  ceiling <- 0.70
  
  logit_part <- 1 / (1 + exp(slope * (titre - ec50)))
  prob_infection <- ceiling * (logit_part * (1 - floor) + floor)
  
  # Generate infection outcomes
  infected <- rbinom(n, 1, prob_infection)
  
  return(data.frame(
    titre = titre,
    infected = infected,
    age_group = group_name
  ))
}

# Simulate three age groups with different CoP strengths
# Note: Higher slope values mean steeper curves (stronger correlate)
young <- simulate_group(n_per_group, "Young", slope = 3.0, ec50 = 1.5)
middle <- simulate_group(n_per_group, "Middle", slope = 1.5, ec50 = 1.8)
old <- simulate_group(n_per_group, "Old", slope = 0.3, ec50 = 2.0)

# Combine all groups
all_data <- rbind(young, middle, old)

cat(sprintf("Simulated data with age-specific correlates:\n"))
cat(sprintf("  Young (n=%d): slope=3.0, ec50=1.5 (steep - strong protection)\n", n_per_group))
cat(sprintf("  Middle (n=%d): slope=1.5, ec50=1.8 (moderate protection)\n", n_per_group))
cat(sprintf("  Old (n=%d): slope=0.3, ec50=2.0 (flat - weak correlate)\n", n_per_group))
cat(sprintf("\nTotal sample size: %d\n", nrow(all_data)))
cat(sprintf("Overall infection rate: %.1f%%\n", 100 * mean(all_data$infected)))

## ----plot-raw-data------------------------------------------------------------
ggplot(all_data, aes(x = titre, y = infected, color = age_group)) +
  geom_point(alpha = 0.4, position = position_jitter(height = 0.02)) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1) +
  facet_wrap(~age_group) +
  labs(
    title = "Infection vs Titre by Age Group",
    x = "Antibody Titre (log scale)",
    y = "Infected (0/1)",
    color = "Age Group"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

## ----fit-hierarchical, results='hide'-----------------------------------------
# Create SeroCOP object with group variable
hier_model <- SeroCOP$new(
  titre = all_data$titre,
  infected = all_data$infected,
  group = all_data$age_group  # Add group variable for hierarchical modeling
)

# Fit the hierarchical model
hier_model$fit_model(
  chains = 4,
  iter = 2000,
  warmup = 1000,
  cores = 4,
  refresh = 0
)

## ----show-hierarchical-fit----------------------------------------------------
cat("✓ Hierarchical model fitted successfully!\n")

## ----extract-parameters-------------------------------------------------------
# Extract group-specific estimates
group_params <- hier_model$extract_group_parameters()
print(group_params)

## ----plot-group-curves, fig.width=12, fig.height=5----------------------------
hier_model$plot_group_curves(title = "Age-Specific Correlates of Protection")

## ----fit-pooled, results='hide'-----------------------------------------------
# Fit without group variable (pooled model)
pooled_model <- SeroCOP$new(
  titre = all_data$titre,
  infected = all_data$infected
)

pooled_model$fit_model(
  chains = 4,
  iter = 2000,
  warmup = 1000,
  cores = 4,
  refresh = 0
)

## ----show-pooled-fit----------------------------------------------------------
cat("✓ Non-hierarchical (pooled) model fitted successfully!\n")

## ----compare-models-----------------------------------------------------------
# Extract LOO estimates
hier_loo <- hier_model$loo$estimates["elpd_loo", c("Estimate", "SE")]
pooled_loo <- pooled_model$loo$estimates["elpd_loo", c("Estimate", "SE")]

# Calculate difference
loo_diff <- hier_loo["Estimate"] - pooled_loo["Estimate"]
loo_se <- sqrt(hier_loo["SE"]^2 + pooled_loo["SE"]^2)

cat("\n=== Model Comparison (LOO-CV) ===\n")
cat(sprintf("Hierarchical model ELPD: %.2f (SE: %.2f)\n", 
           hier_loo["Estimate"], hier_loo["SE"]))
cat(sprintf("Pooled model ELPD:       %.2f (SE: %.2f)\n", 
           pooled_loo["Estimate"], pooled_loo["SE"]))
cat(sprintf("\nDifference: %.2f (SE: %.2f)\n", loo_diff, loo_se))
cat(sprintf("Z-score: %.2f\n", loo_diff / loo_se))

if (loo_diff > 2 * loo_se) {
  cat("\n✓ Strong evidence favoring hierarchical model\n")
} else if (loo_diff > 0) {
  cat("\n→ Hierarchical model preferred but evidence is weak\n")
} else {
  cat("\n→ No clear advantage for hierarchical model\n")
}

