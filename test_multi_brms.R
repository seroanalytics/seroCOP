#!/usr/bin/env Rscript
# Test SeroCOPMulti with brms implementation

library(seroCOP)

cat("Testing SeroCOPMulti with brms...\n\n")

# Generate multi-biomarker data
set.seed(123)
n <- 80

# Create 3 biomarkers with different properties
titre1 <- rnorm(n, mean = 2, sd = 1.5)
titre2 <- rnorm(n, mean = 3, sd = 1.2)
titre3 <- rnorm(n, mean = 1.5, sd = 1.8)

titres <- cbind(IgG = titre1, IgA = titre2, Neutralizing = titre3)

# Generate infection outcomes based on first biomarker
prob <- 0.7 * (1 / (1 + exp(-2 * (titre1 - 1.5))) * 0.95 + 0.05)
infected <- rbinom(n, 1, prob)

cat("Created synthetic data:\n")
cat(sprintf("  %d observations\n", n))
cat(sprintf("  %d biomarkers: %s\n", ncol(titres), paste(colnames(titres), collapse=", ")))
cat(sprintf("  Infection rate: %.1f%%\n\n", mean(infected) * 100))

# Create multi-biomarker model
multi_model <- SeroCOPMulti$new(
  titre = titres,
  infected = infected
)

# Fit all models (small chains for speed)
cat("\nFitting models...\n")
multi_model$fit_all(
  chains = 2,
  iter = 500,
  warmup = 250,
  cores = 2,
  refresh = 0
)

# Compare biomarkers
cat("\nComparing biomarkers...\n")
comparison <- multi_model$compare_biomarkers()
print(comparison)

cat("\n✓ SeroCOPMulti with brms working correctly!\n")
