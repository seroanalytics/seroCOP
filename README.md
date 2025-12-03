# seroCOP

<!-- badges: start -->
<!-- badges: end -->

**seroCOP** is an R package for analysing correlates of protection using Bayesian methods. It fits generalized four-parameter logistic functions to antibody titre and infection outcome data using Stan, providing comprehensive model diagnostics and performance metrics.

## Installation

You can install the development version of seroCOP from GitHub:

```r
# Install devtools if needed
install.packages("devtools")

# Install seroCOP
devtools::install_github("seroanalytics/seroCOP")
```

**Note**: This package requires a working installation of brms and Stan. The brms package will handle Stan installation automatically when you first use it.

## Quick Start

```r
library(seroCOP)

# Simulate data
sim_data <- simulate_serocop_data(
  n = 200,
  floor = 0.05,
  ceiling = 0.90,
  ec50 = 1.5,
  slope = 2.0,
  seed = 123
)

# Create SeroCOP object
model <- SeroCOP$new(
  titre = sim_data$titre,
  infected = sim_data$infected
)

# Fit the model (uses brms with sensible defaults)
model$fit_model(chains = 4, iter = 2000)

# Get performance metrics
metrics <- model$get_metrics()

# Visualize results
model$plot_curve()
model$plot_roc()

# Extract parameter estimates
model$summary()
```

## Hierarchical Group Effects

Model group-level heterogeneity (e.g., age groups with different correlates):

```r
# Simulate age-stratified data
age_group <- sample(c("Young", "Middle", "Old"), 200, replace = TRUE)

# Fit hierarchical model
hier_model <- SeroCOP$new(
  titre = sim_data$titre,
  infected = sim_data$infected,
  group = age_group  # Adds random effects on slope and ec50
)

hier_model$fit_model(chains = 4, iter = 2000)

# Extract group-specific parameters
group_params <- hier_model$extract_group_parameters()

# Plot group-specific curves
hier_model$plot_group_curves()
```

## Multi-Biomarker Analysis

Compare multiple biomarkers simultaneously with optional hierarchical effects:

```r
# Prepare multi-biomarker data (matrix with columns = biomarkers)
titre_matrix <- cbind(
  IgG = rnorm(200, 2, 1),
  IgA = rnorm(200, 1.5, 1),
  Neutralization = rnorm(200, 3, 1.2)
)

# Initialize and fit (without groups)
multi_model <- SeroCOPMulti$new(
  titre = titre_matrix,
  infected = infected
)
multi_model$fit_all(chains = 4, iter = 2000)

# Compare biomarkers
multi_model$compare_biomarkers()
multi_model$plot_comparison()  # AUC vs LOO-ELPD plot
multi_model$plot_all_curves()

# With hierarchical effects across biomarkers
hier_multi <- SeroCOPMulti$new(
  titre = titre_matrix,
  infected = infected,
  group = age_group
)
hier_multi$fit_all(chains = 4, iter = 2000)
hier_multi$plot_group_curves()  # Shows all biomarkers × groups
hier_multi$extract_group_parameters()
```

## Model

The package fits a four-parameter logistic model using brms (Bayesian Regression Models using Stan):

$$P(\text{infection} | \text{titre}) = \text{ceiling} \times \left[\frac{1}{1 + e^{\text{slope} \times (\text{titre} - \text{ec50})}} \times (1 - \text{floor}) + \text{floor}\right]$$

Where:
- **floor**: Lower asymptote (baseline infection probability at high titres)
- **ceiling**: Upper asymptote (maximum infection probability at low titres)  
- **ec50**: Titre at 50% between floor and ceiling (inflection point)
- **slope**: Steepness of the dose-response curve

### Hierarchical Extension

When group variables are provided, the model adds random intercepts:

$$\text{ec50}_g \sim \mathcal{N}(\text{ec50}_{\text{population}}, \sigma_{\text{ec50}})$$
$$\text{slope}_g \sim \mathcal{N}(\text{slope}_{\text{population}}, \sigma_{\text{slope}})$$

This allows for group-specific parameters while borrowing strength across groups.

## Priors

The package automatically sets weakly informative priors based on your data:

- **floor** ~ Beta(1, 9): Favors low baseline infection (mean ≈ 0.1)
- **ceiling** ~ Beta(9, 1): Favors high maximum infection (mean ≈ 0.9)
- **ec50** ~ Normal(μ, σ): Centered at midpoint of titre range
- **slope** ~ Normal(0, 2): Weakly informative, truncated at 0
- **Group-level SDs** ~ student_t(3, 0, 2.5): Weakly informative for hierarchical models


## Documentation

See the package vignettes for detailed examples:

- **SeroCOP Overview**: Basic usage and model details
- **Simulation and Recovery**: Parameter recovery validation
- **Multi-Biomarker Analysis**: Comparing multiple antibody types
- **Hierarchical Group Effects**: Age-stratified and group-specific modeling

```r
# View available vignettes
browseVignettes("seroCOP")
```

## Performance Metrics

The package calculates several performance metrics:

- **ROC AUC**: Area under the ROC curve (discrimination)
- **Brier Score**: Mean squared error of predictions (calibration)
- **LOO-CV**: Leave-one-out cross-validation (out-of-sample performance)
- **ELPD**: Expected log pointwise predictive density

## Citation

If you use this package in your research, please cite:

```
Hodgson, D. (2025). seroCOP: Correlates of Protection Analysis Using Bayesian Methods.
R package version 0.1.0.
```

## Development

This package is under active development. Contributions and feedback are welcome.


## Support

For questions, issues, or feature requests, please open an issue on GitHub.
