# SeroCOP R6 Class for Correlates of Protection Analysis

An R6 class for analyzing correlates of protection using Bayesian
methods. Fits a four-parameter logistic model to antibody titre and
infection outcome data.

## Details

This class provides a complete workflow for correlates of protection
analysis:

1.  Data input and validation

2.  Bayesian model fitting using Stan

3.  Model diagnostics and validation

4.  Performance metrics (ROC AUC, Brier score, LOO-CV)

5.  Visualization

## Public fields

- `titre`:

  Numeric vector of antibody titres (log scale recommended)

- `infected`:

  Binary vector of infection status (0/1)

- `fit`:

  Stan fit object (after fitting)

- `loo`:

  LOO-CV object (after fitting)

- `priors`:

  List of prior distributions for model parameters

## Methods

### Public methods

- [`SeroCOP$new()`](#method-SeroCOP-new)

- [`SeroCOP$definePrior()`](#method-SeroCOP-definePrior)

- [`SeroCOP$fit_model()`](#method-SeroCOP-fit_model)

- [`SeroCOP$predict()`](#method-SeroCOP-predict)

- [`SeroCOP$summary()`](#method-SeroCOP-summary)

- [`SeroCOP$get_metrics()`](#method-SeroCOP-get_metrics)

- [`SeroCOP$plot_curve()`](#method-SeroCOP-plot_curve)

- [`SeroCOP$plot_roc()`](#method-SeroCOP-plot_roc)

- [`SeroCOP$plot_parameters()`](#method-SeroCOP-plot_parameters)

- [`SeroCOP$clone()`](#method-SeroCOP-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new SeroCOP object

#### Usage

    SeroCOP$new(titre, infected)

#### Arguments

- `titre`:

  Numeric vector of antibody titres

- `infected`:

  Binary vector (0/1) of infection outcomes

#### Returns

A new SeroCOP object

------------------------------------------------------------------------

### Method `definePrior()`

Define custom prior distributions for model parameters

#### Usage

    SeroCOP$definePrior(
      floor_alpha = NULL,
      floor_beta = NULL,
      ceiling_alpha = NULL,
      ceiling_beta = NULL,
      ec50_mean = NULL,
      ec50_sd = NULL,
      slope_mean = NULL,
      slope_sd = NULL
    )

#### Arguments

- `floor_alpha`:

  Alpha parameter for floor beta prior (default: 1)

- `floor_beta`:

  Beta parameter for floor beta prior (default: 9)

- `ceiling_alpha`:

  Alpha parameter for ceiling beta prior (default: 9)

- `ceiling_beta`:

  Beta parameter for ceiling beta prior (default: 1)

- `ec50_mean`:

  Mean for ec50 normal prior (default: midpoint of titre range)

- `ec50_sd`:

  SD for ec50 normal prior (default: titre range / 4)

- `slope_mean`:

  Mean for slope normal prior (default: 0)

- `slope_sd`:

  SD for slope normal prior (default: 2)

#### Returns

Self (invisibly)

#### Examples

    \dontrun{
    model <- SeroCOP$new(titre = titre, infected = infected)

    # Use default priors centered on data
    model$definePrior()

    # Custom priors
    model$definePrior(
      floor_alpha = 2, floor_beta = 18,    # Stronger prior for low floor
      ceiling_alpha = 18, ceiling_beta = 2, # Stronger prior for high ceiling
      ec50_mean = 2.0, ec50_sd = 1.0,      # Custom ec50 prior
      slope_mean = 0, slope_sd = 1         # More conservative slope
    )
    }

------------------------------------------------------------------------

### Method `fit_model()`

Fit the Bayesian logistic model

#### Usage

    SeroCOP$fit_model(
      chains = 4,
      iter = 2000,
      warmup = floor(iter/2),
      thin = 1,
      ...
    )

#### Arguments

- `chains`:

  Number of MCMC chains (default: 4)

- `iter`:

  Number of iterations per chain (default: 2000)

- `warmup`:

  Number of warmup iterations (default: iter/2)

- `thin`:

  Thinning interval (default: 1)

- `...`:

  Additional arguments passed to rstan::sampling

#### Returns

Self (invisibly)

------------------------------------------------------------------------

### Method [`predict()`](https://rdrr.io/r/stats/predict.html)

Get posterior predictions for infection probability

#### Usage

    SeroCOP$predict(newdata = NULL)

#### Arguments

- `newdata`:

  Optional new titre values for prediction

#### Returns

Matrix of posterior predictions (iterations x observations)

------------------------------------------------------------------------

### Method [`summary()`](https://rdrr.io/r/base/summary.html)

Get summary statistics for model parameters

#### Usage

    SeroCOP$summary()

#### Returns

Data frame with parameter summaries

------------------------------------------------------------------------

### Method `get_metrics()`

Calculate performance metrics

#### Usage

    SeroCOP$get_metrics()

#### Returns

List containing ROC AUC, Brier score, and LOO-CV metrics

------------------------------------------------------------------------

### Method `plot_curve()`

Plot the fitted curve with uncertainty

#### Usage

    SeroCOP$plot_curve(title = "Correlates of Risk Curve", ...)

#### Arguments

- `title`:

  Plot title

- `...`:

  Additional arguments passed to ggplot2

#### Returns

A ggplot object

------------------------------------------------------------------------

### Method `plot_roc()`

Plot ROC curve

#### Usage

    SeroCOP$plot_roc(title = "ROC Curve")

#### Arguments

- `title`:

  Plot title

#### Returns

A ggplot object

------------------------------------------------------------------------

### Method `plot_parameters()`

Plot posterior distributions of parameters

#### Usage

    SeroCOP$plot_parameters()

#### Returns

A ggplot object

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SeroCOP$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
if (FALSE) { # \dontrun{
# Create synthetic data
set.seed(123)
n <- 200
titre <- rnorm(n, mean = 2, sd = 1.5)
prob <- 0.05 + 0.9 / (1 + exp(2 * (titre - 1.5)))
infected <- rbinom(n, 1, prob)

# Initialize and fit
model <- SeroCOP$new(titre = titre, infected = infected)
model$fit(chains = 4, iter = 2000)

# Get metrics
model$get_metrics()

# Plot results
model$plot_curve()
model$plot_roc()
} # }

## ------------------------------------------------
## Method `SeroCOP$definePrior`
## ------------------------------------------------

if (FALSE) { # \dontrun{
model <- SeroCOP$new(titre = titre, infected = infected)

# Use default priors centered on data
model$definePrior()

# Custom priors
model$definePrior(
  floor_alpha = 2, floor_beta = 18,    # Stronger prior for low floor
  ceiling_alpha = 18, ceiling_beta = 2, # Stronger prior for high ceiling
  ec50_mean = 2.0, ec50_sd = 1.0,      # Custom ec50 prior
  slope_mean = 0, slope_sd = 1         # More conservative slope
)
} # }
```
