# Simulate data for correlates of protection analysis

Generate synthetic data from a four-parameter logistic model for testing
and validation purposes.

## Usage

``` r
simulate_serocop_data(
  n = 200,
  floor = 0.05,
  ceiling = 0.9,
  ec50 = 1.5,
  slope = 2,
  titre_mean = 2,
  titre_sd = 1.5,
  seed = NULL
)
```

## Arguments

- n:

  Number of observations

- floor:

  Lower asymptote (minimum infection probability)

- ceiling:

  Upper asymptote (maximum infection probability)

- ec50:

  Titre at 50% between floor and ceiling

- slope:

  Steepness of the curve

- titre_mean:

  Mean of titre distribution (normal)

- titre_sd:

  Standard deviation of titre distribution

- seed:

  Random seed for reproducibility

## Value

A list containing:

- titre:

  Vector of simulated titres

- infected:

  Vector of simulated infection outcomes

- prob_true:

  True infection probabilities

- params:

  List of true parameters used for simulation

## Examples

``` r
# Simulate data with default parameters
data <- simulate_serocop_data(n = 200, seed = 123)

# Custom parameters
data <- simulate_serocop_data(
  n = 500,
  floor = 0.02,
  ceiling = 0.95,
  ec50 = 2.0,
  slope = 3.0,
  seed = 456
)
```
