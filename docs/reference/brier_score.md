# Calculate Brier score

Calculate the Brier score for probabilistic predictions

## Usage

``` r
brier_score(observed, predicted)
```

## Arguments

- observed:

  Binary vector of observed outcomes (0/1)

- predicted:

  Vector of predicted probabilities

## Value

Numeric Brier score (lower is better)

## Examples

``` r
observed <- c(0, 1, 1, 0, 1)
predicted <- c(0.1, 0.8, 0.7, 0.2, 0.9)
brier_score(observed, predicted)
#> [1] 0.038
```
