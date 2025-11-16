# SeroCOPMulti R6 Class for Multi-Biomarker Correlates of Protection Analysis

An R6 class for analyzing multiple biomarkers simultaneously. Fits
separate models for each biomarker and provides comparison tools.

## Public fields

- `titre`:

  Matrix of antibody titres (columns = biomarkers)

- `infected`:

  Binary vector of infection status (0/1)

- `biomarker_names`:

  Names of biomarkers

- `models`:

  List of fitted SeroCOP objects

## Methods

### Public methods

- [`SeroCOPMulti$new()`](#method-SeroCOPMulti-new)

- [`SeroCOPMulti$fit_all()`](#method-SeroCOPMulti-fit_all)

- [`SeroCOPMulti$compare_biomarkers()`](#method-SeroCOPMulti-compare_biomarkers)

- [`SeroCOPMulti$plot_comparison()`](#method-SeroCOPMulti-plot_comparison)

- [`SeroCOPMulti$plot_all_curves()`](#method-SeroCOPMulti-plot_all_curves)

- [`SeroCOPMulti$clone()`](#method-SeroCOPMulti-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new SeroCOPMulti object

#### Usage

    SeroCOPMulti$new(titre, infected, biomarker_names = NULL)

#### Arguments

- `titre`:

  Matrix of antibody titres (rows = samples, cols = biomarkers)

- `infected`:

  Binary vector (0/1) of infection outcomes

- `biomarker_names`:

  Optional vector of biomarker names

#### Returns

A new SeroCOPMulti object

------------------------------------------------------------------------

### Method `fit_all()`

Fit models for all biomarkers

#### Usage

    SeroCOPMulti$fit_all(
      chains = 4,
      iter = 2000,
      warmup = floor(iter/2),
      cores = 1,
      ...
    )

#### Arguments

- `chains`:

  Number of MCMC chains (default: 4)

- `iter`:

  Number of iterations per chain (default: 2000)

- `warmup`:

  Number of warmup iterations (default: iter/2)

- `cores`:

  Number of cores for parallel processing (default: 1)

- `...`:

  Additional arguments passed to rstan::sampling

#### Returns

Self (invisibly)

------------------------------------------------------------------------

### Method `compare_biomarkers()`

Compare biomarkers on AUC and LOO metrics

#### Usage

    SeroCOPMulti$compare_biomarkers()

#### Returns

Data frame with comparison metrics

------------------------------------------------------------------------

### Method `plot_comparison()`

Plot biomarkers on AUC vs LOO plane

#### Usage

    SeroCOPMulti$plot_comparison(add_labels = TRUE)

#### Arguments

- `add_labels`:

  Logical, add biomarker labels (default: TRUE)

#### Returns

A ggplot object

------------------------------------------------------------------------

### Method `plot_all_curves()`

Plot all fitted curves

#### Usage

    SeroCOPMulti$plot_all_curves()

#### Returns

A ggplot object

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SeroCOPMulti$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
if (FALSE) { # \dontrun{
# Create multi-biomarker data
titres <- matrix(rnorm(600), ncol = 3)
colnames(titres) <- c("IgG", "IgA", "Neutralization")
infected <- rbinom(200, 1, 0.3)

# Initialize and fit
multi_model <- SeroCOPMulti$new(titre = titres, infected = infected)
multi_model$fit_all(chains = 4, iter = 2000)

# Compare biomarkers
multi_model$compare_biomarkers()
multi_model$plot_comparison()
} # }
```
