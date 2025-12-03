# Conversion from Stan to brms

## Summary

The `seroCOP` package has been successfully converted from using custom Stan code to using the `brms` package. This provides several advantages while maintaining identical statistical functionality.

## Changes Made

### 1. Core Implementation (`R/SeroCOP.R`)

**Before (Stan):**
- Custom Stan model file (`inst/stan/logistic_model.stan`)
- Used `rstan::stan()` to fit models
- Manual parameter extraction with `rstan::extract()`

**After (brms):**
- Non-linear formula using `brms::bf()`:
  ```r
  brms::bf(
    infected ~ ceiling * (inv_logit(-slope * (titre - ec50)) * (1 - floor) + floor),
    floor ~ 1,
    ceiling ~ 1,
    ec50 ~ 1,
    slope ~ 1,
    nl = TRUE
  )
  ```
- Uses `brms::brm()` for model fitting
- Parameter extraction via `brms::as_draws_df()`

### 2. Dependencies (`DESCRIPTION`)

**Added:**
- `brms (>= 2.19.0)` to Imports

**Updated:**
- Description now mentions "brms" instead of "Stan"
- All Stan dependencies (BH, RcppEigen, etc.) remain in LinkingTo since brms uses Stan backend

### 3. Method Updates

#### `fit_model()`
- Now uses `brms::brm()` with non-linear formula
- Priors specified using `brms::set_prior()`
- Handles `cores` argument properly via `do.call()`

#### `predict()`
- Uses `fitted()` function (generic S3 method for brmsfit objects)
- Returns posterior samples directly

#### `predict_protection()`
- Extracts parameters from brms using `brms::as_draws_df()`
- Parameter names: `b_floor_Intercept`, `b_ceiling_Intercept`, etc.

#### `summary()`
- Uses `brms::fixef()` instead of `rstan::summary()`

#### `plot_posteriors()`
- Extracts samples using `brms::as_draws_df()`

### 4. Utility Functions (`R/utils.R`)

**`extract_parameters()`:**
- Updated to use `brms::as_draws_df()` for parameter extraction
- Same output format maintained

### 5. Multi-Biomarker (`R/SeroCOPMulti.R`)

- No code changes needed (automatically uses brms via SeroCOP objects)
- Documentation updated to mention brms

## Advantages of brms

1. **Higher-level interface**: Don't need to write Stan code manually
2. **Automatic Stan code generation**: brms generates optimized Stan code
3. **Better integration**: Works seamlessly with tidyverse and related packages
4. **LOO-CV built-in**: `brms::loo()` works directly on brmsfit objects
5. **Extensibility**: Easy to add random effects, group-level effects, etc.
6. **Same statistical model**: Identical results to custom Stan code

## Mathematical Model

The 4-parameter logistic model remains unchanged:

```
P(infection | titre) = ceiling * (inv_logit(-slope * (titre - ec50)) * (1 - floor) + floor)

P(protection | titre) = 1 - (P(infection) / ceiling)
```

Where:
- **floor**: Proportion of maximum risk remaining at high titre
- **ceiling**: Maximum infection probability at low titre
- **ec50**: Titre at inflection point (50% reduction from ceiling to ceiling×floor)
- **slope**: Steepness of the protective curve

## Testing

All functionality tested and working:
- ✅ Single biomarker models (`SeroCOP`)
- ✅ Multi-biomarker models (`SeroCOPMulti`)
- ✅ Parameter extraction
- ✅ Predictions and protection probabilities
- ✅ LOO-CV
- ✅ Plotting functions
- ✅ Performance metrics (AUC, Brier score)

## Backwards Compatibility

The API remains identical - all existing code using `seroCOP` will work without modification. The only change users will notice is that models now print as "brmsfit" objects instead of "stanfit" objects.

## Next Steps

1. Update vignettes to mention brms (vignette content/results unchanged)
2. Update GitHub Actions to handle brms dependencies
3. Consider adding brms-specific features (e.g., posterior predictive checks with `pp_check()`)
