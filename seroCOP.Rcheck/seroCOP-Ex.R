pkgname <- "seroCOP"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('seroCOP')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("SeroCOP")
### * SeroCOP

flush(stderr()); flush(stdout())

### Name: SeroCOP
### Title: SeroCOP R6 Class
### Aliases: SeroCOP

### ** Examples

## Not run: 
##D # Create synthetic data
##D set.seed(123)
##D n <- 200
##D titre <- rnorm(n, mean = 2, sd = 1.5)
##D prob <- 0.05 + 0.9 / (1 + exp(2 * (titre - 1.5)))
##D infected <- rbinom(n, 1, prob)
##D 
##D # Initialize and fit
##D model <- SeroCOP$new(titre = titre, infected = infected)
##D model$fit(chains = 4, iter = 2000)
##D 
##D # Get metrics
##D model$get_metrics()
##D 
##D # Plot results
##D model$plot_curve()
##D model$plot_roc()
## End(Not run)

## ------------------------------------------------
## Method `SeroCOP$definePrior`
## ------------------------------------------------

## Not run: 
##D model <- SeroCOP$new(titre = titre, infected = infected)
##D 
##D # Use default priors centered on data
##D model$definePrior()
##D 
##D # Custom priors
##D model$definePrior(
##D   floor_alpha = 2, floor_beta = 18,    # Stronger prior for low floor
##D   ceiling_alpha = 18, ceiling_beta = 2, # Stronger prior for high ceiling
##D   ec50_mean = 2.0, ec50_sd = 1.0,      # Custom ec50 prior
##D   slope_mean = 0, slope_sd = 1         # More conservative slope
##D )
## End(Not run)



cleanEx()
nameEx("SeroCOPMulti")
### * SeroCOPMulti

flush(stderr()); flush(stdout())

### Name: SeroCOPMulti
### Title: SeroCOPMulti R6 Class for Multi-Biomarker Correlates of
###   Protection Analysis
### Aliases: SeroCOPMulti

### ** Examples

## Not run: 
##D # Create multi-biomarker data
##D titres <- matrix(rnorm(600), ncol = 3)
##D colnames(titres) <- c("IgG", "IgA", "Neutralization")
##D infected <- rbinom(200, 1, 0.3)
##D 
##D # Initialize and fit
##D multi_model <- SeroCOPMulti$new(titre = titres, infected = infected)
##D multi_model$fit_all(chains = 4, iter = 2000)
##D 
##D # Compare biomarkers
##D multi_model$compare_biomarkers()
##D multi_model$plot_comparison()
##D 
##D # With hierarchical effects
##D age_group <- sample(c("Young", "Middle", "Old"), 200, replace = TRUE)
##D hier_model <- SeroCOPMulti$new(titre = titres, infected = infected, group = age_group)
##D hier_model$fit_all(chains = 4, iter = 2000)
##D hier_model$plot_group_curves()
##D hier_model$extract_group_parameters()
## End(Not run)



cleanEx()
nameEx("brier_score")
### * brier_score

flush(stderr()); flush(stdout())

### Name: brier_score
### Title: Calculate Brier score
### Aliases: brier_score

### ** Examples

observed <- c(0, 1, 1, 0, 1)
predicted <- c(0.1, 0.8, 0.7, 0.2, 0.9)
brier_score(observed, predicted)



cleanEx()
nameEx("simulate_serocop_data")
### * simulate_serocop_data

flush(stderr()); flush(stdout())

### Name: simulate_serocop_data
### Title: Simulate data for correlates of protection analysis
### Aliases: simulate_serocop_data

### ** Examples

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



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
