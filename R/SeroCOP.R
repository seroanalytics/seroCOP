#' SeroCOP R6 Class for Correlates of Protection Analysis
#'
#' @description
#' An R6 class for analyzing correlates of protection using Bayesian methods.
#' Fits a four-parameter logistic model to antibody titre and infection outcome data.
#'
#' @details
#' This class provides a complete workflow for correlates of protection analysis:
#' 1. Data input and validation
#' 2. Bayesian model fitting using brms (which uses Stan backend)
#' 3. Model diagnostics and validation
#' 4. Performance metrics (ROC AUC, Brier score, LOO-CV)
#' 5. Visualization
#'
#' @concept r6-classes
#' @export
#' @examples
#' \dontrun{
#' # Create synthetic data
#' set.seed(123)
#' n <- 200
#' titre <- rnorm(n, mean = 2, sd = 1.5)
#' prob <- 0.05 + 0.9 / (1 + exp(2 * (titre - 1.5)))
#' infected <- rbinom(n, 1, prob)
#'
#' # Initialize and fit
#' model <- SeroCOP$new(titre = titre, infected = infected)
#' model$fit(chains = 4, iter = 2000)
#' 
#' # Get metrics
#' model$get_metrics()
#' 
#' # Plot results
#' model$plot_curve()
#' model$plot_roc()
#' }
#' @title SeroCOP R6 Class
#' @description This class implements a Bayesian logistic model for analyzing correlates of protection.
#' It provides methods for model fitting, prediction, performance evaluation, and visualization.
SeroCOP <- R6::R6Class(
  "SeroCOP",
  
  public = list(
    #' @field titre Numeric vector of antibody titres (log scale recommended)
    titre = NULL,
    
    #' @field infected Binary vector of infection status (0/1)
    infected = NULL,
    
    #' @field fit Stan fit object (after fitting)
    fit = NULL,
    
    #' @field loo LOO-CV object (after fitting)
    loo = NULL,
    
    #' @field priors List of prior distributions for model parameters
    priors = NULL,
    
    #' @description
    #' Create a new SeroCOP object
    #' @param titre Numeric vector of antibody titres
    #' @param infected Binary vector (0/1) of infection outcomes
    #' @return A new SeroCOP object
    initialize = function(titre, infected) {
      # Validate inputs
      if (!is.numeric(titre)) {
        stop("titre must be numeric")
      }
      if (!all(infected %in% c(0, 1))) {
        stop("infected must be binary (0/1)")
      }
      if (length(titre) != length(infected)) {
        stop("titre and infected must have the same length")
      }
      if (any(is.na(titre)) || any(is.na(infected))) {
        stop("Missing values are not allowed")
      }

      self$titre <- titre
      self$infected <- infected

      # Set default priors
      self$priors <- private$get_default_priors()

      message(sprintf("SeroCOP initialized with %d observations", length(titre)))
      message(sprintf("  Infection rate: %.1f%%", mean(infected) * 100))
      message(sprintf("  Titre range: [%.2f, %.2f]", min(titre), max(titre)))
    },
    
    #' @description
    #' Define custom prior distributions for model parameters
    #' @param floor_alpha Alpha parameter for floor beta prior (default: 1)
    #' @param floor_beta Beta parameter for floor beta prior (default: 9)
    #' @param ceiling_alpha Alpha parameter for ceiling beta prior (default: 9)
    #' @param ceiling_beta Beta parameter for ceiling beta prior (default: 1)
    #' @param ec50_mean Mean for ec50 normal prior (default: midpoint of titre range)
    #' @param ec50_sd SD for ec50 normal prior (default: titre range / 4)
    #' @param slope_mean Mean for slope normal prior (default: 0)
    #' @param slope_sd SD for slope normal prior (default: 2)
    #' @return Self (invisibly)
    #' @examples
    #' \dontrun{
    #' model <- SeroCOP$new(titre = titre, infected = infected)
    #' 
    #' # Use default priors centered on data
    #' model$definePrior()
    #' 
    #' # Custom priors
    #' model$definePrior(
    #'   floor_alpha = 2, floor_beta = 18,    # Stronger prior for low floor
    #'   ceiling_alpha = 18, ceiling_beta = 2, # Stronger prior for high ceiling
    #'   ec50_mean = 2.0, ec50_sd = 1.0,      # Custom ec50 prior
    #'   slope_mean = 0, slope_sd = 1         # More conservative slope
    #' )
    #' }
    definePrior = function(floor_alpha = NULL, floor_beta = NULL,
                          ceiling_alpha = NULL, ceiling_beta = NULL,
                          ec50_mean = NULL, ec50_sd = NULL,
                          slope_mean = NULL, slope_sd = NULL) {
      
      # Get defaults
      defaults <- private$get_default_priors()
      
      # Update with user-specified values
      self$priors <- list(
        floor_alpha = if (!is.null(floor_alpha)) floor_alpha else defaults$floor_alpha,
        floor_beta = if (!is.null(floor_beta)) floor_beta else defaults$floor_beta,
        ceiling_alpha = if (!is.null(ceiling_alpha)) ceiling_alpha else defaults$ceiling_alpha,
        ceiling_beta = if (!is.null(ceiling_beta)) ceiling_beta else defaults$ceiling_beta,
        ec50_mean = if (!is.null(ec50_mean)) ec50_mean else defaults$ec50_mean,
        ec50_sd = if (!is.null(ec50_sd)) ec50_sd else defaults$ec50_sd,
        slope_mean = if (!is.null(slope_mean)) slope_mean else defaults$slope_mean,
        slope_sd = if (!is.null(slope_sd)) slope_sd else defaults$slope_sd
      )
      
      message("Prior distributions updated:")
      message(sprintf("  floor ~ Beta(%.1f, %.1f)", self$priors$floor_alpha, self$priors$floor_beta))
      message(sprintf("  ceiling ~ Beta(%.1f, %.1f)", self$priors$ceiling_alpha, self$priors$ceiling_beta))
      message(sprintf("  ec50 ~ Normal(%.2f, %.2f)", self$priors$ec50_mean, self$priors$ec50_sd))
      message(sprintf("  slope ~ Normal(%.2f, %.2f) [truncated at 0]", self$priors$slope_mean, self$priors$slope_sd))
      
      invisible(self)
    },
    
    #' @description
    #' Fit the Bayesian logistic model
    #' @param chains Number of MCMC chains (default: 4)
    #' @param iter Number of iterations per chain (default: 2000)
    #' @param warmup Number of warmup iterations (default: iter/2)
    #' @param thin Thinning interval (default: 1)
    #' @param ... Additional arguments passed to brms::brm
    #' @return Self (invisibly)
    fit_model = function(chains = 4, iter = 2000, warmup = 1000, thin = 1, ...) {
      message("Fitting Bayesian logistic model using brms...")
      
      # Prepare data frame
      model_data <- data.frame(
        infected = self$infected,
        titre = self$titre
      )
      
      message("\nUsing prior distributions:")
      message(sprintf("  floor ~ Beta(%.1f, %.1f)", self$priors$floor_alpha, self$priors$floor_beta))
      message(sprintf("  ceiling ~ Beta(%.1f, %.1f)", self$priors$ceiling_alpha, self$priors$ceiling_beta))
      message(sprintf("  ec50 ~ Normal(%.2f, %.2f)", self$priors$ec50_mean, self$priors$ec50_sd))
      message(sprintf("  slope ~ Normal(%.2f, %.2f)\n", self$priors$slope_mean, self$priors$slope_sd))
      
      # Define priors for brms
      priors <- c(
        brms::set_prior(sprintf("beta(%g, %g)", self$priors$floor_alpha, self$priors$floor_beta), 
                       nlpar = "floor", lb = 0, ub = 1),
        brms::set_prior(sprintf("beta(%g, %g)", self$priors$ceiling_alpha, self$priors$ceiling_beta),
                       nlpar = "ceiling", lb = 0, ub = 1),
        brms::set_prior(sprintf("normal(%g, %g)", self$priors$ec50_mean, self$priors$ec50_sd),
                       nlpar = "ec50"),
        brms::set_prior(sprintf("normal(%g, %g)", self$priors$slope_mean, self$priors$slope_sd),
                       nlpar = "slope", lb = 0)
      )
      
      # Define the non-linear formula
      # prob_infection = ceiling * (inv_logit(-slope * (titre - ec50)) * (1 - floor) + floor)
      formula <- brms::bf(
        infected ~ ceiling * (inv_logit(-slope * (titre - ec50)) * (1 - floor) + floor),
        floor ~ 1,
        ceiling ~ 1,
        ec50 ~ 1,
        slope ~ 1,
        nl = TRUE
      )
      
      # Fit the model
      # Handle cores argument - either from ... or default to chains
      dots <- list(...)
      if (!"cores" %in% names(dots)) {
        dots$cores <- min(chains, parallel::detectCores())
      }
      
      self$fit <- do.call(
        brms::brm,
        c(list(
          formula = formula,
          data = model_data,
          family = brms::bernoulli(link = "identity"),
          prior = priors,
          chains = chains,
          iter = iter,
          warmup = warmup,
          thin = thin,
          backend = "rstan"
        ), dots)
      )
      
      # Compute LOO-CV
      message("Computing LOO-CV...")
      self$loo <- brms::loo(self$fit)
      
      message("Model fitting complete!")
      invisible(self)
    },
    
    #' @description
    #' Get posterior predictions for infection probability
    #' @param newdata Optional new titre values for prediction
    #' @return Matrix of posterior predictions (iterations x observations)
    predict = function(newdata = NULL) {
      if (is.null(self$fit)) {
        stop("Model has not been fitted yet. Run fit_model() first.")
      }
      
      if (is.null(newdata)) {
        # Extract fitted probabilities using brms
        # Get posterior predictions (summary = FALSE gives us all samples)
        fitted_vals <- fitted(self$fit, summary = FALSE)
        return(fitted_vals)
      } else {
        # Predict for new data
        new_df <- data.frame(titre = newdata, infected = NA)
        predictions <- fitted(self$fit, newdata = new_df, summary = FALSE)
        return(predictions)
      }
    },
    
    #' @description
    #' Extract probability of protection from the fitted model
    #' @param newdata Optional vector of new titre values for prediction
    #' @return Matrix of protection probabilities (rows = MCMC samples, cols = observations)
    predict_protection = function(newdata = NULL) {
      if (is.null(self$fit)) {
        stop("Model has not been fitted yet. Run fit_model() first.")
      }
      
      # Get infection probabilities
      prob_infection <- self$predict(newdata = newdata)
      
      # Extract ceiling parameter samples from brms fit
      posterior_samples <- brms::as_draws_df(self$fit)
      ceiling_samples <- posterior_samples$b_ceiling_Intercept
      
      # Calculate protection: 1 - (prob_infection / ceiling)
      n_iter <- nrow(prob_infection)
      n_new <- ncol(prob_infection)
      prob_protection <- matrix(NA, nrow = n_iter, ncol = n_new)
      
      for (i in 1:n_iter) {
        prob_protection[i, ] <- 1 - (prob_infection[i, ] / ceiling_samples[i])
      }
      
      return(prob_protection)
    },
    
    #' @description
    #' Get summary statistics for model parameters
    #' @return Data frame with parameter summaries
    summary = function() {
      if (is.null(self$fit)) {
        stop("Model has not been fitted yet. Run fit_model() first.")
      }
      
      # Get brms summary for the non-linear parameters
      summary_df <- as.data.frame(brms::fixef(self$fit))
      
      return(summary_df)
    },
    
    #' @description
    #' Calculate performance metrics
    #' @return List containing ROC AUC, Brier score, and LOO-CV metrics
    get_metrics = function() {
      if (is.null(self$fit)) {
        stop("Model has not been fitted yet. Run fit_model() first.")
      }
      
      # Get posterior predictions
      prob_pred <- self$predict()
      prob_mean <- colMeans(prob_pred)
      
      # Calculate ROC AUC
      roc_obj <- pROC::roc(self$infected, prob_mean, quiet = TRUE)
      auc_value <- as.numeric(pROC::auc(roc_obj))
      
      # Calculate Brier score
      brier <- mean((prob_mean - self$infected)^2)
      
      # Get LOO metrics
      loo_elpd <- self$loo$estimates["elpd_loo", "Estimate"]
      loo_se <- self$loo$estimates["elpd_loo", "SE"]
      
      metrics <- list(
        roc_auc = auc_value,
        brier_score = brier,
        loo_elpd = loo_elpd,
        loo_se = loo_se,
        loo_object = self$loo
      )
      
      # Print summary
      cat("Performance Metrics:\n")
      cat(sprintf("  ROC AUC: %.3f\n", auc_value))
      cat(sprintf("  Brier Score: %.3f\n", brier))
      cat(sprintf("  LOO ELPD: %.2f (SE: %.2f)\n", loo_elpd, loo_se))
      
      invisible(metrics)
    },
    
    #' @description
    #' Plot the fitted curve with uncertainty
    #' @param title Plot title
    #' @param ... Additional arguments passed to ggplot2
    #' @return A ggplot object
    plot_curve = function(title = "Correlates of Risk Curve", ...) {
      if (is.null(self$fit)) {
        stop("Model has not been fitted yet. Run fit_model() first.")
      }
      
      # Create prediction grid
      titre_grid <- seq(min(self$titre), max(self$titre), length.out = 200)
      pred_matrix <- self$predict(newdata = titre_grid)
      
      # Calculate credible intervals
      pred_mean <- colMeans(pred_matrix)
      pred_lower <- apply(pred_matrix, 2, quantile, probs = 0.025)
      pred_upper <- apply(pred_matrix, 2, quantile, probs = 0.975)
      
      # Create plot data
      plot_df <- data.frame(
        titre = titre_grid,
        prob = pred_mean,
        lower = pred_lower,
        upper = pred_upper
      )
      
      obs_df <- data.frame(
        titre = self$titre,
        infected = self$infected
      )
      
      # Create plot
      p <- ggplot2::ggplot() +
        ggplot2::geom_ribbon(
          data = plot_df,
          ggplot2::aes(x = titre, ymin = lower, ymax = upper),
          alpha = 0.3, fill = "steelblue"
        ) +
        ggplot2::geom_line(
          data = plot_df,
          ggplot2::aes(x = titre, y = prob),
          color = "steelblue", linewidth = 1
        ) +
        ggplot2::geom_point(
          data = obs_df,
          ggplot2::aes(x = titre, y = infected),
          alpha = 0.3, position = ggplot2::position_jitter(height = 0.02)
        ) +
        ggplot2::labs(
          title = title,
          x = "Antibody Titre (log scale)",
          y = "Probability of Infection"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
          panel.grid.minor = ggplot2::element_blank()
        )
      
      return(p)
    },
    
    #' @description
    #' Plot ROC curve
    #' @param title Plot title
    #' @return A ggplot object
    plot_roc = function(title = "ROC Curve") {
      if (is.null(self$fit)) {
        stop("Model has not been fitted yet. Run fit_model() first.")
      }
      
      prob_pred <- colMeans(self$predict())
      roc_obj <- pROC::roc(self$infected, prob_pred, quiet = TRUE)
      auc_value <- as.numeric(pROC::auc(roc_obj))
      
      roc_df <- data.frame(
        fpr = 1 - roc_obj$specificities,
        tpr = roc_obj$sensitivities
      )
      
      p <- ggplot2::ggplot(roc_df, ggplot2::aes(x = fpr, y = tpr)) +
        ggplot2::geom_line(color = "steelblue", linewidth = 1) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.5) +
        ggplot2::labs(
          title = sprintf("%s (AUC = %.3f)", title, auc_value),
          x = "False Positive Rate",
          y = "True Positive Rate"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
          panel.grid.minor = ggplot2::element_blank()
        ) +
        ggplot2::coord_equal()
      
      return(p)
    },
    
    #' @description
    #' Plot posterior distributions of parameters
    #' @return A ggplot object
    plot_posteriors = function() {
      if (is.null(self$fit)) {
        stop("Model has not been fitted yet. Run fit_model() first.")
      }
      
      # Extract posterior samples from brms
      posterior_samples <- brms::as_draws_df(self$fit)
      
      param_df <- data.frame(
        floor = posterior_samples$b_floor_Intercept,
        ceiling = posterior_samples$b_ceiling_Intercept,
        ec50 = posterior_samples$b_ec50_Intercept,
        slope = posterior_samples$b_slope_Intercept
      )
      
      param_long <- tidyr::pivot_longer(
        param_df,
        cols = everything(),
        names_to = "parameter",
        values_to = "value"
      )
      
      p <- ggplot2::ggplot(param_long, ggplot2::aes(x = value)) +
        ggplot2::geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
        ggplot2::facet_wrap(~parameter, scales = "free") +
        ggplot2::labs(
          title = "Posterior Distributions of Parameters",
          x = "Parameter Value",
          y = "Frequency"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
        )
      
      return(p)
    },
    
    #' @description
    #' Extract the correlate of protection conditional on exposure.
    #' @param correlate_of_risk Numeric vector of correlates of risk.
    #' @param upper_bound Numeric value for the upper bound (default: 0.7).
    #' @return Numeric vector of correlates of protection.
    extract_cop = function(correlate_of_risk, upper_bound = 0.7) {
      if (missing(correlate_of_risk)) {
        stop("correlate_of_risk must be provided")
      }
      if (!is.numeric(correlate_of_risk) || any(correlate_of_risk < 0 | correlate_of_risk > 1)) {
        stop("correlate_of_risk must be a numeric vector with values between 0 and 1")
      }
      if (!is.numeric(upper_bound) || upper_bound <= 0) {
        stop("upper_bound must be a positive numeric value")
      }
      cop <- (1 - correlate_of_risk) / upper_bound
      return(cop)
    }
  ),
  
  private = list(
    # Get default prior distributions based on data
    # Internal method - generates default prior distributions based on the data
    # Returns a list of default prior parameters
    get_default_priors = function() {
      # Calculate data-driven defaults for ec50
      titre_midpoint <- (max(self$titre) + min(self$titre)) / 2
      titre_range <- max(self$titre) - min(self$titre)
      titre_sd <- titre_range / 4  # Cover most of the range
      
      list(
        floor_alpha = 1,            # Weak prior favoring low floor
        floor_beta = 9,
        ceiling_alpha = 9,          # Weak prior favoring high ceiling
        ceiling_beta = 1,
        ec50_mean = titre_midpoint, # Centered on data midpoint
        ec50_sd = titre_sd,         # SD based on data range
        slope_mean = 0,             # Centered at 0 (no directional bias)
        slope_sd = 2                # Weakly informative
      )
    }
  )
)
