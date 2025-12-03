#' SeroCOPMulti R6 Class for Multi-Biomarker Correlates of Protection Analysis
#'
#' @description
#' An R6 class for analyzing multiple biomarkers simultaneously.
#' Fits separate models for each biomarker and provides comparison tools.
#'
#' @concept r6-classes
#' @export
#' @examples
#' \dontrun{
#' # Create multi-biomarker data
#' titres <- matrix(rnorm(600), ncol = 3)
#' colnames(titres) <- c("IgG", "IgA", "Neutralization")
#' infected <- rbinom(200, 1, 0.3)
#'
#' # Initialize and fit
#' multi_model <- SeroCOPMulti$new(titre = titres, infected = infected)
#' multi_model$fit_all(chains = 4, iter = 2000)
#' 
#' # Compare biomarkers
#' multi_model$compare_biomarkers()
#' multi_model$plot_comparison()
#' }
SeroCOPMulti <- R6::R6Class(
  "SeroCOPMulti",
  
  public = list(
    #' @field titre Matrix of antibody titres (columns = biomarkers)
    titre = NULL,
    
    #' @field infected Binary vector of infection status (0/1)
    infected = NULL,
    
    #' @field biomarker_names Names of biomarkers
    biomarker_names = NULL,
    
    #' @field models List of fitted SeroCOP objects
    models = NULL,
    
    #' @description
    #' Create a new SeroCOPMulti object
    #' @param titre Matrix of antibody titres (rows = samples, cols = biomarkers)
    #' @param infected Binary vector (0/1) of infection outcomes
    #' @param biomarker_names Optional vector of biomarker names
    #' @return A new SeroCOPMulti object
    initialize = function(titre, infected, biomarker_names = NULL) {
      # Convert to matrix if needed
      if (is.vector(titre)) {
        titre <- matrix(titre, ncol = 1)
      }
      
      if (!is.matrix(titre) && !is.data.frame(titre)) {
        stop("titre must be a matrix, data frame, or vector")
      }
      
      titre <- as.matrix(titre)
      
      # Validate inputs
      if (!is.numeric(titre)) {
        stop("titre must be numeric")
      }
      if (!all(infected %in% c(0, 1))) {
        stop("infected must be binary (0/1)")
      }
      if (nrow(titre) != length(infected)) {
        stop("Number of rows in titre must equal length of infected")
      }
      if (any(is.na(titre)) || any(is.na(infected))) {
        stop("Missing values are not allowed")
      }
      
      self$titre <- titre
      self$infected <- infected
      
      # Set biomarker names
      if (is.null(biomarker_names)) {
        if (!is.null(colnames(titre))) {
          self$biomarker_names <- colnames(titre)
        } else {
          self$biomarker_names <- paste0("Biomarker", 1:ncol(titre))
        }
      } else {
        if (length(biomarker_names) != ncol(titre)) {
          stop("Length of biomarker_names must equal number of columns in titre")
        }
        self$biomarker_names <- biomarker_names
      }
      
      message(sprintf("SeroCOPMulti initialized with %d observations and %d biomarkers", 
                      nrow(titre), ncol(titre)))
      message(sprintf("  Biomarkers: %s", paste(self$biomarker_names, collapse = ", ")))
      message(sprintf("  Infection rate: %.1f%%", mean(infected) * 100))
    },
    
    #' @description
    #' Fit models for all biomarkers
    #' @param chains Number of MCMC chains (default: 4)
    #' @param iter Number of iterations per chain (default: 2000)
    #' @param warmup Number of warmup iterations (default: iter/2)
    #' @param cores Number of cores for parallel processing (default: 1)
    #' @param ... Additional arguments passed to rstan::sampling
    #' @return Self (invisibly)
    fit_all = function(chains = 4, iter = 2000, warmup = floor(iter/2), 
                      cores = 1, ...) {
      
      n_biomarkers <- ncol(self$titre)
      self$models <- vector("list", n_biomarkers)
      names(self$models) <- self$biomarker_names
      
      message(sprintf("\nFitting models for %d biomarkers...\n", n_biomarkers))
      
      for (i in 1:n_biomarkers) {
        biomarker <- self$biomarker_names[i]
        message(sprintf("=== Fitting %s (%d/%d) ===", biomarker, i, n_biomarkers))
        
        # Create individual SeroCOP model
        model <- SeroCOP$new(
          titre = self$titre[, i],
          infected = self$infected
        )
        
        # Fit the model
        model$fit_model(
          chains = chains,
          iter = iter,
          warmup = warmup,
          ...
        )
        
        self$models[[i]] <- model
        message("")
      }
      
      message("All models fitted successfully!")
      invisible(self)
    },
    
    #' @description
    #' Compare biomarkers on AUC and LOO metrics
    #' @return Data frame with comparison metrics
    compare_biomarkers = function() {
      if (is.null(self$models)) {
        stop("Models have not been fitted yet. Run fit_all() first.")
      }
      
      n_biomarkers <- length(self$models)
      
      # Extract metrics for each biomarker
      results <- data.frame(
        biomarker = character(n_biomarkers),
        auc = numeric(n_biomarkers),
        auc_lower = numeric(n_biomarkers),
        auc_upper = numeric(n_biomarkers),
        brier = numeric(n_biomarkers),
        loo_elpd = numeric(n_biomarkers),
        loo_se = numeric(n_biomarkers),
        stringsAsFactors = FALSE
      )
      
      for (i in 1:n_biomarkers) {
        model <- self$models[[i]]
        
        # Get predictions
        prob_pred <- colMeans(model$predict())
        
        # Calculate AUC with CI
        roc_obj <- pROC::roc(model$infected, prob_pred, quiet = TRUE)
        auc_ci <- pROC::ci.auc(roc_obj, conf.level = 0.95)
        
        results$biomarker[i] <- self$biomarker_names[i]
        results$auc[i] <- as.numeric(pROC::auc(roc_obj))
        results$auc_lower[i] <- auc_ci[1]
        results$auc_upper[i] <- auc_ci[3]
        results$brier[i] <- mean((prob_pred - model$infected)^2)
        results$loo_elpd[i] <- model$loo$estimates["elpd_loo", "Estimate"]
        results$loo_se[i] <- model$loo$estimates["elpd_loo", "SE"]
      }
      
      # Print comparison
      cat("\n=== Biomarker Comparison ===\n\n")
      print(results, row.names = FALSE)
      cat("\n")
      
      invisible(results)
    },
    
    #' @description
    #' Plot biomarkers on AUC vs LOO plane
    #' @param add_labels Logical, add biomarker labels (default: TRUE)
    #' @return A ggplot object
    plot_comparison = function(add_labels = TRUE) {
      if (is.null(self$models)) {
        stop("Models have not been fitted yet. Run fit_all() first.")
      }
      
      comparison <- self$compare_biomarkers()
      
      p <- ggplot2::ggplot(comparison, 
                           ggplot2::aes(x = loo_elpd, y = auc, label = biomarker)) +
        ggplot2::geom_point(size = 4, color = "steelblue") +
        ggplot2::geom_errorbar(
          ggplot2::aes(ymin = auc_lower, ymax = auc_upper),
          width = 0, color = "steelblue", alpha = 0.5
        ) +
        ggplot2::geom_errorbarh(
          ggplot2::aes(xmin = loo_elpd - loo_se, xmax = loo_elpd + loo_se),
          height = 0, color = "steelblue", alpha = 0.5
        ) +
        ggplot2::labs(
          title = "Biomarker Comparison: AUC vs LOO-ELPD",
          x = "LOO-ELPD (higher = better)",
          y = "ROC AUC (higher = better)"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
        )
      
      if (add_labels) {
        p <- p + ggplot2::geom_text(
          vjust = -0.8, 
          hjust = 0.5,
          size = 3.5
        )
      }
      
      return(p)
    },
    
    #' @description
    #' Plot all fitted curves
    #' @return A ggplot object
    plot_all_curves = function() {
      if (is.null(self$models)) {
        stop("Models have not been fitted yet. Run fit_all() first.")
      }
      
      n_biomarkers <- length(self$models)
      
      # Combine data for all biomarkers
      plot_data <- list()
      obs_data <- list()
      
      for (i in 1:n_biomarkers) {
        model <- self$models[[i]]
        
        # Create prediction grid
        titre_grid <- seq(min(model$titre), max(model$titre), length.out = 100)
        pred_matrix <- model$predict(newdata = titre_grid)
        
        pred_mean <- colMeans(pred_matrix)
        pred_lower <- apply(pred_matrix, 2, quantile, probs = 0.025)
        pred_upper <- apply(pred_matrix, 2, quantile, probs = 0.975)
        
        plot_data[[i]] <- data.frame(
          biomarker = self$biomarker_names[i],
          titre = titre_grid,
          prob = pred_mean,
          lower = pred_lower,
          upper = pred_upper
        )
        
        obs_data[[i]] <- data.frame(
          biomarker = self$biomarker_names[i],
          titre = model$titre,
          infected = model$infected
        )
      }
      
      plot_df <- do.call(rbind, plot_data)
      obs_df <- do.call(rbind, obs_data)
      
      p <- ggplot2::ggplot() +
        ggplot2::geom_ribbon(
          data = plot_df,
          ggplot2::aes(x = titre, ymin = lower, ymax = upper, fill = biomarker),
          alpha = 0.2
        ) +
        ggplot2::geom_line(
          data = plot_df,
          ggplot2::aes(x = titre, y = prob, color = biomarker),
          linewidth = 1
        ) +
        ggplot2::geom_point(
          data = obs_df,
          ggplot2::aes(x = titre, y = infected),
          alpha = 0.2, size = 0.8,
          position = ggplot2::position_jitter(height = 0.02)
        ) +
        ggplot2::facet_wrap(~biomarker, scales = "free_x") +
        ggplot2::labs(
          title = "Correlates of Risk Curves by Biomarker",
          x = "Antibody Titre (log scale)",
          y = "Probability of Infection"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
          legend.position = "none"
        )
      
      return(p)
    },
    
    #' @description
    #' Extract the correlate of protection conditional on exposure for all biomarkers.
    #' @param correlate_of_risk_list List of numeric vectors of correlates of risk for each biomarker.
    #' @param upper_bound Numeric value for the upper bound (default: 0.7).
    #' @return List of numeric vectors of correlates of protection for each biomarker.
    extract_cop_multi = function(correlate_of_risk_list, upper_bound = 0.7) {
      if (!is.list(correlate_of_risk_list)) {
        stop("correlate_of_risk_list must be a list of numeric vectors")
      }
      if (!is.numeric(upper_bound) || upper_bound <= 0) {
        stop("upper_bound must be a positive numeric value")
      }
      cop_list <- lapply(correlate_of_risk_list, function(cor) {
        if (!is.numeric(cor) || any(cor < 0 | cor > 1)) {
          stop("Each correlate_of_risk must be a numeric vector with values between 0 and 1")
        }
        1 - (cor / upper_bound)
      })
      names(cop_list) <- names(correlate_of_risk_list)
      return(cop_list)
    }
  )
)
