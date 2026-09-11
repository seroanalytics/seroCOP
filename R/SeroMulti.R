#' SeroMulti R6 Class for Multidimensional Correlates of Protection
#'
#' @description
#' Fits a joint Bayesian four-parameter logistic model to two to four biomarkers
#' measured for each person at the same time point. The model estimates one
#' protective slope and EC50 per biomarker and can display a two-biomarker CoP
#' surface.
#'
#' @concept r6-classes
#' @export
SeroMulti <- R6::R6Class(
  "SeroMulti",

  public = list(
    #' @field titre Numeric matrix of log antibody titres.
    titre = NULL,
    #' @field infected Binary infection outcomes.
    infected = NULL,
    #' @field biomarker_names Names of biomarker columns.
    biomarker_names = NULL,
    #' @field fit Fitted Stan model object after calling fit_model().
    fit = NULL,
    #' @field loo Leave-one-out cross-validation result after fitting.
    loo = NULL,
    #' @field priors Prior distribution parameters.
    priors = NULL,

    #' @description Create a multidimensional CoP model.
    #' @param titre Numeric matrix with one row per person and 2--4 biomarker columns.
    #' @param infected Binary infection outcome vector.
    #' @param biomarker_names Optional names for the biomarker columns.
    #' @return A new SeroMulti object.
    initialize = function(titre, infected, biomarker_names = NULL) {
      if (!is.matrix(titre) && !is.data.frame(titre)) {
        stop("titre must be a numeric matrix or data frame")
      }
      titre <- as.matrix(titre)
      if (!is.numeric(titre)) stop("titre must be numeric")
      if (ncol(titre) < 2 || ncol(titre) > 4) {
        stop("titre must contain between 2 and 4 biomarker columns")
      }
      if (length(infected) != nrow(titre)) {
        stop("Number of rows in titre must equal length of infected")
      }
      if (!all(infected %in% c(0, 1))) stop("infected must be binary (0/1)")
      if (anyNA(titre) || anyNA(infected)) stop("Missing values are not allowed")

      if (is.null(biomarker_names)) {
        biomarker_names <- colnames(titre)
        if (is.null(biomarker_names)) biomarker_names <- paste0("Biomarker", seq_len(ncol(titre)))
      }
      if (length(biomarker_names) != ncol(titre) || anyDuplicated(biomarker_names)) {
        stop("biomarker_names must be unique and match the number of titre columns")
      }
      colnames(titre) <- biomarker_names

      self$titre <- titre
      self$infected <- as.integer(infected)
      self$biomarker_names <- biomarker_names
      self$priors <- private$default_priors()
      message(sprintf("SeroMulti initialized with %d observations and %d biomarkers",
                      nrow(titre), ncol(titre)))
    },

    #' @description Update prior distributions before fitting.
    #' @param floor_alpha,floor_beta Beta prior parameters for floor.
    #' @param ceiling_alpha,ceiling_beta Beta prior parameters for ceiling.
    #' @param ec50_mean,ec50_sd Normal prior parameters shared across EC50 values.
    #' @param slope_scale Median of the log-normal slope prior.
    #' @return Self invisibly.
    definePrior = function(floor_alpha = NULL, floor_beta = NULL,
                           ceiling_alpha = NULL, ceiling_beta = NULL,
                           ec50_mean = NULL, ec50_sd = NULL,
                           slope_scale = NULL) {
      updated <- utils::modifyList(self$priors, Filter(Negate(is.null), list(
        floor_alpha = floor_alpha, floor_beta = floor_beta,
        ceiling_alpha = ceiling_alpha, ceiling_beta = ceiling_beta,
        ec50_mean = ec50_mean, ec50_sd = ec50_sd, slope_scale = slope_scale
      )))
      if (any(unlist(updated[c("floor_alpha", "floor_beta", "ceiling_alpha", "ceiling_beta", "ec50_sd", "slope_scale")]) <= 0)) {
        stop("Prior scale and shape parameters must be positive")
      }
      self$priors <- updated
      invisible(self)
    },

    #' @description Fit the joint multidimensional Stan model.
    #' @param chains Number of MCMC chains.
    #' @param iter Number of iterations per chain.
    #' @param warmup Number of warmup iterations per chain.
    #' @param cores Number of parallel chains.
    #' @param ... Additional arguments passed to rstan::sampling.
    #' @return Self invisibly.
    fit_model = function(chains = 4, iter = 2000, warmup = floor(iter / 2),
                         cores = 1, ...) {
      stan_data <- c(list(N = nrow(self$titre), D = ncol(self$titre),
                          titre = self$titre, infected = self$infected), self$priors)
      self$fit <- rstan::sampling(
        object = private$stan_model(), data = stan_data, chains = chains,
        iter = iter, warmup = warmup, cores = cores, ...
      )
      log_lik <- rstan::extract(self$fit, pars = "log_lik", permuted = FALSE)
      self$loo <- loo::loo(log_lik)
      invisible(self)
    },

    #' @description Predict infection probabilities from posterior draws.
    #' @param newdata Optional numeric matrix with the same biomarker columns.
    #' @return Matrix of posterior draws by observations.
    predict = function(newdata = NULL) {
      if (is.null(self$fit)) stop("Model has not been fitted yet. Run fit_model() first.")
      if (is.null(newdata)) newdata <- self$titre
      newdata <- private$validate_newdata(newdata)
      draws <- rstan::extract(self$fit, pars = c("floor", "ceiling", "ec50", "slope"))
      predictions <- matrix(NA_real_, nrow = length(draws$floor), ncol = nrow(newdata))
      for (draw in seq_len(nrow(predictions))) {
        eta <- (sweep(newdata, 2, draws$ec50[draw, ], "-") *
          rep(draws$slope[draw, ], each = nrow(newdata)))
        eta <- rowSums(eta)
        predictions[draw, ] <- draws$ceiling[draw] *
          (stats::plogis(-eta) * (1 - draws$floor[draw]) + draws$floor[draw])
      }
      predictions
    },

    #' @description Predict correlate-of-protection values from posterior draws.
    #' @param newdata Optional numeric matrix with the same biomarker columns.
    #' @return Matrix of posterior draws by observations.
    predict_protection = function(newdata = NULL) {
      probabilities <- self$predict(newdata)
      ceiling <- rstan::extract(self$fit, pars = "ceiling")$ceiling
      1 - probabilities / ceiling
    },

    #' @description Summarise the two-biomarker CoP surface.
    #' @param grid_size Number of values along each biomarker axis.
    #' @return Data frame containing titres and posterior mean CoP.
    surface_data = function(grid_size = 50) {
      if (ncol(self$titre) != 2) stop("surface_data is available only for two biomarkers")
      if (!is.numeric(grid_size) || length(grid_size) != 1 || grid_size < 2) {
        stop("grid_size must be at least 2")
      }
      grid <- expand.grid(
        seq(min(self$titre[, 1]), max(self$titre[, 1]), length.out = grid_size),
        seq(min(self$titre[, 2]), max(self$titre[, 2]), length.out = grid_size)
      )
      names(grid) <- self$biomarker_names
      protection <- self$predict_protection(as.matrix(grid))
      grid$cop <- colMeans(protection)
      grid
    },

    #' @description Plot the CoP surface for a two-biomarker model.
    #' @param type Either a filled contour plot or a base-R perspective plot.
    #' @param grid_size Number of values along each biomarker axis.
    #' @return A ggplot object for contour plots; invisibly returns surface data for perspective plots.
    plot_surface = function(type = c("contour", "persp"), grid_size = 50) {
      type <- match.arg(type)
      surface <- self$surface_data(grid_size)
      if (type == "contour") {
        return(ggplot2::ggplot(surface, ggplot2::aes_string(
          x = self$biomarker_names[1], y = self$biomarker_names[2], fill = "cop"
        )) +
          ggplot2::geom_raster(interpolate = TRUE) +
          ggplot2::geom_contour(ggplot2::aes_string(z = "cop"), color = "white") +
          ggplot2::scale_fill_viridis_c(name = "CoP") +
          ggplot2::labs(title = "Multidimensional Correlate of Protection",
                         x = self$biomarker_names[1], y = self$biomarker_names[2]) +
          ggplot2::theme_minimal())
      }
      z <- matrix(surface$cop, nrow = grid_size, ncol = grid_size)
      graphics::persp(unique(surface[[1]]), unique(surface[[2]]), z,
                      xlab = self$biomarker_names[1], ylab = self$biomarker_names[2],
                      zlab = "CoP", theta = 35, phi = 25, col = "steelblue")
      invisible(surface)
    }
  ),

  private = list(
    default_priors = function() {
      titre_range <- range(self$titre)
      list(floor_alpha = 1, floor_beta = 9, ceiling_alpha = 9,
           ceiling_beta = 1, ec50_mean = mean(titre_range),
           ec50_sd = diff(titre_range) / 4, slope_scale = 1)
    },
    stan_model = function() {
      stan_file <- system.file("stan", "multi_logistic_model.stan", package = "seroCOP")
      if (!nzchar(stan_file)) stan_file <- file.path("inst", "stan", "multi_logistic_model.stan")
      rstan::stan_model(stan_file)
    },
    validate_newdata = function(newdata) {
      if (!is.matrix(newdata) && !is.data.frame(newdata)) {
        stop("newdata must be a numeric matrix or data frame")
      }
      newdata <- as.matrix(newdata)
      if (!is.numeric(newdata) || ncol(newdata) != ncol(self$titre) || anyNA(newdata)) {
        stop("newdata must be a complete numeric matrix with the same number of biomarker columns")
      }
      newdata
    }
  )
)