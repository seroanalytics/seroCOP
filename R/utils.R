#' Simulate data for correlates of protection analysis
#'
#' @description
#' Generate synthetic data from a four-parameter logistic model for testing
#' and validation purposes.
#'
#' @param n Number of observations
#' @param floor Lower asymptote (minimum infection probability)
#' @param ceiling Upper asymptote (maximum infection probability)
#' @param ec50 Titre at 50% between floor and ceiling
#' @param slope Steepness of the curve
#' @param titre_mean Mean of titre distribution (normal)
#' @param titre_sd Standard deviation of titre distribution
#' @param seed Random seed for reproducibility
#'
#' @return A list containing:
#'   \item{titre}{Vector of simulated titres}
#'   \item{infected}{Vector of simulated infection outcomes}
#'   \item{prob_true}{True infection probabilities}
#'   \item{params}{List of true parameters used for simulation}
#'
#' @export
#' @examples
#' # Simulate data with default parameters
#' data <- simulate_serocop_data(n = 200, seed = 123)
#' 
#' # Custom parameters
#' data <- simulate_serocop_data(
#'   n = 500,
#'   floor = 0.02,
#'   ceiling = 0.95,
#'   ec50 = 2.0,
#'   slope = 3.0,
#'   seed = 456
#' )
simulate_serocop_data <- function(n = 200,
                                  floor = 0.05,
                                  ceiling = 0.9,
                                  ec50 = 1.5,
                                  slope = 2.0,
                                  titre_mean = 2.0,
                                  titre_sd = 1.5,
                                  seed = NULL) {
  
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Validate parameters
  if (floor < 0 || floor > 1) stop("floor must be between 0 and 1")
  if (ceiling < 0 || ceiling > 1) stop("ceiling must be between 0 and 1")
  if (floor >= ceiling) stop("floor must be less than ceiling")
  if (slope <= 0) stop("slope must be positive")
  if (n <= 0) stop("n must be positive")
  
  # Generate titres
  titre <- stats::rnorm(n, mean = titre_mean, sd = titre_sd)
  
  # Calculate infection probabilities
  prob_true <- floor + (ceiling - floor) / (1 + exp(slope * (titre - ec50)))
  
  # Generate infection outcomes
  infected <- stats::rbinom(n, size = 1, prob = prob_true)
  
  # Return results
  list(
    titre = titre,
    infected = infected,
    prob_true = prob_true,
    params = list(
      floor = floor,
      ceiling = ceiling,
      ec50 = ec50,
      slope = slope,
      titre_mean = titre_mean,
      titre_sd = titre_sd
    )
  )
}


#' Calculate Brier score
#'
#' @description
#' Calculate the Brier score for probabilistic predictions
#'
#' @param observed Binary vector of observed outcomes (0/1)
#' @param predicted Vector of predicted probabilities
#'
#' @return Numeric Brier score (lower is better)
#' @export
#' @examples
#' observed <- c(0, 1, 1, 0, 1)
#' predicted <- c(0.1, 0.8, 0.7, 0.2, 0.9)
#' brier_score(observed, predicted)
brier_score <- function(observed, predicted) {
  if (length(observed) != length(predicted)) {
    stop("observed and predicted must have the same length")
  }
  if (!all(observed %in% c(0, 1))) {
    stop("observed must be binary (0/1)")
  }
  if (any(predicted < 0 | predicted > 1)) {
    stop("predicted must be probabilities between 0 and 1")
  }
  
  mean((predicted - observed)^2)
}


#' Extract parameter estimates from SeroCOP fit
#'
#' @description
#' Extract posterior summaries of model parameters
#'
#' @param serocop_obj A fitted SeroCOP object
#' @param prob Probability for credible intervals (default 0.95)
#'
#' @return Data frame with parameter estimates and credible intervals
#' @export
extract_parameters <- function(serocop_obj, prob = 0.95) {
  if (!inherits(serocop_obj, "SeroCOP")) {
    stop("Input must be a SeroCOP object")
  }
  if (is.null(serocop_obj$fit)) {
    stop("Model has not been fitted yet")
  }
  
  alpha <- (1 - prob) / 2
  probs <- c(alpha, 0.5, 1 - alpha)
  
  # Extract posterior samples from brms
  posterior_samples <- brms::as_draws_df(serocop_obj$fit)
  
  params <- list(
    floor = posterior_samples$b_floor_Intercept,
    ceiling = posterior_samples$b_ceiling_Intercept,
    ec50 = posterior_samples$b_ec50_Intercept,
    slope = posterior_samples$b_slope_Intercept
  )
  
  result <- data.frame(
    parameter = c("floor", "ceiling", "ec50", "slope"),
    mean = c(mean(params$floor), mean(params$ceiling), 
             mean(params$ec50), mean(params$slope)),
    median = c(stats::median(params$floor), stats::median(params$ceiling),
               stats::median(params$ec50), stats::median(params$slope)),
    lower = c(stats::quantile(params$floor, probs[1]), 
              stats::quantile(params$ceiling, probs[1]),
              stats::quantile(params$ec50, probs[1]), 
              stats::quantile(params$slope, probs[1])),
    upper = c(stats::quantile(params$floor, probs[3]), 
              stats::quantile(params$ceiling, probs[3]),
              stats::quantile(params$ec50, probs[3]), 
              stats::quantile(params$slope, probs[3]))
  )
  
  rownames(result) <- NULL
  return(result)
}
