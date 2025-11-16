// Four-parameter logistic model for correlates of protection
// Models the probability of infection as a function of antibody titre
// Parameters:
//   - floor: lower asymptote (minimum infection probability)
//   - ceiling: upper asymptote (maximum infection probability)
//   - ec50: titre at 50% between floor and ceiling (inflection point)
//   - slope: steepness of the curve

data {
  int<lower=0> N;                    // number of observations
  vector[N] titre;                   // log antibody titres
  array[N] int<lower=0,upper=1> infected;  // infection status (0/1)
  
  // Prior parameters for floor (beta distribution)
  real<lower=0> floor_alpha;
  real<lower=0> floor_beta;
  
  // Prior parameters for ceiling (beta distribution)
  real<lower=0> ceiling_alpha;
  real<lower=0> ceiling_beta;
  
  // Prior parameters for ec50 (normal distribution)
  real ec50_mean;
  real<lower=0> ec50_sd;
  
  // Prior parameters for slope (normal distribution, truncated at 0)
  real slope_mean;
  real<lower=0> slope_sd;
}

parameters {
  real<lower=0,upper=1> floor;       // lower asymptote
  real<lower=0,upper=1> ceiling;     // upper asymptote
  real ec50;                         // titre at midpoint
  real<lower=0> slope;               // slope parameter (steepness)
}

transformed parameters {
  vector[N] prob_infection;
  
  // Four-parameter logistic function
  for (n in 1:N) {
    prob_infection[n] = floor + (ceiling - floor) / (1 + exp(slope * (titre[n] - ec50)));
  }
}

model {
  // Priors (using parameters passed from R)
  floor ~ beta(floor_alpha, floor_beta);
  ceiling ~ beta(ceiling_alpha, ceiling_beta);
  ec50 ~ normal(ec50_mean, ec50_sd);
  slope ~ normal(slope_mean, slope_sd);
  
  // Likelihood
  infected ~ bernoulli(prob_infection);
}

generated quantities {
  vector[N] log_lik;                 // log-likelihood for LOO-CV
  array[N] int infected_rep;         // posterior predictive samples
  
  for (n in 1:N) {
    log_lik[n] = bernoulli_lpmf(infected[n] | prob_infection[n]);
    infected_rep[n] = bernoulli_rng(prob_infection[n]);
  }
}
