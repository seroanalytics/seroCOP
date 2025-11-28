// Four-parameter logistic model for correlates of protection
// Models the probability of infection as a function of antibody titre
// Parameters:
//   - floor: proportion of maximum risk remaining at high titre (relative protection at high titre)
//   - ceiling: maximum infection probability (at low antibody titre)
//   - ec50: titre at inflection point (50% reduction from ceiling to ceiling*floor)
//   - slope: steepness of the protective curve (higher = steeper decline in risk with titre)

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
  real<lower=0,upper=1> floor;       // proportion of maximum risk at high titre
  real<lower=0,upper=1> ceiling;     // maximum infection probability at low titre
  real ec50;                         // titre at inflection point
  real<lower=0> slope;               // steepness of protective curve
}

transformed parameters {
  vector[N] prob_infection;
  vector[N] prob_protection;
  // Four-parameter logistic function
  // At low titre: prob → ceiling (maximum risk)
  // At high titre: prob → ceiling * floor (minimum risk, as proportion of ceiling)
  for (n in 1:N) {
    prob_infection[n] = ceiling * (inv_logit(-slope * (titre[n] - ec50)) * (1 - floor) + floor);
    prob_protection[n] = 1 - (prob_infection[n] / ceiling);
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
