data {
  int<lower=1> N;
  int<lower=2, upper=4> D;
  matrix[N, D] titre;
  array[N] int<lower=0, upper=1> infected;

  real<lower=0> floor_alpha;
  real<lower=0> floor_beta;
  real<lower=0> ceiling_alpha;
  real<lower=0> ceiling_beta;
  real ec50_mean;
  real<lower=0> ec50_sd;
  real<lower=0> slope_scale;
}

parameters {
  real<lower=0, upper=1> floor;
  real<lower=0, upper=1> ceiling;
  vector[D] ec50;
  vector<lower=0>[D] slope;
}

transformed parameters {
  vector[N] eta;
  vector[N] prob_infection;
  vector[N] prob_protection;

  for (n in 1:N) {
    eta[n] = 0;
    for (d in 1:D) {
      eta[n] += slope[d] * (titre[n, d] - ec50[d]);
    }
    prob_infection[n] = ceiling * (inv_logit(-eta[n]) * (1 - floor) + floor);
    prob_protection[n] = 1 - prob_infection[n] / ceiling;
  }
}

model {
  floor ~ beta(floor_alpha, floor_beta);
  ceiling ~ beta(ceiling_alpha, ceiling_beta);
  ec50 ~ normal(ec50_mean, ec50_sd);
  slope ~ lognormal(log(slope_scale), 0.7);

  infected ~ bernoulli(prob_infection);
}

generated quantities {
  vector[N] log_lik;
  array[N] int infected_rep;

  for (n in 1:N) {
    log_lik[n] = bernoulli_lpmf(infected[n] | prob_infection[n]);
    infected_rep[n] = bernoulli_rng(prob_infection[n]);
  }
}