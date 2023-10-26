data {
  int<lower=1> N;  // number of observations
  int<lower=1> K; // number of predictors
  int<lower=1> M;  // number of unique trees
  vector[N] log_k;  // response variable log(k)
  vector[N] log_k_ref;  // offset
  matrix[N, K] x;  // predictor
  array[N] int<lower=1, upper=M> tree; // integer
}

parameters {
  vector[K] beta;
  real<lower=0> sigma;  // error SD
  vector[M] phi_raw;
  real<lower=0> tau;
}

transformed parameters {
  vector[M] phi = phi_raw * tau;  // tree random effects
}

model {
  beta ~ normal(0, 2.5);
  tau ~ normal(0, 2.5);
  phi_raw ~ normal(0, 1);
  log_k ~ normal(log_k_ref + x * beta + phi[tree], sigma);
}

// generated quantities {
//   vector[N] log_lik;
//   for (n in 1:N) log_lik[n] = normal_lpdf(log_k[n] |
//     log_k_ref + x[n] * beta + phi[tree[n]], sigma);
// }
