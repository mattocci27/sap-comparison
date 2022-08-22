data {
  int<lower=0> N; // number of samples
  int<lower=0> J; // number of species
  int<lower=0> K; // number of pressure types
  int<lower=0> JK; // number of species x pressure types
  array[N] int<lower=1,upper=J> jj; // species
  array[N] int<lower=1,upper=K> kk; // pressure
  array[N] int<lower=1,upper=JK> jk; // species x pressure
  vector[N] y; // pres_calib - tens_calib
}

parameters{
  real mu_hat;
  real<lower=0> sigma;
  vector[J] alpha;
  vector[K] beta;
  // vector[JK] gamma;
}

// transformed parameters{
//   vector[J] alpha;
//   vector[K] beta;
//   vector[JK] gamma;
//   vector<lower=0>[3] tau;
//   for (i in 1:3) tau[i] = 2.5 * tan(tau_unif[i]);
//   alpha = alpha_raw * tau[1];
//   beta = beta_raw * tau[2];
//   gamma = gamma_raw * tau[3];
// }

model {
  vector[N] mu;
  to_vector(alpha) ~ normal(0, 5);
  to_vector(beta) ~ normal(0, 10);
  // to_vector(gamma) ~ normal(0, 5);
  // mu = mu_hat + alpha[jj] + beta[kk] + gamma[jk];
  mu = mu_hat + alpha[jj] + beta[kk];
  y ~ normal(mu, sigma);
}

// generated quantities {
//   vector[N] log_lik;
//   vector[JK] pred;
//   vector[JK] effect;
//   for (n in 1:N) {
//     pred[jk[n]] = mu_hat + alpha[jj[n]] + beta[kk[n]] + gamma[jk[n]];
//     effect[jk[n]] = alpha[jj[n]] + beta[kk[n]] + gamma[jk[n]];
//     log_lik[n] = normal_lpdf(y[n] | pred[jk[n]], sigma);
//   }
// }

