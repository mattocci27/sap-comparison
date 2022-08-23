data {
  int<lower=0> N; // number of samples
  int<lower=0> J; // number of species
  int<lower=0> K; // number of pressure types
  int<lower=0> JK; // number of pressure types
  array[N] int<lower=1,upper=J> jj; // species
  array[N] int<lower=1,upper=K> kk; // pressure
  array[N] int<lower=1,upper=JK> jk; // pressure
  vector[N] y; // pres_calib - tens_calib
}

parameters{
  real mu_hat;
  real<lower=0> sigma;
  vector[J] alpha_raw;
  vector[K] beta_raw;
  vector<lower=0,upper=pi()/2>[2] tau_unif;
}

transformed parameters{
  vector[J] alpha;
  vector[K] beta;
  vector<lower=0>[2] tau;
  for (i in 1:2) tau[i] = 2.5 * tan(tau_unif[i]);
  alpha = alpha_raw * tau[1];
  beta = beta_raw * tau[2];
}

model {
  vector[N] mu;
  to_vector(alpha_raw) ~ std_normal();
  to_vector(beta_raw) ~ std_normal();
  for (n in 1:N) {
    mu[n] = mu_hat + alpha[jj[n]] + beta[kk[n]];
  }
  y ~ normal(mu, sigma);
}

generated quantities {
  vector[N] log_lik;
  vector[JK] pred;
  vector[JK] effect;
  real beta_12 = beta[1] - beta[2];
  real beta_13 = beta[1] - beta[3];
  real beta_23 = beta[2] - beta[3];
  real alpha_12 = alpha[1] - alpha[2];
  real alpha_13 = alpha[1] - alpha[3];
  real alpha_14 = alpha[1] - alpha[4];
  real alpha_15 = alpha[1] - alpha[5];
  real alpha_23 = alpha[2] - alpha[3];
  real alpha_24 = alpha[2] - alpha[4];
  real alpha_25 = alpha[2] - alpha[5];
  real alpha_34 = alpha[3] - alpha[4];
  real alpha_35 = alpha[3] - alpha[5];
  real alpha_45 = alpha[4] - alpha[5];
  for (n in 1:N) {
    pred[jk[n]] = mu_hat + alpha[jj[n]] + beta[kk[n]];
    effect[jk[n]] = alpha[jj[n]] + beta[kk[n]];
    log_lik[n] = normal_lpdf(y[n] | pred[jk[n]], sigma);
  }
}

