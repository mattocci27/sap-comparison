data {
  int<lower=0> N; // number of samples
  int<lower=0> J; // number of species
  int<lower=0> K; // number of pressure types
  int<lower=0> JK; // number of species x pressure types
  array[N] int<lower=1,upper=J> jj; // species
  array[N] int<lower=1,upper=K> kk; // pressure
  array[N] int<lower=1,upper=JK> jk; // species x pressure
  vector[N] y1; // pres_calib
  vector[N] y2; // tens_calib
}

transformed data{
  vector[N] log_y2;
  log_y2 = log(y2);
}

parameters{
  real mu_hat;
  real<lower=0> sigma;
  vector[J] alpha_raw;
  vector[K] beta_raw;
  vector[JK] gamma_raw;
  vector<lower=0,upper=pi()/2>[3] tau_unif;
}

transformed parameters{
  vector[J] alpha;
  vector[K] beta;
  vector[JK] gamma;
  vector<lower=0>[3] tau;
  for (i in 1:3) tau[i] = 2.5 * tan(tau_unif[i]);
  alpha = alpha_raw * tau[1];
  beta = beta_raw * tau[2];
  gamma = gamma_raw * tau[3];
}

model {
  vector[N] mu;
  to_vector(alpha_raw) ~ std_normal();
  to_vector(beta_raw) ~ std_normal();
  to_vector(gamma_raw) ~ std_normal();
  mu = mu_hat + alpha[jj] + beta[kk] + gamma[jk] + log_y2;
  y1 ~ gamma(exp(mu), sigma);
}

generated quantities {
  vector[N] log_lik;
  vector[JK] pred;
  vector[JK] effect;
  for (n in 1:N) {
    pred[jk[n]] = exp(mu_hat + alpha[jj[n]] + beta[kk[n]] + gamma[jk[n]] + log_y2[n]);
    effect[jk[n]] = alpha[jj[n]] + beta[kk[n]] + gamma[jk[n]];
    log_lik[n] = gamma_lpdf(y1[n] | pred[jk[n]], sigma);
  }
}

