data {
  int<lower=0> N; // number of samples
  int<lower=0> J; // number of species
  int<lower=0> K; // number of parameters (int, beta1, beta2, bp)
  array[N] int<lower=1,upper=J> jj; // species
  array[N] int<lower=0> y; // count
  array[N] int<lower=0> total; // total
  matrix[N, K] x; // tree-level predictor
  matrix[1, J] u; // sp-level intercept
}

parameters {
  matrix[K, 1] gamma;
  matrix[K, J] z;
  cholesky_factor_corr[K] L_Omega;
  vector<lower=0,upper=pi()/2>[K] tau_unif;
}

transformed parameters {
  matrix[K, J] beta;
  vector<lower=0>[K] tau;
  for (k in 1:K) tau[k] = 2.5 * tan(tau_unif[k]);
  beta = gamma * u + diag_pre_multiply(tau, L_Omega) * z;
}

model {
  vector[N] alpha;
  for (n in 1:N) {
    alpha[n] = x[n, ] * beta[, jj[n]];
  }
  to_vector(z) ~ std_normal();
  L_Omega ~ lkj_corr_cholesky(2);
  to_vector(gamma) ~ normal(0, 5);
  y ~ binomial_logit(total, alpha);
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = binomial_logit_lpmf(y[n] | total[n], x[n, ] * beta[, jj[n]]);
  }
}
