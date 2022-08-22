data {
  int<lower=1> N; // number of samples
  int<lower=1> J; // number of species
  int<lower=1> K; // number of pressure types
  int<lower=1> JK; // number of species x pressure types
  int<lower=1> L; // number of methods (pres or tens)
  int<lower=1> M; // number of branches
  array[N] int<lower=1,upper=J> jj; // species
  array[N] int<lower=1,upper=K> kk; // pressure
  array[N] int<lower=1,upper=JK> jk; // species x pressure
  array[N] int<lower=1,upper=L> ll; // pres or tens
  array[N] int<lower=1,upper=M> mm; // branch ID
  int<lower=0> model_type; // normal, log-nomral, gamma
  vector[N] y; // pres_calib and tens_calib
}

transformed data {
 vector[N] log_y;
  log_y = log(y);
}

parameters{
  vector[L] mu_hat;
  real<lower=0> sigma;
  matrix[L, J] alpha;
  matrix[L, K] beta;
  matrix[L, JK] gamma;
  matrix[L, M] delta_raw;
  cholesky_factor_corr[L] L_Omega_4;
  vector<lower=0,upper=pi()/2>[L] tau_unif_4;
}


transformed parameters{
  matrix[L, M] delta;
  vector<lower=0>[L] tau_4;
  for (i in 1:L) {
    tau_4[i] = 2.5 * tan(tau_unif_4[i]);
  }
  delta = diag_pre_multiply(tau_4, L_Omega_4) * delta_raw;
}

model {
  vector[N] mu;
  to_vector(alpha) ~ normal(0, 5);
  to_vector(beta) ~ normal(0, 5);
  to_vector(gamma) ~ normal(0, 5);
  to_vector(delta_raw) ~ std_normal();
  L_Omega_4 ~ lkj_corr_cholesky(2);
  for (n in 1:N) {
  mu[n] = mu_hat[ll[n]] + alpha[ll[n], jj[n]] + beta[ll[n], kk[n]] + gamma[ll[n], jk[n]] + delta[ll[n], mm[n]];
  }
  if (model_type == 1) {
    y ~ normal(mu, sigma);
  } else if (model_type == 2) {
    log_y ~ normal(mu, sigma);
  } else if (model_type == 3) {
    y ~ gamma(exp(mu), sigma);
  }
}

generated quantities {
  real mu_diff;
  vector[J] alpha_diff;
  vector[K] beta_diff;
  vector[JK] gamma_diff;
  vector[JK] effect_diff;
  matrix[L, JK] effect;
  mu_diff = mu_hat[1] - mu_hat[2];
  for (n in 1:N) {
    effect[ll[n], jk[n]] = alpha[ll[n], jj[n]] + beta[ll[n], kk[n]] + gamma[ll[n], jk[n]];
  }
  for (j in 1:J) alpha_diff[j] = alpha[1, j] - alpha[2, j];
  for (k in 1:K) beta_diff[k] = beta[1, k] - beta[2, k];
  for (i in 1:JK) gamma_diff[i] = gamma[1, i] - gamma[2, i];
  for (i in 1:JK) effect_diff[i] = effect[1, i] - effect[2, i];
}

