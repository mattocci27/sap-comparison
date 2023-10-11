data {
  int<lower=0> N; // number of data
  int<lower=0> J; // number of tree samples
  int<lower=0> K; // number of species
  array[N] int<lower=1,upper=J> jj; // tree sample
  array[N] int<lower=1,upper=K> kk; // species
  vector[N] y;
  matrix[N, 2] x;
  matrix[K, J] uj;
  matrix[1, K] uk;
}

parameters {
  real<lower=0> sigma;
  real log_a_tilde;
  real b_tilde;
  matrix[2, J] zj;
  matrix[2, K] zk;
  cholesky_factor_corr[2] L_Omega_j;
  cholesky_factor_corr[2] L_Omega_k;
  vector<lower=0>[2] tau_j;
  vector<lower=0>[2] tau_k;
}

transformed parameters {
  real log_a = 4.78 + 2.5 * log_a_tilde;
  real b = 1.23 + 2.5 * b_tilde;
  matrix[2, 1] gamma;
  gamma[1, 1] = log_a;
  gamma[2, 1] = b;
  matrix[2, K] alpha = gamma * uk + diag_pre_multiply(tau_k, L_Omega_k) * zk;
  matrix[2, J] A = alpha * uj + diag_pre_multiply(tau_j, L_Omega_j) * zj;
  vector[N] mu;
    for(n in 1:N) {
      mu[n] = x[n, ] * A[, jj[n]];
   }
}

model {
  tau_j ~ normal(0, 2.5);
  tau_k ~ normal(0, 2.5);
  to_vector(zk) ~ std_normal();
  to_vector(zj) ~ std_normal();
  L_Omega_k ~ lkj_corr_cholesky(2);
  L_Omega_j ~ lkj_corr_cholesky(2);
  log_a_tilde ~ std_normal();
  b_tilde ~ std_normal();
  y ~ normal(mu, sigma);
}

generated quantities {
  vector[N] log_lik;
  cov_matrix[2] Sigma_j = diag_pre_multiply(tau_j, L_Omega_j) *
    diag_post_multiply(L_Omega_j', tau_j);
  cov_matrix[2] Sigma_k = diag_pre_multiply(tau_k, L_Omega_k) *
    diag_post_multiply(L_Omega_k', tau_k);
  real<lower=-1, upper=1> rho_j = Sigma_j[1, 2] * inv(tau_j[1] * tau_j[2]);
  real<lower=-1, upper=1> rho_k = Sigma_k[1, 2] * inv(tau_k[1] * tau_k[2]);
  vector[K] a = exp(to_vector(alpha[1, ]) + pow(sigma, 2) / 2);
  for (n in 1:N) {
    log_lik[n] = normal_lpdf(y[n] | mu[n], sigma);
  }
}
