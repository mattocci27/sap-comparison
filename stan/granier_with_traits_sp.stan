data {
  int<lower=0> N; // number of data
  int<lower=0> J; // number of tree samples
  int<lower=0> K; // number of species
  int<lower=0> L; // number of xylem types
  int<lower=0> T; // number of trait predcitors
  array[N] int<lower=1, upper=J> jj; // tree sample
  vector[N] y;
  matrix[N, 2] x;
  matrix[J, T] xj; // trait - segment
  matrix[L, T] ul;
  matrix[T, K] uk;
  matrix[K, J] uj;
}

parameters {
  real<lower=0> sigma;
  matrix[2, L] gamma;
  matrix[2, T] z_l;
  matrix[2, K] z_k;
  matrix[2, J] z_j;
  cholesky_factor_corr[2] L_Omega_k;
  cholesky_factor_corr[2] L_Omega_l;
  cholesky_factor_corr[2] L_Omega_j;
  vector<lower=0>[2] tau_k;
  vector<lower=0>[2] tau_l;
  vector<lower=0>[2] tau_j;
}

transformed parameters {
  matrix[2, T] beta = gamma * ul + diag_pre_multiply(tau_l, L_Omega_l) * z_l;
  matrix[2, K] alpha = beta * uk + diag_pre_multiply(tau_k, L_Omega_k) * z_k;
  matrix[2, J] A = alpha * uj + diag_pre_multiply(tau_j, L_Omega_j) * z_j;
}

model {
  vector[N] mu;
  for (n in 1:N) {
    mu[n] = x[n, ] * A[, jj[n]];
  }

  to_vector(zl) ~ std_normal();
  to_vector(zk) ~ std_normal();
  to_vector(zj) ~ std_normal();
  tau_l ~ normal(0, 2.5);
  tau_k ~ normal(0, 2.5);
  tau_j ~ normal(0, 2.5);
  L_Omega_l ~ lkj_corr_cholesky(2);
  L_Omega_k ~ lkj_corr_cholesky(2);
  L_Omega_j ~ lkj_corr_cholesky(2);
  to_vector(gamma) ~ normal(0, 5);
  y ~ normal(mu, sigma);
}

generated quantities {
  vector[N] log_lik;
  vector[N] mu;
  for (n in 1:N) {
    mu[n] = x[n, ] * A[, jj[n]];
    log_lik[n] = normal_lpdf(y[n] | mu[n], sigma);
  }
}
