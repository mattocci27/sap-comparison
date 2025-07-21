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
  matrix[K, J] uj;
  matrix[L, K] uk;
  matrix[1, L] ul;
}

parameters {
  real<lower=0> sigma;
  vector[4] gamma;
  matrix[2, J] zj;
  matrix[4, K] zk;
  matrix[4, L] zl;
  cholesky_factor_corr[2] L_Omega_j;
  cholesky_factor_corr[4] L_Omega_k;
  cholesky_factor_corr[4] L_Omega_l;
  vector<lower=0>[2] tau_j;
  vector<lower=0>[4] tau_k;
  vector<lower=0>[4] tau_l;
}

// simple
transformed parameters {
  // Multiplying to apply hyper-xylem to xylem mapping and adding errors
  matrix[4, L] beta = gamma * to_row_vector(ul) + diag_pre_multiply(tau_l, L_Omega_l) * zl;
  // Multiplying to apply hyper-species to species mapping and adding errors
  matrix[4, K] alpha = beta * uk + diag_pre_multiply(tau_k, L_Omega_k) * zk;
  // Multiplying to apply species to tree sample mapping
  matrix[T, J] alpha_a = alpha[1:2, ] * uj;
  matrix[T, J] alpha_b = alpha[3:4, ] * uj;

  // Calculate the final effects per tree sample, incorporating trait predictors
  row_vector[J] mu_a = to_row_vector(diagonal(xj * alpha_a));
  row_vector[J] mu_b = to_row_vector(diagonal(xj * alpha_b));

  matrix[2, J] A = to_matrix(append_row(mu_a, mu_b)) + diag_pre_multiply(tau_j, L_Omega_j) * zj;

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
