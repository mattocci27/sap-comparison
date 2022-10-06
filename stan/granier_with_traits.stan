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
  matrix[T, 1] gamma_a;
  matrix[T, 1] gamma_b;
  matrix[2, J] zj_A;
  matrix[T, J] zj_a;
  matrix[T, K] zk_a;
  matrix[T, L] zl_a;
  matrix[T, J] zj_b;
  matrix[T, K] zk_b;
  matrix[T, L] zl_b;
  cholesky_factor_corr[2] L_Omega_A;
  cholesky_factor_corr[T] L_Omega_j_a;
  cholesky_factor_corr[T] L_Omega_k_a;
  cholesky_factor_corr[T] L_Omega_l_a;
  cholesky_factor_corr[T] L_Omega_j_b;
  cholesky_factor_corr[T] L_Omega_k_b;
  cholesky_factor_corr[T] L_Omega_l_b;
  vector<lower=0>[2] tau_A;
  vector<lower=0>[T] tau_j_a;
  vector<lower=0>[T] tau_k_a;
  vector<lower=0>[T] tau_l_a;
  vector<lower=0>[T] tau_j_b;
  vector<lower=0>[T] tau_k_b;
  vector<lower=0>[T] tau_l_b;
}

transformed parameters {
  matrix[T, L] beta_a = gamma_a * ul + diag_pre_multiply(tau_l_a, L_Omega_l_a) * zl_a;
  matrix[T, K] alpha_a = beta_a * uk + diag_pre_multiply(tau_k_a, L_Omega_k_a) * zk_a;
  vector[J] a_hat = diagonal(xj * (alpha_a * uj + diag_pre_multiply(tau_j_a, L_Omega_j_a) * zj_a));
  matrix[T, L] beta_b = gamma_b * ul + diag_pre_multiply(tau_l_b, L_Omega_l_b) * zl_b;
  matrix[T, K] alpha_b = beta_b * uk + diag_pre_multiply(tau_k_b, L_Omega_k_b) * zk_b;
  vector[J] b_hat = diagonal(xj * (alpha_b * uj + diag_pre_multiply(tau_j_b, L_Omega_j_b) * zj_b));
  matrix[2, J] A_hat = append_row(to_row_vector(a_hat), to_row_vector(b_hat));
  matrix[2, J] A = A_hat + diag_pre_multiply(tau_A, L_Omega_A) * zj_A;
}

model {
  vector[N] mu;
  for (n in 1:N) {
    mu[n] = x[n, ] * A[, jj[n]];
  }

  to_vector(zl_a) ~ std_normal();
  to_vector(zk_a) ~ std_normal();
  to_vector(zj_a) ~ std_normal();
  tau_l_a ~ normal(0, 2.5);
  tau_k_a ~ normal(0, 2.5);
  tau_j_a ~ normal(0, 2.5);
  L_Omega_l_a ~ lkj_corr_cholesky(2);
  L_Omega_k_a ~ lkj_corr_cholesky(2);
  L_Omega_j_a ~ lkj_corr_cholesky(2);

  to_vector(zl_b) ~ std_normal();
  to_vector(zk_b) ~ std_normal();
  to_vector(zj_b) ~ std_normal();
  tau_l_b ~ normal(0, 2.5);
  tau_k_b ~ normal(0, 2.5);
  tau_j_b ~ normal(0, 2.5);
  L_Omega_l_b ~ lkj_corr_cholesky(2);
  L_Omega_k_b ~ lkj_corr_cholesky(2);
  L_Omega_j_b ~ lkj_corr_cholesky(2);

  to_vector(zj_A) ~ std_normal();
  to_vector(gamma_a) ~ normal(0, 5);
  to_vector(gamma_b) ~ normal(0, 5);
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
