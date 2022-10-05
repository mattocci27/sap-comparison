data {
  int<lower=0> N; // number of data
  int<lower=0> J; // number of tree samples
  int<lower=0> K; // number of species
  int<lower=0> L; // number of xylem types
  int<lower=0> T; // number of trait predcitors
  array[N] int<lower=1, upper=J> jj; // tree sample
  array[N] int<lower=1, upper=K> kk; // species
  array[N] int<lower=1, upper=L> ll; // xylem type
  // array[N] int<lower=1, upper=T> tt; // traits
  vector[N] y;
  matrix[N, 2] x;
  matrix[J, T] xj; // trait - segment
  matrix[K, J] uj;
  matrix[L, K] uk;
  matrix[1, L] ul;
}

parameters {
  real<lower=0> sigma;
  // real log_a_tilde;
  // real b_tilde;
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
  // real log_a = 4.78 + 2.5 * log_a_tilde;
  // real b  = 1.23 + 2.5 * b_tilde;
  // matrix[2, 1] gamma;
  // gamma[1, 1] = log_a;
  // gamma[2, 1] = b;
  matrix[T, L] beta_a = gamma_a * ul + diag_pre_multiply(tau_l_a, L_Omega_l_a) * zl_a;
  matrix[T, K] alpha_a = beta_a * uk + diag_pre_multiply(tau_k_a, L_Omega_k_a) * zk_a;
  // matrix[J, J] a_tmp = xj * (alpha_a * uj + diag_pre_multiply(tau_j_a, L_Omega_j_a) * zj_a);
  // vector[J] a_hat = diagonal(a_tmp);
  // vector[J] a_hat = diagonal(xj * alpha_a * uj + diag_pre_multiply(tau_j_a, L_Omega_j_a) * zj_a);
  vector[J] a_hat = diagonal(xj * (alpha_a * uj + diag_pre_multiply(tau_j_a, L_Omega_j_a) * zj_a));
  matrix[T, L] beta_b = gamma_b * ul + diag_pre_multiply(tau_l_b, L_Omega_l_b) * zl_b;
  matrix[T, K] alpha_b = beta_b * uk + diag_pre_multiply(tau_k_b, L_Omega_k_b) * zk_b;

  vector[J] b_hat = diagonal(xj * (alpha_b * uj + diag_pre_multiply(tau_j_b, L_Omega_j_b) * zj_b));
  matrix[2, J] A_hat = append_row(to_row_vector(a_hat), to_row_vector(b_hat));
  matrix[2, J] A = A_hat + diag_pre_multiply(tau_A, L_Omega_A) * zj_A;
}

model {
 vector[N] mu;
  for(n in 1:N) {
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

  to_vector(gamma_a) ~ normal(0, 5);
  to_vector(gamma_b) ~ normal(0, 5);
  // log_a_tilde ~ std_normal();
  // b_tilde ~ std_normal();
  y ~ normal(mu, sigma);
}
