data {
  int<lower=0> N; // number of data
  int<lower=0> J; // number of tree samples
  int<lower=0> K; // number of species
  int<lower=0> T; // number of trait predcitors
  array[N] int<lower=1, upper=J> jj; // tree sample
  vector[N] y;
  matrix[N, 2] x;
  matrix[J, T] xj; // trait - segment
  matrix[K, J] uj;
  matrix[T, K] uk;
}

parameters {
  real<lower=0> sigma;
  matrix[2*T, T] beta;
  matrix[2, J] zj;
  matrix[2*T, K] zk;
  cholesky_factor_corr[2] L_Omega_j;
  cholesky_factor_corr[2*T] L_Omega_k;
  vector<lower=0>[2] tau_j;
  vector<lower=0>[2*T] tau_k;
}

// simple
transformed parameters {
  // Multiplying to apply hyper-species to species mapping and adding errors
  matrix[2*T, K] alpha = beta * uk + diag_pre_multiply(tau_k, L_Omega_k) * zk;
  // Multiplying to apply species to tree sample mapping
  matrix[T, J] alpha_a = alpha[1:T, ] * uj;
  matrix[T, J] alpha_b = alpha[(T+1):(2*T), ] * uj;

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

  to_vector(zk) ~ std_normal();
  to_vector(zj) ~ std_normal();
  tau_k ~ normal(0, 2.5);
  tau_j ~ normal(0, 2.5);
  L_Omega_k ~ lkj_corr_cholesky(2);
  L_Omega_j ~ lkj_corr_cholesky(2);

  to_vector(beta) ~ normal(0, 5);
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
