data {
  int<lower=0> N; // number of data
  int<lower=0> J; // number of tree samples
  int<lower=0> K; // number of species
  int<lower=0> T; // number of trait predcitors
  array[N] int<lower=1, upper=J> jj; // tree sample
  array[J] int<lower=1, upper=K> kk;  // Mapping of segments to species
  vector[N] y;
  matrix[N, 2] x;
  matrix[J, T] xj; // trait - segment
  matrix[K, J] uj;
  matrix[1, K] uk;
  matrix[K, T] xk; // trait - species level
}

parameters {
  real<lower=0> sigma;
  vector[2] alpha;
  vector[2] beta;
  real<lower=0> sig_alpha;
  real<lower=0> sig_beta;
  matrix[T, 1] beta_a;
  matrix[T, 1] beta_b;
  matrix[2, J] zj_A;
  matrix[T, J] zj_a;
  matrix[T, K] zk_a;
  matrix[T, J] zj_b;
  matrix[T, K] zk_b;
  cholesky_factor_corr[2] L_Omega_A;
  cholesky_factor_corr[T] L_Omega_j_a;
  cholesky_factor_corr[T] L_Omega_k_a;
  cholesky_factor_corr[T] L_Omega_j_b;
  cholesky_factor_corr[T] L_Omega_k_b;
  vector<lower=0>[2] tau_A;
  vector<lower=0>[T] tau_j_a;
  vector<lower=0>[T] tau_k_a;
  vector<lower=0>[T] tau_j_b;
  vector<lower=0>[T] tau_k_b;
}

transformed parameters {
  matrix[T, K] alpha_a = beta_a * uk + diag_pre_multiply(tau_k_a, L_Omega_k_a) * zk_a;
  matrix[T, K] alpha_b = beta_b * uk + diag_pre_multiply(tau_k_b, L_Omega_k_b) * zk_b;
  matrix[T, J] species_adjusted_alpha_a = alpha_a * uj;
  matrix[T, J] random_effects_a = diag_pre_multiply(tau_j_a, L_Omega_j_a) * zj_a;
  matrix[T, J] species_adjusted_alpha_b = alpha_b * uj;
  matrix[T, J] random_effects_b = diag_pre_multiply(tau_j_b, L_Omega_j_b) * zj_b;
  vector[J] a_hat;
  vector[J] b_hat;
    for (j in 1:J) {
        a_hat[j] = dot_product(xj[j], species_adjusted_alpha_a[, j] + random_effects_a[, j]);
        b_hat[j] = dot_product(xj[j], species_adjusted_alpha_b[, j] + random_effects_b[, j]);
    }
  matrix[2, J] A_hat = append_row(to_row_vector(a_hat), to_row_vector(b_hat));
  matrix[2, J] A = A_hat + diag_pre_multiply(tau_A, L_Omega_A) * zj_A;
  vector[N] mu;
  for (n in 1:N) {
    mu[n] = x[n, ] * A[, jj[n]];
  }
  // species-level a and b
  vector[K] species_a; //
  vector[K] species_b; //
  vector[K] count; // Counts of segments for each species
  // Initialize species-level effects and counts to zero
  for (k in 1:K) {
    species_a[k] = 0;
    species_b[k] = 0;
    count[k] = 0;
  }

  // Aggregate segment-level effects for each species
  for (j in 1:J) {
    int k = kk[j]; // Species index for segment j
    species_a[k] += a_hat[j]; // Add the effect a for the species
    species_b[k] += b_hat[j]; // Add the effect b for the species
    count[k] += 1; // Increment count for the species
  }

  // Calculate mean effects for each species
  for (k in 1:K) {
    if (count[k] > 0) {
      species_a[k] /= count[k]; // Calculate mean effect a
      species_b[k] /= count[k]; // Calculate mean effect b
    }
    // Else, you might want to handle the case where count[k] == 0, if necessary
  }
}

model {
  to_vector(zk_a) ~ std_normal();
  to_vector(zj_a) ~ std_normal();
  tau_k_a ~ normal(0, 2.5);
  tau_j_a ~ normal(0, 2.5);
  L_Omega_k_a ~ lkj_corr_cholesky(2);
  L_Omega_j_a ~ lkj_corr_cholesky(2);

  to_vector(zk_b) ~ std_normal();
  to_vector(zj_b) ~ std_normal();
  tau_k_b ~ normal(0, 2.5);
  tau_j_b ~ normal(0, 2.5);
  L_Omega_k_b ~ lkj_corr_cholesky(2);
  L_Omega_j_b ~ lkj_corr_cholesky(2);

  to_vector(zj_A) ~ std_normal();
  to_vector(beta_a) ~ normal(0, 5);
  to_vector(beta_b) ~ normal(0, 5);
  y ~ normal(mu, sigma);

  species_a ~ normal(alpha[1] + alpha[2] * to_vector(xk[,1]), sig_alpha);
  species_b ~ normal(beta[1] + beta[2] * to_vector(xk[,2]), sig_beta);

}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = normal_lpdf(y[n] | mu[n], sigma);
  }
}
