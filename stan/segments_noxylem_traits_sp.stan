data {
  int<lower=0> N; // number of data
  int<lower=0> J; // number of tree samples
  int<lower=0> K; // number of species
  int<lower=0> T; // number of trait predcitors
  array[N] int<lower=1,upper=J> jj; // tree sample
  array[J] int<lower=1, upper=K> kk;  // Mapping of segments to species
  vector[N] y;
  matrix[N, 2] x;
  matrix[K, T] xk; // trait - species level
}

parameters {
  real<lower=0> sigma;
  matrix[T, 1] beta_a;  // species-level traits effects on a
  matrix[T, 1] beta_b;  // species-level traits effects on b
  matrix[2, J] zj_A;
  cholesky_factor_corr[2] L_Omega_A;
  vector<lower=0>[2] tau_A;
}

transformed parameters {
  vector[K] a_species = to_vector(xk * beta_a);  // species-level a for each species
  vector[J] a_hat;  // this will store the a_hat for each segment
  vector[K] b_species = to_vector(xk * beta_b);  // species-level b for each species
  vector[J] b_hat;  // this will store the b_hat for each segment
  for (j in 1:J) {
    a_hat[j] = a_species[kk[j]];
    b_hat[j] = b_species[kk[j]];
  }
  matrix[2, J] A_hat = append_row(to_row_vector(a_hat), to_row_vector(b_hat));
  matrix[2, J] A = A_hat + diag_pre_multiply(tau_A, L_Omega_A) * zj_A;
  vector[N] mu;
  for (n in 1:N) {
    mu[n] = x[n, ] * A[, jj[n]];
  }
}

model {
  to_vector(beta_a) ~ normal(0, 5);
  to_vector(beta_b) ~ normal(0, 5);
  to_vector(zj_A) ~ std_normal();
  tau_A ~ normal(0, 2.5);
  L_Omega_A ~ lkj_corr_cholesky(2);
  y ~ normal(mu, sigma);
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = normal_lpdf(y[n] | mu[n], sigma);
  }
}
