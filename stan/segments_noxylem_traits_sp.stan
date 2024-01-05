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
  vector[T] beta_a;  // species-level traits effects on a
  vector[T] beta_b;  // species-level traits effects on b
  matrix[2, J] zj_A;
  cholesky_factor_corr[2] L_Omega_A;
  vector<lower=0>[2] tau_A;
}

transformed parameters {
  vector[K] a_species;
  vector[K] b_species;
  for (k in 1:K) {
    a_species[k] = dot_product(xk[k], beta_a);
    b_species[k] = dot_product(xk[k], beta_b);
  }
  vector[J] a_hat;  // this will store the a_hat for each segment
  vector[J] b_hat;  // this will store the b_hat for each segment
  for (j in 1:J) {
    a_hat[j] = a_species[kk[j]];
    b_hat[j] = b_species[kk[j]];
  }
  matrix[2, J] A_hat = append_row(to_row_vector(a_hat), to_row_vector(b_hat));
  matrix[2, J] A = A_hat + diag_pre_multiply(tau_A, L_Omega_A) * zj_A;

  matrix[2, K] A_species; // 2 rows for 'a' and 'b', K columns for each species
  vector[K] count; // To count the number of segments per species

  // Initialize the matrix for species-level effects and counts to zero
  for (k in 1:K) {
    A_species[1, k] = 0;
    A_species[2, k] = 0;
    count[k] = 0;
  }

  // Sum the segment-level effects for each species
  for (j in 1:J) {
    int k = kk[j]; // Species index for segment j
    A_species[1, k] += A[1, j]; // Sum 'a' effects for the species
    A_species[2, k] += A[2, j]; // Sum 'b' effects for the species
    count[k] += 1; // Increment count for the species
  }

  // Average the effects for each species
  for (k in 1:K) {
    if (count[k] > 0) {
      A_species[1, k] /= count[k]; // Average 'a' effects for the species
      A_species[2, k] /= count[k]; // Average 'b' effects for the species
    }
  }

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
