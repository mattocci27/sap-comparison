data {
  int<lower=0> N;  // number of observations
  int<lower=0> J;  // number of tree samples
  int<lower=0> K;  // number of species
  int<lower=0> T;  // number of tree‐level trait predictors
  array[N] int<lower=1, upper=J> jj; // tree‐sample index for each obs
  vector[N]               y;        // outcome: log(fd)
  matrix[N, 2]            x;        // (intercept, log(k)) per observation
  matrix[J, T] xj; // trait -sp
  matrix[K, J] uj;
}

parameters {
  real<lower=0> sigma;
  // (1) Species‐level random effects: 2 rows × K species
  matrix[2, K] theta;
  matrix[2, K] zk;
  cholesky_factor_corr[2] L_Omega_k;
  vector<lower=0>[2] tau_k;
  // (2) Tree‐level trait slopes: a 2×T matrix B
  matrix[2, T] B; // one row for “a”, one row for “b”
  matrix[2, J] zj;
  cholesky_factor_corr[2] L_Omega_j;
  vector<lower=0>[2]      tau_j;
}

transformed parameters {
  // STEP 1: Build the 2×K species‐level coefficient matrix θ_centered:
  //   - Non‐centered decomposition: θ = (τ_k × L_Omega_k) × zk
  //──────────────────────────────────────────────────────────────────────────────
  matrix[2, K] theta_centered = diag_pre_multiply(tau_k, L_Omega_k) * zk;
  // STEP 2: Map species → tree to get a 2×J matrix of “species‐level intercept & slope” for each tree sample:
  //   α_species_tree = θ_centered (2×K) × uj (K×J)  → (2×J)
  //
  matrix[2, J] alpha_species_tree;
  alpha_species_tree = theta_centered * uj;
  //                      (2×K) × (K×J)  = (2×J)

  // STEP 3: Build the 2×J “tree‐level trait adjustment” from the T traits:
  //   • B is 2×T, xj' is T×J → B × xj' is 2×J
  //
  matrix[2, J] alpha_trait_tree;
  alpha_trait_tree = B * xj';
  //                  (2×T) × (T×J) = (2×J)

  // STEP 4: Now combine species‐level + trait‐level + tree‐noise to get final 2×J “A”:
  //
  matrix[2, J] A = alpha_species_tree
                 + alpha_trait_tree
                 + diag_pre_multiply(tau_j, L_Omega_j) * zj;
  //   (2×J) +   (2×J)  +   (2×2)×(2×J)  = (2×J)

  // ─────────────────────────────────────────────────────────────────────────────
  // That 2×J matrix “A” is:
  //   A[1, j] = (species intercept for j’s species)
  //             + row1_of_B × xj[j,·]
  //             + tree‐noise_intercept[j]
  //   A[2, j] = (species slope for j’s species)
  //             + row2_of_B × xj[j,·]
  //             + tree‐noise_slope[j]
  //
  // In the model below, we’ll use A[, jj[n]] to find which (intercept, slope) to multiply
  // against x[n, ] for observation n.
  // ─────────────────────────────────────────────────────────────────────────────
}

model {
  //──────────────────────────────────────────────────────────────────────────────
  // (1) Priors on non‐centered species‐level errors:
  to_vector(zk) ~ std_normal();

  // (2) Priors on species‐level scales + correlation:
  tau_k     ~ normal(0, 2.5);
  L_Omega_k ~ lkj_corr_cholesky(2);

  // (3) Priors on B (tree‐level trait slopes):
  to_vector(B) ~ normal(0, 5);

  // (4) Priors on tree‐level noise:
  to_vector(zj) ~ std_normal();
  tau_j     ~ normal(0, 2.5);
  L_Omega_j ~ lkj_corr_cholesky(2);

  // (5) Likelihood:  for each observation n,
  //      μ[n] = x[n, ] (1×2) × A[, jj[n]] (2×1)
  vector[N] mu;
  for (n in 1:N) {
    mu[n] = x[n, ] * A[, jj[n]];
  }
  y ~ normal(mu, sigma);
}

generated quantities {
  vector[N] log_lik;
  vector[N] mu_rep;
  for (n in 1:N) {
    mu_rep[n]  = x[n, ] * A[, jj[n]];
    log_lik[n] = normal_lpdf(y[n] | mu_rep[n], sigma);
  }
}
