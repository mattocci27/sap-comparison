// -----------------------------------------------------------------------------
//  segments_noxylem_traits_sp_simple_all.stan
//
//  – Uses tree‐level predictors (xj) to explain tree‐level intercept (Co‐a) & slope (Co‐b).
//  – Keeps a species‐level random‐effects layer (θ_k) with no species‐level covariates.
//
//  Data inputs:
//
//    N   = total # of log(fd) observations across all trees
//    J   = # of distinct tree‐samples (each “j” is one tree)
//    K   = # of distinct species
//    T   = # of tree‐level traits (e.g. VAF, SWC, etc.)
//
//    jj[n]      = integer in {1,…,J}, “which tree sample does obs n belong to?”
//    y[n]       = observed log(fd)[n]
//    x[n, 1:2]  = [1,  log(k) ] for each obs n  (the “intercept” and the log(k) predictor in the final regression)
//    xj[j, 1:T] = T tree‐level covariates for tree sample j  (e.g. log(VAF_j), log(SWC_j), …)
//    uj[k, j]   = 1 if tree j belongs to species k; 0 otherwise  (i.e. a K×J one‐hot “species → tree” map)
//
//  Parameters:
//
//    – Species‐level random effects (θ_k) for (Co‐a, Co‐b) of each species k:
//        θ_centered[k,1:2] = diag_pre_multiply(τ_k, L_Omega_k) · zk[k,1:2]
//        → these are then mapped to each tree by uj[k,j] to give α_species_tree[1:2,j].
//
//    – Tree‐level trait slopes B[1:2, 1:T]:
//        That multiplies xj′ to give α_trait_tree[1:2, j] = B[1:2,·] · xj[j,·]′.
//
//    – Tree‐level noise zj[1:2, j], with covariance τ_j and L_Omega_j.
//
//  Final “A” for tree j is a 2×1 vector = [Co‐a; Co‐b] for that tree, used in the regression
//     μ[n] = x[n, ] (1×2) · A[1:2, jj[n]] (2×1).
//
// -----------------------------------------------------------------------------
data {
  int<lower=0> N;              // total number of log(fd) observations
  int<lower=0> J;              // number of tree samples
  int<lower=0> K;              // number of species
  int<lower=0> T;              // number of tree‐level trait predictors

  array[N] int<lower=1, upper=J> jj; // which tree each observation n belongs to
  vector[N] y;                      // outcome: log(fd) for each obs

  matrix[N, 2] x;                   // each obs: (1, log(k)) for the final regression

  matrix[J, T] xj;                  // tree‐level traits  (each row j has T covariates)
  matrix[K, J] uj;                  // one‐hot species→tree:  if tree j belongs to species k, uj[k,j] = 1
}

parameters {
  //  (1) Species‐level random effects: 2 rows × K species
  //      θ_centered[k, 1:2] = diag_pre_multiply(τ_k, L_Omega_k) · zk[k, 1:2]
  vector<lower=0>[2]    tau_k;          // species‐level SDs (Co‐a, Co‐b)
  cholesky_factor_corr[2] L_Omega_k;    // 2×2 Cholesky for species correlation
  matrix[2, K]         zk;              // non‐centered species noise (2×K standard normals)

  //  (2) Tree‐level trait slopes: a 2×T matrix B
  //      Row 1 = slope for Co‐a vs. each of the T traits
  //      Row 2 = slope for Co‐b vs. each of the T traits
  matrix[2, T] B;                       // 2×T

  //  (3) Tree‐level random noise: 2 rows × J trees
  vector<lower=0>[2]     tau_j;         // tree‐level SDs (Co‐a noise, Co‐b noise)
  cholesky_factor_corr[2] L_Omega_j;    // 2×2 Cholesky for tree correlation
  matrix[2, J]          zj;             // non‐centered tree noise (2×J standard normals)

  //  (4) Population‐level residual for each observation n in y[n] ~ Normal(mu[n], sigma):
  real<lower=0> sigma;                  // residual SD
}

transformed parameters {
  // STEP 1: Build species‐level “centered” random effects θ_centered[k,1:2]
  //         (2×K = diag_pre_multiply(τ_k, L_Omega_k) × zk).
  matrix[2, K] theta_centered;
  theta_centered = diag_pre_multiply(tau_k, L_Omega_k) * zk;
  //    → for each species k, θ_centered[1,k] = “species intercept (Co‐a)”;
  //      θ_centered[2,k] = “species slope (Co‐b)”.

  // STEP 2: Map species‐level random effects to each tree j:
  //   α_species_tree[1:2, j] = sum_{k=1}^K [ θ_centered[1:2, k] * uj[k,j] ].
  //
  matrix[2, J] alpha_species_tree;
  alpha_species_tree = theta_centered * uj;
  //   dims: (2×K) × (K×J) = 2×J
  //
  //   → For tree j, alpha_species_tree[1,j] is “species‐intercept contribution”,
  //     alpha_species_tree[2,j] is “species‐slope contribution”.

  // STEP 3: Build tree‐level “trait adjustment” from xj:
  //   α_trait_tree[1:2, j] = B[1:2, 1:T] × xj[j, 1:T]′.
  //   i.e. each row of B is a 1×T row, each row of xj is a 1×T row → result 2×J.
  //
  matrix[2, J] alpha_trait_tree;
  alpha_trait_tree = B * xj';
  //   dims: (2×T) × (T×J) = 2×J
  //   → For tree j, alpha_trait_tree[1,j] is trait‐contribution to Co‐a,
  //     and alpha_trait_tree[2,j] is trait‐contribution to Co‐b.

  // STEP 4: Build tree‐level noise = (2×2) × (2×J) = 2×J
  //   We do diag_pre_multiply(τ_j, L_Omega_j) * zj  → (2×J) per draw.
  matrix[2, J] alpha_tree_noise;
  alpha_tree_noise = diag_pre_multiply(tau_j, L_Omega_j) * zj;

  // STEP 5: Final 2×J “A” = [ Co‐a_j ; Co‐b_j ] for each tree j:
  //
  matrix[2, J] A;
  A = alpha_species_tree
    + alpha_trait_tree
    + alpha_tree_noise;
  //   dims: (2×J)
  //   → A[1,j] = final tree‐intercept (Co‐a) for tree j
  //      A[2,j] = final tree‐slope     (Co‐b) for tree j
}

model {
  // (1) Species‐level non‐centered standard normals:
  to_vector(zk) ~ std_normal();

  // (2) Priors on species‐level SDs and correlation:
  tau_k     ~ normal(0, 2.5);
  L_Omega_k ~ lkj_corr_cholesky(2);

  // (3) Priors on tree‐level trait slopes B:
  to_vector(B) ~ normal(0, 5);

  // (4) Tree‐level non‐centered standard normals:
  to_vector(zj) ~ std_normal();

  // (5) Priors on tree‐level SDs and correlation:
  tau_j     ~ normal(0, 2.5);
  L_Omega_j ~ lkj_corr_cholesky(2);

  // (6) Residual SD:
  sigma ~ normal(0, 2.5);

  // (7) Likelihood: each observation n selects which tree j = jj[n],
  //     then μ[n] = x[n, ] (1×2) · A[1:2, j] (2×1) = intercept_j + k‐slope_j * log(k_n).
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
