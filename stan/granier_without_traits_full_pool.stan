data {
  int<lower=0> N; // number of data
  int<lower=0> K; // number of species
  int<lower=0> L; // number of xylem types
  array[N] int<lower=1,upper=K> kk; // species
  array[N] int<lower=1,upper=L> ll; // xylem type
  vector[N] y;
  matrix[N, 2] x;
  matrix[L, K] uk;
  matrix[1, L] ul;
}

parameters {
  real<lower=0> sigma;
  real log_a_tilde;
  real b_tilde;
  matrix[2, K] zk;
  matrix[2, L] zl;
  cholesky_factor_corr[2] L_Omega_k;
  cholesky_factor_corr[2] L_Omega_l;
  vector<lower=0>[2] tau_k;
  vector<lower=0>[2] tau_l;
}

transformed parameters {
  real log_a = 4.78 + 2.5 * log_a_tilde;
  real b = 1.23 + 2.5 * b_tilde;
  matrix[2, 1] gamma;
  gamma[1, 1] = log_a;
  gamma[2, 1] = b;
  matrix[2, L] beta = gamma * ul + diag_pre_multiply(tau_l, L_Omega_l) * zl;
  matrix[2, K] alpha = beta * uk + diag_pre_multiply(tau_k, L_Omega_k) * zk;
}

model {
 vector[N] mu;
  for(n in 1:N) {
    mu[n] = x[n, ] * alpha[, kk[n]];
  }
  tau_k ~ normal(0, 2.5);
  tau_l ~ normal(0, 2.5);
  to_vector(zk) ~ std_normal();
  to_vector(zl) ~ std_normal();
  L_Omega_k ~ lkj_corr_cholesky(2);
  L_Omega_l ~ lkj_corr_cholesky(2);
  log_a_tilde ~ std_normal();
  b_tilde ~ std_normal();
  y ~ normal(mu, sigma);
}

