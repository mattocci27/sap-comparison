
data {
  int<lower=0> N; // number of data
  int<lower=0> J; // number of tree samples
  int<lower=0> K; // number of species
  int<lower=0> L; // number of xylem types
  array[N] int<lower=1,upper=J> jj; // tree sample
  array[N] int<lower=1,upper=K> kk; // species
  array[N] int<lower=1,upper=L> ll; // xylem type
  vector[N] y;
  matrix[N,2] x;
  matrix[K,J] uj;
  matrix[L,K] uk;
  matrix[1,L] ul;
}

parameters {
  real<lower=0> sigma;
  matrix[2, 1] gamma;
  matrix[2, J] zj;
  matrix[2, K] zk;
  matrix[2, L] zl;
  cholesky_factor_corr[2] L_Omega_j;
  cholesky_factor_corr[2] L_Omega_k;
  cholesky_factor_corr[2] L_Omega_l;
  vector<lower=0,upper=pi()/2>[2] tau_unif_j;
  vector<lower=0,upper=pi()/2>[2] tau_unif_k;
  vector<lower=0,upper=pi()/2>[2] tau_unif_l;
}

transformed parameters {
  vector<lower=0>[2] tau_j = 2.5 * tan(tau_unif_j);
  vector<lower=0>[2] tau_k = 2.5 * tan(tau_unif_k);
  vector<lower=0>[2] tau_l = 2.5 * tan(tau_unif_l);
  matrix[2, L] beta = gamma * ul + diag_pre_multiply(tau_l, L_Omega_l) * zl;
  matrix[2, K] alpha = beta * uk + diag_pre_multiply(tau_k, L_Omega_k) * zk;
  matrix[2, J] A = alpha * uj + diag_pre_multiply(tau_j, L_Omega_j) * zj;
}

model {
 vector[N] mu;
  for(n in 1:N) {
    mu[n] = x[n, ] * A[, jj[n]];
  }
  to_vector(zl) ~ std_normal();
  to_vector(zk) ~ std_normal();
  to_vector(zj) ~ std_normal();
  L_Omega_l ~ lkj_corr_cholesky(2);
  L_Omega_k ~ lkj_corr_cholesky(2);
  L_Omega_j ~ lkj_corr_cholesky(2);
  to_vector(gamma) ~ normal(0, 5);
  y ~ normal(mu, sigma);
}

