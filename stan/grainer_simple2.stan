data {
  int<lower=0> N; // number of data
  int<lower=0> J; // number of segments
  array[N] int<lower=1,upper=J> jj; // segements
  vector[N] y;
  matrix[N, 2] x;
  matrix[1, J] u;
}

parameters {
  real<lower=0> sigma;
  matrix[2, 1] alpha;
  matrix[2, J] z;
  cholesky_factor_corr[2] L_Omega;
  vector<lower=0,upper=pi()/2>[2] tau;
}

transformed parameters {
  matrix[2, J] A = alpha * u + diag_pre_multiply(tau, L_Omega) * z;
}

model {
 vector[N] mu;
  for(n in 1:N) {
    mu[n] = x[n, ] * A[, jj[n]];
  }
  tau ~ normal(0, 2.5);
  to_vector(z) ~ std_normal();
  L_Omega ~ lkj_corr_cholesky(2);
  to_vector(alpha) ~ normal(0, 5);
  y ~ normal(mu, sigma);
}

