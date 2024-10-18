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
  real log_a_tilde;
  real b_tilde;
  matrix[2, J] z;
  cholesky_factor_corr[2] L_Omega;
  vector<lower=0>[2] tau;
}

transformed parameters {
  matrix[2, 1] alpha;
  real log_a = 4.78 + 2.5 * log_a_tilde;
  real b = 1.23 + 2.5 * b_tilde;
  alpha[1, 1] = log_a;
  alpha[2, 1] = b;
  matrix[2, J] A = alpha * u + diag_pre_multiply(tau, L_Omega) * z;
}

model {
 vector[N] mu;
  for(n in 1:N) {
    mu[n] = x[n, ] * A[, jj[n]];
  }
  tau ~ normal(0, 2.5);
  to_vector(z) ~ std_normal();
  log_a_tilde ~ std_normal();
  b_tilde ~ std_normal();
  L_Omega ~ lkj_corr_cholesky(2);
  y ~ normal(mu, sigma);
}

generated quantities {
  real a = exp(log_a + pow(sigma, 2) / 2);
}
