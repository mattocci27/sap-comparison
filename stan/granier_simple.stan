data {
  int<lower=0> N; // number of data
  vector[N] log_fd;
  vector[N] log_k;
}

parameters {
  real<lower=0> sigma;
  real log_a_tilde;
  real b_tilde;
}

transformed parameters {
  real log_a = 4.78 + 2.5 * log_a_tilde;
  real b = 1.23 + 2.5 * b_tilde;
}

model {
  vector[N] mu;
  sigma ~ normal(0, 2.5);
  log_a_tilde ~ std_normal();
  b_tilde ~ std_normal();
  log_fd ~ normal(log_a + b * log_k, sigma);
}

generated quantities {
  matrix[2, 1] alpha;
  alpha[1, 1] = log_a;
  alpha[2, 1] = b;
}
