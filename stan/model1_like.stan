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
  sigma ~ normal(0, 2.5);
  log_a_tilde ~ std_normal();
  b_tilde ~ std_normal();
  log_fd ~ normal(log_a + b * log_k, sigma);
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = normal_lpdf(log_fd[n] | log_a + b * log_k[n], sigma);
  }
}
