data {
  int<lower=0> N; // number of data
  vector[N] log_fd;
  vector[N] log_k;
}

parameters {
  real<lower=0> sigma;
  real log_a;
  real b;
}

model {
  vector[N] mu;
  sigma ~ normal(0, 2.5);
  log_a ~ normal(0, 5);
  b ~ normal(0, 5);
  log_fd ~ normal(log_a + b * log_k, sigma);
}

