data {
  int<lower=1> N;  // number of observations
  int<lower=1> K; // number of predictors
  int<lower=1> M;  // number of unique trees
  int<lower=1> T;  // number of unique times
  vector[N] log_ks;  // response variable log(ks)
  vector[N] log_ks_ref;  // offset
  matrix[N, K] x;  // predictor
  array[N] int<lower=1, upper=M> tree; // integer
  array[N] int<lower=1, upper=T> time; // integer
}

parameters {
  vector[K] beta;
  real<lower=0> sigma;  // error SD
  vector[M] phi_raw;
  vector[T] epsilon;
  real<lower=-1, upper=1> rho;  // autocorrelation coefficient
  real<lower=0> tau;
}

transformed parameters {
  vector[T] theta;  // AR1-structured error terms
  vector[M] phi = phi_raw * tau;  // tree random effects
  theta[1] = epsilon[1];
  for (t in 2:T) {
    theta[t] = rho * theta[t - 1] + epsilon[t];
  }
}

model {
  vector[N] mu = log_ks_ref + x * beta + phi[tree] + theta[time];
  beta ~ normal(0, 2.5);
  epsilon ~ normal(0, 2.5);
  tau ~ normal(0, 2.5);
  phi_raw ~ normal(0, 1);
  log_ks ~ normal(mu, sigma);
  // log_ks ~ normal(log_ks_ref + x * beta + phi[tree] + theta[time], sigma);
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) log_lik[n] = normal_lpdf(log_ks[n] |
    log_ks_ref + x[n] * beta + phi[tree[n]] + theta[time[n]], sigma);
}
