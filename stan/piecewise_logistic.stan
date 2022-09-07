data {
  int<lower=0> N; // number of samples
  int<lower=0> J; // number of species
  int<lower=0> K; // number of parameters (int, beta1, beta2, bp)
  array[N] int<lower=1,upper=J> jj; // species
  array[N] int<lower=0> y; // count
  array[N] int<lower=0> total; // total
  vector[N] x; // presure
  matrix[1, J] u; // sp-level intercept
}

parameters {
  real gamma1;
  real gamma2;
  real gamma3;
  real<lower=0,upper=0.08> gamma4;
  matrix[K, J] z;
  cholesky_factor_corr[K] L_Omega;
  vector<lower=0,upper=pi()/2>[K] tau_unif;
}

transformed parameters{
  matrix[K, J] beta;
  matrix[K, 1] gamma;
  vector<lower=0>[K] tau;
  vector[N] alpha;
  gamma[1, 1] = gamma1;
  gamma[2, 1] = gamma2;
  gamma[3, 1] = gamma3;
  gamma[4, 1] = gamma4;

  for (k in 1:K) tau[k] = 2.5 * tan(tau_unif[k]);
  beta = gamma * u + diag_pre_multiply(tau, L_Omega) * z;

  for (n in 1:N) {
    if (x[n] < beta[4, jj[n]]) {
      alpha[n] = beta[1, jj[n]] + beta[2, jj[n]] * x[n];
    } else {
      alpha[n] = beta[1, jj[n]] + beta[2, jj[n]] * x[n] +
         beta[3, jj[n]] * (x[n] - beta[4, jj[n]]);
    }
  }
}

model {
  to_vector(z) ~ std_normal();
  L_Omega ~ lkj_corr_cholesky(2);
  to_vector(gamma) ~ normal(0, 5);
  y ~ binomial_logit(total, alpha);
}
