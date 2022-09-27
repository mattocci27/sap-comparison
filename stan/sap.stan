data {
  int<lower=0> Ni; // number of data
  int<lower=0> Nj; // number of tree samples
  int<lower=0> Nk; // number of species
  int<lower=0> Nl; // number of xylem types
  int<lower=0> Mi; // number of predcitors for data level
  int<lower=0> Mj; // number of predcitors for tree level
  int<lower=0> Mk; // number of predcitors for species level
  int<lower=0> Ml; // number of predcitors for xylem level
  array[Ni] int<lower=1,upper=Nj> jj; // tree sample
  array[Nj] int<lower=1,upper=Nk> jj2; // sp lab from tree
  array[Ni] int<lower=1,upper=Nk> kk; // species
  array[Nk] int<lower=1,upper=Nl> kk2; // xylem lab from sp
  array[Ni] int<lower=1,upper=Nl> ll; // xylem type
  vector[Ni] y;
  vector[Ni] p_g;
  vector[Ni] rho;
  vector[Ni] x;
  // matrix[Mi,Ni] u; // data-level predictor
  // matrix[Ni,Mi] xi1; // data-level predcitor (p_g)
  // matrix[Mi,Mj] xj; // sample-level predictor (tarits)
  // matrix[1,Mj] xk; // sp-level intercept
  // matrix[1,Mj] xl; // xylem-level intercept
}

parameters {
  real<lower=0> sigma;
  vector[Mj] mu_alpha_l_tilde;
  vector<lower=0,upper=pi()/2>[Mj] tau_unif_alpha_l;
  vector<lower=0,upper=pi()/2>[Mj] tau_unif_alpha_k;
  vector<lower=0,upper=pi()/2>[Mj-1] tau_unif_alpha_j;

  vector[Mj] mu_beta_l_tilde;
  vector<lower=0,upper=pi()/2>[Mj] tau_unif_beta_l;
  vector<lower=0,upper=pi()/2>[Mj] tau_unif_beta_k;
  vector<lower=0,upper=pi()/2>[Mj-1] tau_unif_beta_j;

  array[Mj] vector[Nl] theta_alpha_l;
  array[Mj] vector[Nk] theta_alpha_k;
  array[Mj-1] vector[Nj] theta_alpha_j;
  array[Mj] vector[Nl] theta_beta_l;
  array[Mj] vector[Nk] theta_beta_k;
  array[Mj-1] vector[Nj] theta_beta_j;

  cholesky_factor_corr[Mi] L_Omega;
  vector<lower=0,upper=pi()/2>[Mi] tau_unif;
  matrix[Mi, Ni] z;
}

transformed parameters{
  matrix[Ni,Mi] gamma;
  vector[Ni] log_a;
  vector[Ni] b;
  vector<lower=0>[Mi] tau;

  array[Mj] vector[Nl] alpha_l;
  array[Mj] vector[Nk] alpha_k;
  array[Mj-1] vector[Nj] alpha_j;
  vector<lower=0>[Mj] tau_alpha_l;
  vector<lower=0>[Mj] tau_alpha_k;
  vector<lower=0>[Mj-1] tau_alpha_j;
  vector[Ni] log_mu_a;

  array[Mj] vector[Nl] beta_l;
  array[Mj] vector[Nk] beta_k;
  array[Mj-1] vector[Nj] beta_j;
  vector<lower=0>[Mj] tau_beta_l;
  vector<lower=0>[Mj] tau_beta_k;
  vector<lower=0>[Mj-1] tau_beta_j;
  vector[Ni] mu_b;

  for (i in 1:Mi) tau[i] = 2.5 * tan(tau_unif[i]);

  for (j in 1:Mj) {
    tau_alpha_l[j] = 2.5 * tan(tau_unif_alpha_l[j]);
    tau_alpha_k[j] = 2.5 * tan(tau_unif_alpha_k[j]);
    alpha_l[j] = mu_alpha_l_tilde[j] + tau_alpha_l[j] * theta_alpha_l[j];
    for (k in 1:Nk) {
      alpha_k[j][k] = alpha_l[j][kk2[k]] +
        tau_alpha_k[j] * theta_alpha_k[j][k];
        }
    }
  for (j in 1:(Mj-1)) {
    tau_alpha_j[j] = 2.5 * tan(tau_unif_alpha_j[j]);
    for (nj in 1:Nj) {
      alpha_j[j][nj] = alpha_k[j][jj2[nj]] +
       tau_alpha_j[j] * theta_alpha_j[j][nj];
     }
  }
  for (j in 1:Mj) {
    tau_beta_l[j] = 2.5 * tan(tau_unif_beta_l[j]);
    tau_beta_k[j] = 2.5 * tan(tau_unif_beta_k[j]);
    beta_l[j] = mu_beta_l_tilde[j] + tau_beta_l[j] * theta_beta_l[j];
    for (k in 1:Nk) {
      beta_k[j][k] = beta_l[j][kk2[k]] +
        tau_beta_k[j] * theta_beta_k[j][k];
    }
  }
  for (j in 1:(Mj-1)) {
    tau_beta_j[j] = 2.5 * tan(tau_unif_beta_j[j]);
    for (nj in 1:Nj) {
      beta_j[j][nj] = beta_k[j][jj2[nj]] +
       tau_beta_j[j] * theta_beta_j[j][nj];
     }
  }
  for (n in 1:Ni) {
    log_mu_a[n] = alpha_j[1][jj[n]] +
    (alpha_j[2][jj[n]]) * p_g[n] +
    (alpha_k[3][kk[n]]) * rho[n];

    mu_b[n] = beta_j[1][jj[n]] +
    (beta_j[2][jj[n]]) * p_g[n] +
    (beta_k[3][kk[n]]) * rho[n];
  }
  gamma = append_col(log_mu_a, mu_b) + (diag_pre_multiply(tau, L_Omega) * z)';
  log_a = gamma[,1];
  b = gamma[,2];
}
model {
  mu_alpha_l_tilde ~ normal(0, 5);
  mu_beta_l_tilde ~ normal(0, 5);
  for (i in 1:Mj) {
    theta_alpha_l[i] ~ std_normal();
    theta_beta_l[i] ~ std_normal();
    theta_alpha_k[i] ~ std_normal();
    theta_beta_k[i] ~ std_normal();
  }
  for (i in 1:2) {
    theta_alpha_j[i] ~ std_normal();
    theta_beta_j[i] ~ std_normal();
  }

  for (n in 1:Ni) {
    y[n]  ~ normal(log_a[n] + b[n] * x, sigma);
  }
}



