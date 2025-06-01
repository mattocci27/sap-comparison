# ────────────────────────────────────────────────────────────────────────────────
# Function 1: Compute species‐level R² (based on xk)
# ────────────────────────────────────────────────────────────────────────────────
compute_species_r2 <- function(draws_df, stan_data) {
  # draws_df  : posterior‐draws tibble (from fit_all_draws_segments_noxylem_traits_sp_simple_all)
  # stan_data : list containing xk (K×T), uj (K×J), plus J, K, T

  xk <- stan_data$xk       # K × T matrix of species‐level predictors
  uj <- stan_data$uj       # K × J mapping from species → tree
  K  <- nrow(xk)
  T  <- ncol(xk)
  J  <- stan_data$J
  M  <- nrow(draws_df)     # number of posterior draws

  # 1) Extract beta[1,·] and beta[2,·] into (M × T) matrices:
  beta1_mat <- draws_df %>%
    select(matches("^beta\\[1,[0-9]+\\]$")) %>%
    as.matrix()
  stopifnot(ncol(beta1_mat) == T)

  beta2_mat <- draws_df %>%
    select(matches("^beta\\[2,[0-9]+\\]$")) %>%
    as.matrix()
  stopifnot(ncol(beta2_mat) == T)

  # 2) Extract tau_k[1], tau_k[2] into (M × 2):
  tau_k_mat <- draws_df %>%
    select(matches("^tau_k\\[[12]\\]$")) %>%
    as.matrix()
  stopifnot(ncol(tau_k_mat) == 2)

  # 3) Extract L_Omega_k into an array (M × 2 × 2):
  L11_k <- draws_df %>% pull(`L_Omega_k[1,1]`)
  L21_k <- draws_df %>% pull(`L_Omega_k[2,1]`)
  L22_k <- draws_df %>% pull(`L_Omega_k[2,2]`)

  Lk_array <- array(0.0, dim = c(M, 2, 2))
  Lk_array[,1,1] <- L11_k
  Lk_array[,1,2] <- 0
  Lk_array[,2,1] <- L21_k
  Lk_array[,2,2] <- L22_k

  # 4) Extract zk[1,1..K], zk[2,1..K] into (M × 2 × K):
  zk1_mat <- draws_df %>%
    select(matches("^zk\\[1,[0-9]+\\]$")) %>%
    as.matrix()  # M × K
  zk2_mat <- draws_df %>%
    select(matches("^zk\\[2,[0-9]+\\]$")) %>%
    as.matrix()  # M × K
  stopifnot(ncol(zk1_mat) == K, ncol(zk2_mat) == K)

  zk_array <- array(0.0, dim = c(M, 2, K))
  zk_array[,1, ] <- zk1_mat
  zk_array[,2, ] <- zk2_mat

  # 5) Build species_noise_array (M × 2 × K):
  species_noise_array <- array(0.0, dim = c(M, 2, K))
  for (m in seq_len(M)) {
    Dm  <- diag(tau_k_mat[m, ], 2, 2)          # diag(tau_k[m,1], tau_k[m,2])
    LkD <- Lk_array[m, , ] %*% Dm              # (2×2) × (2×2) → 2×2
    species_noise_array[m, , ] <- LkD %*% zk_array[m, , ]  # (2×2) × (2×K) → 2×K
  }

  # 6) Extract tau_j[1], tau_j[2] into (M × 2):
  tau_j_mat <- draws_df %>%
    select(matches("^tau_j\\[[12]\\]$")) %>%
    as.matrix()
  stopifnot(ncol(tau_j_mat) == 2)

  # 7) Extract L_Omega_j into array (M × 2 × 2):
  L11_j <- draws_df %>% pull(`L_Omega_j[1,1]`)
  L21_j <- draws_df %>% pull(`L_Omega_j[2,1]`)
  L22_j <- draws_df %>% pull(`L_Omega_j[2,2]`)

  Lj_array <- array(0.0, dim = c(M, 2, 2))
  Lj_array[,1,1] <- L11_j
  Lj_array[,1,2] <- 0
  Lj_array[,2,1] <- L21_j
  Lj_array[,2,2] <- L22_j

  # 8) Extract zj[1,1..J], zj[2,1..J] into (M × 2 × J):
  zj1_mat <- draws_df %>%
    select(matches("^zj\\[1,[0-9]+\\]$")) %>%
    as.matrix()  # M × J
  zj2_mat <- draws_df %>%
    select(matches("^zj\\[2,[0-9]+\\]$")) %>%
    as.matrix()  # M × J
  stopifnot(ncol(zj1_mat) == J, ncol(zj2_mat) == J)

  zj_array <- array(0.0, dim = c(M, 2, J))
  zj_array[,1, ] <- zj1_mat
  zj_array[,2, ] <- zj2_mat

  # 9) Loop over draws to build A_fixed, A_total and compute R²_intercept and R²_slope:
  A_fixed_array  <- array(0.0, dim = c(M, 2, J))
  A_total_array  <- array(0.0, dim = c(M, 2, J))
  R2_int_draws   <- numeric(M)
  R2_slo_draws   <- numeric(M)

  for (m in seq_len(M)) {
    # (a) Build beta_hat_fixed_species[m,, ] = β[m,, ] × t(xk):
    beta_m      <- rbind(beta1_mat[m, ], beta2_mat[m, ])  # 2×T
    beta_hat_f  <- beta_m %*% t(xk)                       # 2×K
    # (b) Add species noise:
    beta_hat_t  <- beta_hat_f + species_noise_array[m, , ]   # 2×K
    # (c) A_fixed[m,, ] = beta_hat_f (2×K) × uj (K×J) = 2×J
    A_fixed_array[m, , ] <- beta_hat_f %*% uj
    # (d) A_total[m,, ] = beta_hat_t (2×K) × uj + tree_noise (2×J)
    A_sp_total <- beta_hat_t %*% uj                            # 2×J
    tree_noise <- (Lj_array[m, , ] %*% diag(tau_j_mat[m, ], 2, 2)) %*% zj_array[m, , ]
    A_total_array[m, , ] <- A_sp_total + tree_noise            # 2×J

    # (e) R²_intercept[m]:
    intercept_fixed <- A_fixed_array[m, 1, ]  # length J
    intercept_total <- A_total_array[m,  1, ]
    vfit0  <- var(intercept_fixed)
    vtot0  <- var(intercept_total)
    R2_int_draws[m] <- if (vtot0 > 0) vfit0 / vtot0 else NA_real_

    # (f) R²_slope[m]:
    slope_fixed <- A_fixed_array[m, 2, ]
    slope_total <- A_total_array[m,  2, ]
    vfit1  <- var(slope_fixed)
    vtot1  <- var(slope_total)
    R2_slo_draws[m] <- if (vtot1 > 0) vfit1 / vtot1 else NA_real_
  }

  # 10) Summarize:
  result <- list(
    r2_intercept_draws = R2_int_draws,
    r2_slope_draws     = R2_slo_draws,
    r2_intercept = tibble(
      mean   = mean(R2_int_draws,   na.rm = TRUE),
      median = median(R2_int_draws, na.rm = TRUE),
      ci2_5  = quantile(R2_int_draws, 0.025, na.rm = TRUE),
      ci97_5 = quantile(R2_int_draws, 0.975, na.rm = TRUE)
    ),
    r2_slope = tibble(
      mean   = mean(R2_slo_draws,   na.rm = TRUE),
      median = median(R2_slo_draws, na.rm = TRUE),
      ci2_5  = quantile(R2_slo_draws, 0.025, na.rm = TRUE),
      ci97_5 = quantile(R2_slo_draws, 0.975, na.rm = TRUE)
    )
  )

  return(result)
}

# ────────────────────────────────────────────────────────────────────────────────
# Function 2: Compute segment‐level R² (based on xj)
# ────────────────────────────────────────────────────────────────────────────────
compute_segment_r2 <- function(draws_df, stan_data) {
  # draws_df  : posterior‐draws tibble (from fit_all_draws_segments_noxylem_traits_simple_all)
  # stan_data : list containing xj (J×T), uj (K×J), plus J, K, T

  xj <- stan_data$xj    # J × T matrix of tree‐level predictors
  uj <- stan_data$uj    # K × J mapping species → tree
  J  <- stan_data$J
  K  <- stan_data$K
  T  <- stan_data$T     # number of tree‐level traits
  M  <- nrow(draws_df)

  # 1) Find the correct column names for B[1,*] and B[2,*]:
  all_names <- colnames(draws_df)
  # Try loose pattern with optional space:
  b1_cols <- grep("^B\\[1, *[0-9]+\\]$", all_names, value = TRUE)
  b2_cols <- grep("^B\\[2, *[0-9]+\\]$", all_names, value = TRUE)
  if (length(b1_cols) != T || length(b2_cols) != T) {
    stop("Could not find exactly T columns for B[1,*] or B[2,*].\n",
         "Please inspect colnames(draws_df) for the correct naming pattern.")
  }

  # 2) Build B1_mat, B2_mat (M × T):
  B1_mat <- draws_df %>% select(all_of(b1_cols)) %>% as.matrix()
  stopifnot(ncol(B1_mat) == T)
  B2_mat <- draws_df %>% select(all_of(b2_cols)) %>% as.matrix()
  stopifnot(ncol(B2_mat) == T)

  # 3) Reconstruct species_noise_array (M × 2 × K) exactly as in Stan:
  tau_k_mat <- draws_df %>% select(matches("^tau_k\\[[12]\\]$")) %>% as.matrix()
  stopifnot(ncol(tau_k_mat) == 2)
  L11_k <- draws_df %>% pull(`L_Omega_k[1,1]`)
  L21_k <- draws_df %>% pull(`L_Omega_k[2,1]`)
  L22_k <- draws_df %>% pull(`L_Omega_k[2,2]`)
  Lk_array <- array(0.0, dim = c(M, 2, 2))
  Lk_array[,1,1] <- L11_k;  Lk_array[,1,2] <- 0
  Lk_array[,2,1] <- L21_k;  Lk_array[,2,2] <- L22_k
  zk1_mat <- draws_df %>% select(matches("^zk\\[1,[0-9]+\\]$")) %>% as.matrix()
  zk2_mat <- draws_df %>% select(matches("^zk\\[2,[0-9]+\\]$")) %>% as.matrix()
  stopifnot(ncol(zk1_mat) == K, ncol(zk2_mat) == K)
  zk_array <- array(0.0, dim = c(M, 2, K))
  zk_array[,1, ] <- zk1_mat
  zk_array[,2, ] <- zk2_mat
  species_noise_array <- array(0.0, dim = c(M, 2, K))
  for (m in seq_len(M)) {
    Dk   <- diag(tau_k_mat[m, ], 2, 2)
    LkD  <- Lk_array[m, , ] %*% Dk
    species_noise_array[m, , ] <- LkD %*% zk_array[m, , ]
  }

  # 4) Reconstruct tree_noise_array (M × 2 × J):
  tau_j_mat <- draws_df %>% select(matches("^tau_j\\[[12]\\]$")) %>% as.matrix()
  stopifnot(ncol(tau_j_mat) == 2)
  L11_j <- draws_df %>% pull(`L_Omega_j[1,1]`)
  L21_j <- draws_df %>% pull(`L_Omega_j[2,1]`)
  L22_j <- draws_df %>% pull(`L_Omega_j[2,2]`)
  Lj_array <- array(0.0, dim = c(M, 2, 2))
  Lj_array[,1,1] <- L11_j;  Lj_array[,1,2] <- 0
  Lj_array[,2,1] <- L21_j;  Lj_array[,2,2] <- L22_j
  zj1_mat <- draws_df %>% select(matches("^zj\\[1,[0-9]+\\]$")) %>% as.matrix()
  zj2_mat <- draws_df %>% select(matches("^zj\\[2,[0-9]+\\]$")) %>% as.matrix()
  stopifnot(ncol(zj1_mat) == J, ncol(zj2_mat) == J)
  zj_array <- array(0.0, dim = c(M, 2, J))
  zj_array[,1, ] <- zj1_mat
  zj_array[,2, ] <- zj2_mat
  tree_noise_array <- array(0.0, dim = c(M, 2, J))
  for (m in seq_len(M)) {
    Dj   <- diag(tau_j_mat[m, ], 2, 2)
    LjD  <- Lj_array[m, , ] %*% Dj
    tree_noise_array[m, , ] <- LjD %*% zj_array[m, , ]
  }

  # 5) Loop over draws to build A and compute R²_intercept & R²_slope:
  R2_int_draws <- numeric(M)
  R2_slo_draws <- numeric(M)

  for (m in seq_len(M)) {
    # (a) θ_centered (2×K) = species_noise_array[m, , ]
    theta_c <- species_noise_array[m, , ]
    # (b) α_species_tree (2×J) = theta_c %*% uj
    alpha_species_tree <- theta_c %*% uj
    # (c) α_trait_tree (2×J) = B[m,, ] %*% t(xj)
    Bm               <- rbind(B1_mat[m, ], B2_mat[m, ])  # 2×T
    alpha_trait_tree <- Bm %*% t(xj)                     # 2×J
    # (d) tree_noise (2×J) = tree_noise_array[m, , ]
    tree_noise       <- tree_noise_array[m, , ]
    # (e) A = species + trait + noise:
    A_m <- alpha_species_tree + alpha_trait_tree + tree_noise  # 2×J

    # (f) R²_intercept:
    vfit0      <- var(alpha_trait_tree[1, ])
    vtotal_int <- var(A_m[1, ])
    R2_int_draws[m] <- if (vtotal_int > 0) vfit0 / vtotal_int else NA_real_

    # (g) R²_slope:
    vfit1      <- var(alpha_trait_tree[2, ])
    vtotal_sl  <- var(A_m[2, ])
    R2_slo_draws[m] <- if (vtotal_sl > 0) vfit1 / vtotal_sl else NA_real_
  }

  # 6) Summarize:
  result <- list(
    r2_intercept_draws = R2_int_draws,
    r2_slope_draws     = R2_slo_draws,
    r2_intercept = tibble(
      mean   = mean(R2_int_draws,   na.rm = TRUE),
      median = median(R2_int_draws, na.rm = TRUE),
      ci2_5  = quantile(R2_int_draws, 0.025, na.rm = TRUE),
      ci97_5 = quantile(R2_int_draws, 0.975, na.rm = TRUE)
    ),
    r2_slope = tibble(
      mean   = mean(R2_slo_draws,   na.rm = TRUE),
      median = median(R2_slo_draws, na.rm = TRUE),
      ci2_5  = quantile(R2_slo_draws, 0.025, na.rm = TRUE),
      ci97_5 = quantile(R2_slo_draws, 0.975, na.rm = TRUE)
    )
  )

  return(result)
}

# ────────────────────────────────────────────────────────────────────────────────
# Example usage:
#    species_res <- compute_species_R2(fit_all_draws_segments_noxylem_traits_sp_simple_all, stan_data_noxylem_all)
#    print(species_res$R2_intercept)
#    print(species_res$R2_slope)
#
#    segment_res <- compute_segment_R2(fit_all_draws_segments_noxylem_traits_simple_all, stan_data_noxylem_all)
#    print(segment_res$R2_intercept)
#    print(segment_res$R2_slope)
# ────────────────────────────────────────────────────────────────────────────────
