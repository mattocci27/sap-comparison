# ────────────────────────────────────────────────────────────────────────────────
compute_species_all_r2 <- function(draws_df, stan_data) {
  beta_int_col <- "beta_1"
  beta_slope_col <- "beta_2"

  draws_cleaned <- draws_df %>%
    janitor::clean_names()

  xk <- stan_data$xk

  beta_int_pred <- draws_cleaned %>%
    select(matches(beta_int_col)) |>
    as.matrix()

  beta_slope_pred <- draws_cleaned %>%
    select(matches(beta_slope_col)) |>
    as.matrix()

  log_int_pred <- xk %*% t(beta_int_pred)
  log_slope_pred <- xk %*% t(beta_slope_pred)

  log_int_obs <- draws_cleaned %>%
    select(starts_with("beta_hat_1")) %>%
    as.matrix() %>%
    t()

  log_slope_obs <- draws_cleaned %>%
    select(starts_with("beta_hat_2")) %>%
    as.matrix() %>%
    t()

  dim(log_int_pred)
  dim(log_int_obs)

  r2_int <- apply(log_int_pred, 2, var) / (apply(log_int_pred, 2, var) + apply(log_int_obs - log_int_pred, 2, var))
  r2_int_q <- quantile(r2_int, c(0.025, 0.5, 0.975)) |> round(3)

  r2_slope <- apply(log_slope_pred, 2, var) / (apply(log_slope_pred, 2, var) + apply(log_slope_obs - log_slope_pred, 2, var))
  r2_slope_q <- quantile(r2_slope, c(0.025, 0.5, 0.975)) |> round(3)

  list(
    r2_co_a = r2_int_q,
    r2_co_b = r2_slope_q
  )
}


compute_species_single_r2 <- function(draws_df, stan_data) {

  draws_cleaned <- draws_df %>%
    janitor::clean_names()

  xk <- stan_data$xk

  beta_int_pred <- draws_cleaned %>%
    select(beta_1_1, beta_2_1) |>
    as.matrix()

  beta_slope_pred <- draws_cleaned %>%
    select(beta_1_2, beta_2_2) |>
    as.matrix()

  log_int_pred <- xk %*% t(beta_int_pred)
  log_slope_pred <- xk %*% t(beta_slope_pred)

  log_int_obs <- draws_cleaned %>%
    select(starts_with("beta_hat_1")) %>%
    as.matrix() %>%
    t()

  log_slope_obs <- draws_cleaned %>%
    select(starts_with("beta_hat_2")) %>%
    as.matrix() %>%
    t()

  dim(log_int_pred)
  dim(log_int_obs)

  r2_int <- apply(log_int_pred, 2, var) / (apply(log_int_pred, 2, var) + apply(log_int_obs - log_int_pred, 2, var))
  r2_int_q <- quantile(r2_int, c(0.025, 0.5, 0.975)) |> round(3)

  r2_slope <- apply(log_slope_pred, 2, var) / (apply(log_slope_pred, 2, var) + apply(log_slope_obs - log_slope_pred, 2, var))
  r2_slope_q <- quantile(r2_slope, c(0.025, 0.5, 0.975)) |> round(3)

  list(
    r2_co_a = r2_int_q,
    r2_co_b = r2_slope_q
  )
}


compute_segments_single_r2 <- function(draws_df, stan_data) {

  draws_cleaned <- draws_df %>%
    janitor::clean_names()

  xj <- stan_data$xj

  beta_int_pred <- draws_cleaned %>%
    select(beta_1, beta_2) |>
    as.matrix()

  beta_slope_pred <- draws_cleaned %>%
    select(beta_3, beta_4) |>
    as.matrix()

  log_int_pred <- xj %*% t(beta_int_pred)
  log_slope_pred <- xj %*% t(beta_slope_pred)

  log_int_obs <- draws_cleaned %>%
    select(starts_with("a_1")) %>%
    as.matrix() %>%
    t()

  log_slope_obs <- draws_cleaned %>%
    select(starts_with("a_2")) %>%
    as.matrix() %>%
    t()

  dim(log_int_pred)
  dim(log_int_obs)

  r2_int <- apply(log_int_pred, 2, var) / (apply(log_int_pred, 2, var) + apply(log_int_obs - log_int_pred, 2, var))
  r2_int_q <- quantile(r2_int, c(0.025, 0.5, 0.975)) |> round(3)

  r2_slope <- apply(log_slope_pred, 2, var) / (apply(log_slope_pred, 2, var) + apply(log_slope_obs - log_slope_pred, 2, var))
  r2_slope_q <- quantile(r2_slope, c(0.025, 0.5, 0.975)) |> round(3)

  list(
    r2_co_a = r2_int_q,
    r2_co_b = r2_slope_q
  )
}

compute_segments_all_r2 <- function(draws_df, stan_data) {

  draws_cleaned <- draws_df %>%
    janitor::clean_names()

  alpha_species_tree_1 <- draws_cleaned %>%
    select(matches("alpha_species_tree_1")) |>
    as.matrix()

  alpha_trait_tree_1 <- draws_cleaned %>%
    select(matches("alpha_trait_tree_1")) |>
    as.matrix()

  log_int_pred <- alpha_species_tree_1 + alpha_trait_tree_1

  alpha_species_tree_2 <- draws_cleaned %>%
    select(matches("alpha_species_tree_2")) |>
    as.matrix()

  alpha_trait_tree_2 <- draws_cleaned %>%
    select(matches("alpha_trait_tree_2")) |>
    as.matrix()

  log_slope_pred <- alpha_species_tree_2 + alpha_trait_tree_2

  log_int_obs <- draws_cleaned %>%
    select(starts_with("a_1")) %>%
    as.matrix() %>%
    t()

  log_slope_obs <- draws_cleaned %>%
    select(starts_with("a_2")) %>%
    as.matrix() %>%
    t()

  log_slope_pred <- log_slope_pred |> t()
  log_int_pred <- log_int_pred |> t()

  dim(log_int_obs) |> print()
  dim(log_int_pred) |> print()


  r2_int <- apply(log_int_pred, 2, var) / (apply(log_int_pred, 2, var) + apply(log_int_obs - log_int_pred, 2, var))
  r2_int_q <- quantile(r2_int, c(0.025, 0.5, 0.975)) |> round(3)

  r2_slope <- apply(log_slope_pred, 2, var) / (apply(log_slope_pred, 2, var) + apply(log_slope_obs - log_slope_pred, 2, var))
  r2_slope_q <- quantile(r2_slope, c(0.025, 0.5, 0.975)) |> round(3)

  list(
    r2_co_a = r2_int_q,
    r2_co_b = r2_slope_q
  )
}
