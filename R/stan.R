
generate_anova_data <- function(data, log = FALSE, err = FALSE) {
  data <- read_csv(data)
  list_data <- list(
    N = nrow(data),
    jj = data$species |>  as.factor() |> as.numeric(),
    kk = data$pressure |>  as.factor() |> as.numeric(),
    mm = data$tree_id |>  as.factor() |> as.numeric(),
    mk = paste(data$tree_id, data$pressure, sep = "_")
      |> as.factor() |> as.numeric(),
    jk = paste(data$species, data$pressure, sep = "_")
      |> as.factor() |> as.numeric(),
    y = data$pres_calib_mean - data$tens_calib_mean,
    y1 = data$pres_calib_mean,
    y2 = data$tens_calib_mean
  )

  if (log) {
    list_data$y <- log(data$pres_calib_mean / data$tens_calib_mean)
    list_data$y1 <- log(data$pres_calib_mean)
    list_data$y2 <- log(data$tens_calib_mean)
  }
  if (err) {
    list_data$sig1 <- data$pres_calib_sd
    list_data$sig2 <- data$tens_calib_sd
  }

  list_data$J <- unique(data$species) |> length()
  list_data$K <- unique(data$pressure) |> length()
  list_data$M <- unique(data$tree_id) |> length()
  list_data$JK <- unique(list_data$jk) |> length()
  list_data$MK <- unique(list_data$mk) |> length()
  list_data
}

generate_dummy_data <- function(n = 30, n_alpha = 5, n_beta = 5, sigma_alpha = 0.9, sigma_beta = 0.9, sigma_gamma = 0.5, mu_hat = 1.3, sigma = 2, seed = 123) {
  set.seed(seed)
  # set.seed(123)
  # sigma_alpha <- 3
  # sigma_beta <- 3
  # sigma_gamma <- 0.2
  # sigma <- 1
  # mu_hat <- 1
  # n <- 3
  # n_alpha <- 5
  # n_beta <- 3
  n_gamma <- n_alpha * n_beta
  alpha <- rnorm(n_alpha, 0, sigma_alpha)
  beta <- rnorm(n_beta, 0, sigma_beta)
  gamma <- rnorm(n_alpha * n_beta, 0, sigma_gamma)

  effect_raw <- rep(alpha, each = n_beta) + rep(beta, n_alpha) + gamma
  # effect <- mapply(rnorm, n, effect_raw, sigma) |> as.numeric() |> round(2)
  effect <- mapply(rnorm, n, effect_raw, sigma) |> as.numeric()
  mu <- mu_hat + effect

  list_data <- list()
  list_data$y <- mu
  list_data$J <- n_alpha
  list_data$K <- n_beta
  list_data$JK <- n_alpha * n_beta
  list_data$N <- length(mu)
  list_data$jj <- rep(1:n_alpha, each = n_beta) |> rep(each = n)
  list_data$kk <- rep(1:n_beta, n_alpha) |> rep(each = n)
  list_data$jk <- rep(1:n_gamma, each = n)
  list_data$alpha <- alpha
  list_data$beta <- beta
  list_data$gamma <- gamma
  list_data$sigma_alpha <- sigma_alpha
  list_data$sigma_beta <- sigma_beta
  list_data$sigma_gamma <- sigma_gamma
  list_data$sigma <- sigma

  alpha_ <- rep(alpha, each = n_beta) |> rep(each = n)
  beta_ <- rep(beta, n_alpha) |> rep(each = n)
  gamma_ <- rep(gamma, each = n)
  # alpha_ <- rep(alpha, each = 3) |> rep(each = n) |> round(2)
  # beta_ <- rep(beta, 5) |> rep(each = n) |> round(2)
  # gamma_ <- rep(gamma, each = n) |> round(2)

  # tibble(alpha_, beta_, gamma_, effect, jj = list_data$jj, kk = list_data$kk, jk = list_data$jk) |>
  #   DT::datatable()

  list_data$data <- tibble(
    y = mu,
    species = list_data$jj,
    pressure = list_data$kk,
    inter = list_data$jk,
    effect = alpha_ + beta_ + gamma_
  )
  list_data
}

#' @title Create summary table for posteriors
create_stan_tab <- function(draws) {
  tmp <- draws |>
    janitor::clean_names() |>
    dplyr::select(contains(c("alpha", "beta", "gamma","delta", "effect", "pred", "mu_hat")))
  mean_ <- apply(tmp, 2, mean)
  q2_5 <- apply(tmp, 2, \(x)(quantile(x, 0.025)))
  q5 <- apply(tmp, 2, \(x)(quantile(x, 0.05)))
  q97_5 <- apply(tmp, 2, \(x)(quantile(x, 0.975)))
  q95 <- apply(tmp, 2, \(x)(quantile(x, 0.9)))
  tibble(para = names(mean_), mean_, q2_5, q5, q95, q97_5)
}

# targets::tar_load(ks_trees)
# data <- read_csv(ks_trees)

generate_anova_mvn_data <- function(data, model_type = c("normal", "log-normal", "gamma")) {
  data <- read_csv(data) |>
    filter(!is.na(ks))
  list_data <- list(
    N = nrow(data),
    jj = data$species |>  as.factor() |> as.numeric(),
    kk = data$pressure |>  as.factor() |> as.numeric(),
    jk = paste(data$species, data$pressure, sep = "_")
      |> as.factor() |> as.numeric(),
    ll = data$pres_type |>  as.factor() |> as.numeric(),
    mm = data$tree_id |>  as.factor() |> as.numeric(),
    y = data$ks
  )
    list_data$J <- unique(data$species) |> length()
    list_data$K <- unique(data$pressure) |> length()
    list_data$JK <- unique(list_data$jk) |> length()
    list_data$L <- unique(list_data$ll) |> length()
    list_data$M <- unique(list_data$mm) |> length()

    list_data$model_type <- case_when(
      model_type == "normal" ~ 1,
      model_type == "log-normal" ~ 2,
      model_type == "gamma" ~ 3,
    )

  list_data
}
# hoge <- data |>
#   group_by(species, tree_id, pres_type, pressure) |>
#   summarize(mean_ = mean(ks), var_ = var(ks))

# plot(mean_ ~ var_ , hoge, log = "xy")

#tar_read(anova_mvn_data)


write_anova_yml <- function(output, draws, ll = 0.025, hh = 0.975) {
  tmp <- draws |>
      janitor::clean_names()  |>
      dplyr::select(tau_1, tau_2, sigma) |>
      mutate(var_ = tau_1^2 + tau_2^2 + sigma^2) |>
      mutate(var_tau_1 = tau_1^2 / var_) |>
      mutate(var_tau_2 = tau_2^2 / var_) |>
      mutate(var_sigma = sigma^2 / var_)

  median_ <- apply(tmp * 100, 2, median) |> round(1)
  lwr <- apply(tmp, 2, \(x)quantile(x * 100, ll)) |> round(1)
  upr <- apply(tmp, 2, \(x)quantile(x * 100, hh)) |> round(1)

  out <- file(paste(output), "w") # write
  writeLines(
    paste0("species: ",
      median_["var_tau_1"], "% [",
      lwr["var_tau_1"],
      ", ",
      upr["var_tau_1"],
      "]"
    ),
    out,
    sep = "\n")
  writeLines(
    paste0("pressure: ",
      median_["var_tau_2"], "% [",
      lwr["var_tau_2"],
      ", ",
      upr["var_tau_2"],
      "]"
    ),
    out,
    sep = "\n")
  writeLines(
    paste0("residuals: ",
      median_["var_sigma"], "% [",
      lwr["var_sigma"],
      ", ",
      upr["var_sigma"],
      "]"
    ),
    out,
    sep = "\n")

   close(out)
   # The return value must be a vector of paths to the files we write:
   paste(output)
}


generate_logistic_stan_data <- function(data) {
  # library(tidyverse)
  # library(targets)
  tar_load(cond_count_csv)
  data <- cond_count_csv

  d <- read_csv(data) |>
    filter(!is.na(count))

  list(
    N = nrow(d),
    J = unique(d$species) |> length(),
    K = 4,
    jj = as.factor(d$species) |> as.numeric(),
    y = d$count,
    total = d$total,
    x = d$pressure,
    u = t(as.matrix(rep(1, 5)))
  )

}
