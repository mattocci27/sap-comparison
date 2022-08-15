
generate_anova_data <- function(data, log = FALSE) {
  data <- read_csv(data)
  list_data <- list(
    N = nrow(data),
    jj = data$species |>  as.factor() |> as.numeric(),
    kk = data$pressure |>  as.factor() |> as.numeric(),
    jk = paste(data$species, data$pressure, sep = "_")
      |> as.factor() |> as.numeric(),
    y = data$pres_calib - data$tens_calib,
    y1 = data$pres_calib,
    y2 = data$tens_calib
  )

  if (log) {
    list_data$y <- log(data$pres_calib / data$tens_calib)
    list_data$y1 <- log(data$pres_calib)
    list_data$y2 <- log(data$tens_calib)
  }

  list_data$J <- unique(data$species) |> length()
  list_data$K <- unique(data$pressure) |> length()
  list_data$JK <- unique(list_data$jk) |> length()
  list_data
}

generate_dummy_data <- function(n = 30, sigma_alpha = 0.9, sigma_beta = 0.9, sigma_gamma = 0.5, mu_hat = 1.3, sigma = 2, seed = 123) {
  set.seed(seed)
  # sigma_alpha <- 3
  # sigma_beta <- 3
  # sigma_gamma <- 2
  # sigma <- 0.01
  # mu_hat <- 1
  # n <- 3
  alpha <- rnorm(5, 0, sigma_alpha)
  beta <- rnorm(3, 0, sigma_beta)
  gamma <- rnorm(15, 0, sigma_gamma)

  effect_raw <- rep(alpha, each = 3) + rep(beta, 5) + gamma
  effect <- mapply(rnorm, n, effect_raw, sigma) |> as.numeric()
  mu <- mu_hat + effect

  list_data <- list()
  list_data$y <- mu
  list_data$J <- 5
  list_data$K <- 3
  list_data$JK <- 15
  list_data$N <- length(mu)
  list_data$jj <- rep(1:5, each = 3) |> rep(each = n)
  list_data$kk <- rep(1:3, 5) |> rep(each = n)
  list_data$jk <- rep(1:15, each = n)
  list_data$alpha <- alpha
  list_data$beta <- beta
  list_data$gamma <- gamma
  list_data
}

#' @title Create summary table for posteriors
create_stan_tab <- function(draws) {
  tmp <- draws |>
    janitor::clean_names() |>
    dplyr::select(contains(c("alpha", "beta", "gamma", "effect", "pred", "mu_hat")))
  mean_ <- apply(tmp, 2, mean)
  q2_5 <- apply(tmp, 2, \(x)(quantile(x, 0.025)))
  q5 <- apply(tmp, 2, \(x)(quantile(x, 0.05)))
  q97_5 <- apply(tmp, 2, \(x)(quantile(x, 0.975)))
  q95 <- apply(tmp, 2, \(x)(quantile(x, 0.9)))
  tibble(para = names(mean_), mean_, q2_5, q5, q95, q97_5)
}
