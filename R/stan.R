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

generate_piecewise_logistic_stan_data <- function(data) {
  # library(tidyverse)
  # library(targets)
  # tar_load(cond_count_csv)
  # data <- cond_count_csv
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

generate_logistic_stan_data <- function(data, quad = FALSE) {
  d <- read_csv(data) |>
    filter(!is.na(count))
  scale_d <- scale(d$pressure)
  stan_data <- list(
    N = nrow(d),
    J = unique(d$species) |> length(),
    K = 2,
    jj = as.factor(d$species) |> as.numeric(),
    y = d$count,
    total = d$total,
    x = cbind(1, scale_d),
    u = t(as.matrix(rep(1, 5)))
  )
  if (quad) {
    stan_data$K <- 3
    stan_data$x <- cbind(1, scale_d, scale_d^2)
  }
  stan_data
}

generate_dummy_data_ab <- function(n_measure = 6, n_tree = 9, n_sp = 20, n_xy = 4, seed = 123) {
  set.seed(seed)
  delta_p <- seq(0, 1, length = n_measure)
  # n_measure <- 6
  # n_tree <- 9
  # n_sp <- 20
  # n_xy <- 4

  tmp <- n_sp / n_xy
  tmp_end <- seq(tmp, 4*tmp, length = n_xy)
  tmp_start <- seq(1, 3*tmp + 1, length = n_xy)

  alpha0_xy <- rnorm(n_xy, 0, 1)
  alpha1_xy <- rnorm(n_xy, 0, 1)
  alpha2_xy <- rnorm(n_xy, 0, 1)

  beta0_xy <- rnorm(n_xy, 0, 1)
  beta1_xy <- rnorm(n_xy, 0, 1)
  beta2_xy <- rnorm(n_xy, 0, 1)

  rho_xy <- rnorm(n_xy, 1)
  alpha0_sp_raw <- mapply(rnorm, n_sp, alpha0_xy, 0.8)
  alpha1_sp_raw <- mapply(rnorm, n_sp, alpha1_xy, 0.8)
  alpha2_sp_raw <- mapply(rnorm, n_sp, alpha2_xy, 0.8)
  beta0_sp_raw <- mapply(rnorm, n_sp, beta0_xy, 0.8)
  beta1_sp_raw <- mapply(rnorm, n_sp, beta1_xy, 0.8)
  beta2_sp_raw <- mapply(rnorm, n_sp, beta2_xy, 0.8)
  rho_sp_raw <- mapply(rnorm, n_sp, rho_xy, 0.8)

  alpha0_sp <- c(
   alpha0_sp_raw[tmp_start[1]:tmp_end[1], 1],
   alpha0_sp_raw[tmp_start[2]:tmp_end[2], 2],
   alpha0_sp_raw[tmp_start[3]:tmp_end[3], 3],
   alpha0_sp_raw[tmp_start[4]:tmp_end[4], 4]
  )
  alpha1_sp <- c(
   alpha1_sp_raw[tmp_start[1]:tmp_end[1], 1],
   alpha1_sp_raw[tmp_start[2]:tmp_end[2], 2],
   alpha1_sp_raw[tmp_start[3]:tmp_end[3], 3],
   alpha1_sp_raw[tmp_start[4]:tmp_end[4], 4]
  )
  alpha2_sp <- c(
   alpha2_sp_raw[tmp_start[1]:tmp_end[1], 1],
   alpha2_sp_raw[tmp_start[2]:tmp_end[2], 2],
   alpha2_sp_raw[tmp_start[3]:tmp_end[3], 3],
   alpha2_sp_raw[tmp_start[4]:tmp_end[4], 4]
  )

  beta0_sp <- c(
   beta0_sp_raw[tmp_start[1]:tmp_end[1], 1],
   beta0_sp_raw[tmp_start[2]:tmp_end[2], 2],
   beta0_sp_raw[tmp_start[3]:tmp_end[3], 3],
   beta0_sp_raw[tmp_start[4]:tmp_end[4], 4]
  )
  beta1_sp <- c(
   beta1_sp_raw[tmp_start[1]:tmp_end[1], 1],
   beta1_sp_raw[tmp_start[2]:tmp_end[2], 2],
   beta1_sp_raw[tmp_start[3]:tmp_end[3], 3],
   beta1_sp_raw[tmp_start[4]:tmp_end[4], 4]
  )
  beta2_sp <- c(
   beta2_sp_raw[tmp_start[1]:tmp_end[1], 1],
   beta2_sp_raw[tmp_start[2]:tmp_end[2], 2],
   beta2_sp_raw[tmp_start[3]:tmp_end[3], 3],
   beta2_sp_raw[tmp_start[4]:tmp_end[4], 4]
  )

  rho_sp <- c(
   rho_sp_raw[tmp_start[1]:tmp_end[1], 1],
   rho_sp_raw[tmp_start[2]:tmp_end[2], 2],
   rho_sp_raw[tmp_start[3]:tmp_end[3], 3],
   rho_sp_raw[tmp_start[4]:tmp_end[4], 4]
  )

  alpha0_samp <- mapply(rnorm, n_tree, alpha0_sp, 0.2) |> as.numeric()
  alpha1_samp <- mapply(rnorm, n_tree, alpha1_sp, 0.2) |> as.numeric()

  beta0_samp <- mapply(rnorm, n_tree, beta0_sp, 0.2) |> as.numeric()
  beta1_samp <- mapply(rnorm, n_tree, beta1_sp, 0.2) |> as.numeric()

  rho_samp <- mapply(rnorm, n_tree, rho_sp, 0.2) |> as.numeric()

  l_samp <- runif(n_tree * n_sp, 0.8, 1.2)

  group_data <- tibble(
    tree = rep(1:n_tree, n_sp),
    sp = rep(1:n_sp, each = n_tree),
    xy = rep(1:n_xy, each = n_sp * n_tree / n_xy)) |>
    mutate(l_samp = l_samp) |>
    mutate(rho_samp = rho_samp) |>
    mutate(alpha0_samp = alpha0_samp) |>
    mutate(alpha1_samp = alpha1_samp) |>
    mutate(beta0_samp = beta0_samp) |>
    mutate(beta1_samp = beta1_samp) |>
    mutate(alpha2_sp = rep(alpha2_sp, each = n_tree)) |>
    mutate(beta2_sp = rep(beta2_sp, each = n_tree))


  data <- NULL
  for (i in 1:n_measure) {
    data_tmp <- group_data |>
      mutate(measure_id = i) |>
      mutate(p = delta_p[i])
    data <- bind_rows(data, data_tmp)
  }

  k <- rnorm(nrow(data), 0, 1)
  data <- data |>
    mutate(delta_p = p / l_samp) |>
    mutate(k = k) |>
    mutate(mu_a = alpha0_samp + alpha1_samp * delta_p + alpha2_sp * rho_samp) |>
    mutate(mu_b = beta0_samp + beta1_samp * delta_p + beta2_sp * rho_samp) |>
    mutate(a = rnorm(nrow(data), mu_a, 0.1)) |>
    mutate(b = rnorm(nrow(data), mu_b, 0.1)) |>
    mutate(fd = rnorm(nrow(data), a + b * k, 0.1)) |>
    mutate(sp_tree = str_c(sp, "_", tree))

  data2 <- data |>
    mutate(sp_tree = str_c(sp, "_", tree))  |>
    filter(!duplicated(sp_tree))



  sp_lab <- data |>
    dplyr::select(sp, xy) |>
    unique()

  tree_lab <- data |>
    dplyr::select(sp, tree) |>
    unique()


  list(
    data = data,
    Ni = nrow(data),
    Nj = data$sp_tree |> unique() |> length(),
    Nk = data$sp |> unique() |> length(),
    Nl = data$xy |> unique() |> length(),
    Mi = 2,
    Mj = 3,
    Ml = 1,
    Mk = 1,
    jj = data$sp_tree |> as.factor() |> as.numeric(),
    jj2 = tree_lab$sp |> as.factor() |> as.numeric(),
    kk = data$sp |> as.factor() |> as.numeric(),
    kk2 = sp_lab$xy |> as.factor() |> as.numeric(),
    ll = data$xy,
    # Xi = rbind(1, data$k),
    xj = rbind(1, data2$delta_p, data2$rho_samp),
    # Xj1 = cbind(data2$rho_samp),
    # xk = cbind(rep(1, 3)),
    # xl = cbind(rep(1, 3)),
    y = data$fd,
    x = data$k,
    delta_p = data$delta_p,
    rho = data$rho_samp
  ) #|>
  # str()
}

clean_sap_data <- function(data, file) {
  d <- read_csv(data) |>
    janitor::clean_names() |>
    rename(species = species_name)

  d <- d |>
    mutate(fd = ifelse(is.na(fd), removed_fd, fd)) |>
    mutate(k = ifelse(is.na(k), removed_k, k)) |>
    filter(!is.na(k)) |>
    filter(!is.na(fd)) |>
    filter(k > 0) |>
    # remove two acacia species
    filter(k < 2) |>
    filter(fd > 0) |>
    mutate(sample_id = str_c(species, "_", sample_number))

  my_write_csv(d, file)

}

generate_sap_stan_data <- function(data, remove_abnormal_values = FALSE, upper_pressure = FALSE) {
  # d <- read_csv("data-raw/calibration_raw_data.csv") |>
  #   janitor::clean_names() |>
  #   rename(species = species_name)
  # d <- read_csv("data/fd_k_traits.csv")
  d <- read_csv(data)
  if (remove_abnormal_values) {
    d <- d |>
      filter(is.na(removed_k))
  }

  if (upper_pressure) {
   d <- d |>
    filter(p_2 <= upper_pressure)
  }

  tmp <- d |>
    group_by(sample_id, species) |>
    nest() |>
    ungroup() |>
    arrange(sample_id)

  uj <- model.matrix(~ species, tmp)
  uj[apply(uj, 1, sum) == 2, 1] <- 0

  tmp <- d |>
    group_by(species, xylem_type) |>
    nest() |>
    ungroup() |>
    arrange(xylem_type) |>
    arrange(species)

  uk <- model.matrix(~ xylem_type, tmp)
  uk[apply(uk, 1, sum) == 2, 1] <- 0

  ul <- matrix(rep(1, 4), ncol = 4)

  stan_data <- list(
    N = nrow(d),
    J = unique(d$sample_id) |> length(),
    K = unique(d$species) |> length(),
    L = unique(d$xylem_type) |> length(),
    jj = as.factor(d$sample_id) |> as.numeric(),
    kk = as.factor(d$species) |> as.numeric(),
    ll = as.factor(d$xylem_type) |> as.numeric(),
    uj = t(uj),
    uk = t(uk),
    ul = ul,
    x = cbind(1, log(d$k)),
    y = log(d$fd)
  )

  stan_data
}


# alpha <- rbind(rep(1, 31), 1:31)
# alpha %*% sap_stan_data$uj

# beta <- rbind(rep(1, 4), 1:4)
# beta %*% sap_stan_data$uk

generate_sap_stan_data_sp <- function(data, remove_abnormal_values = FALSE, upper_pressure = FALSE) {
  # d <- read_csv("data-raw/calibration_raw_data.csv") |>
  #   janitor::clean_names() |>
  #   rename(species = species_name)
  d <- read_csv(data)
  if (remove_abnormal_values) {
    d <- d |>
      filter(is.na(removed_k))
  }

  if (upper_pressure) {
   d <- d |>
    filter(p_2 <= upper_pressure)
  }

  nd <- d |>
    group_by(species) |>
    nest()

  nd2 <- nd |>
    mutate(segment_number = map_dbl(data, \(x)length(unique(x$sample_id)))) |>
    mutate(stan_data = map2(data, segment_number, \(x, y) {
      list(
        N = nrow(x),
        log_fd = log(x$fd),
        log_k = log(x$k),
        u = matrix(rep(1, y), ncol = y),
        jj = as.factor(x$sample_id) |> as.numeric(),
        J = y,
        y = log(x$fd),
        x = cbind(1, log(x$k))
      )
    }))
  nd2
}

generate_sap_stan_data_segment <- function(data, remove_abnormal_values = FALSE, upper_pressure = FALSE) {
  # d <- read_csv("data-raw/calibration_raw_data.csv") |>
  #   janitor::clean_names() |>
  #   rename(species = species_name)
  # d <- read_csv("data/fd_k_traits.csv")
  d <- read_csv(data)
  if (remove_abnormal_values) {
    d <- d |>
      filter(is.na(removed_k))
  }

  if (upper_pressure) {
   d <- d |>
    filter(p_2 <= upper_pressure)
  }


  nd <- d |>
    group_by(species, sample_id) |>
    nest()

  nd2 <- nd |>
    mutate(stan_data = map(data, \(x) {
      list(
        N = nrow(x),
        log_fd = log(x$fd),
        log_k = log(x$k)
      )
    }))

  nd2
}



fit_model <- function(data, model_file,
                            iter_warmup = 1,
                            iter_sampling = 1,
                            adapt_delta = 0.9,
                            max_treedepth = 15,
                            chains = 4,
                            parallel_chains = 1,
                            refresh = 0,
                            seed = 123) {
  model <- cmdstan_model(model_file)
  fit <- model$sample(
    data = data,
    seed = seed,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    chains = chains,
    parallel_chains = parallel_chains,
    refresh = refresh)

  summary_ <- posterior::summarise_draws(fit,
  mean = ~mean(.x),
  sd = ~sd(.x),
  mad = ~mad(.x),
  ~posterior::quantile2(.x, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)),
  posterior::default_convergence_measures())

  list(summary = summary_, draws = fit$draws())
 }

#' @title Compile a Stan model and return a path to the compiled model output.
#' @description We return the paths to the Stan model specification
#'   and the compiled model file so `targets` can treat them as
#'   dynamic files and compile the model if either file changes.
#' @return Path to the compiled Stan model, which is just an RDS file.
#'   To run the model, you can read this file into a new R session with
#'   `readRDS()` and feed it to the `object` argument of `sampling()`.
#' @param model_file Path to a Stan model file.
#'   This is a text file with the model spceification.
#' @references https://github.com/wlandau/targets-stan
#' @examples
#' library(cmdstanr)
#' compile_model("stan/model.stan")
compile_model <- function(model_file) {
  quiet(cmdstan_model(model_file))
  model_file
}

# alpha <- rbind(rep(1, 31), 1:31)
# alpha %*% sap_stan_data$uj

# beta <- rbind(rep(1, 4), 1:4)
# beta %*% sap_stan_data$uk

#' @title Suppress output and messages for code.
#' @description Used in the pipeline.
#' @return The result of running the code.
#' @param code Code to run quietly.
#' @references https://github.com/wlandau/targets-stan
#' @examples
#' library(cmdstanr)
#' library(tidyverse)
#' compile_model("stan/model.stan")
#' quiet(fit_model("stan/model.stan", simulate_data_discrete()))
#' out
quiet <- function(code) {
  sink(nullfile())
  on.exit(sink())
  suppressMessages(code)
}
