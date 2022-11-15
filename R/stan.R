#' @title Get posterior estimates mcmc summary
#' @param data data frame, summary of mcmc
#' @param row variable name (e.g., "theta")
#' @param col summary name (e.g., "mean", "q50")
#' @param digits integer indicating the number of decimal places
#' @param nsmall the minimum number of digits to the right of the decimal point
get_post_para <- function(data, row, col, digits = 2, nsmall = 2) {
  data |>
    mutate_if(is.numeric, \(x) round(x, digits = digits)) |>
    mutate_if(is.numeric, \(x) format(x, nsmall = nsmall)) |>
    filter(variable == {{row}}) |>
    pull({{col}})
}

#' @title write_csv for targets
#' @inheritParams readr::write_csv
my_write_csv <- function(x, path, append = FALSE, col_names = !append) {
    write_csv(x, path, append = FALSE, col_names = !append)
    paste(path)
}

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
  p_g <- seq(0, 1, length = n_measure)
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
      mutate(p = p_g[i])
    data <- bind_rows(data, data_tmp)
  }

  k <- rnorm(nrow(data), 0, 1)
  data <- data |>
    mutate(p_g = p / l_samp) |>
    mutate(k = k) |>
    mutate(mu_a = alpha0_samp + alpha1_samp * p_g + alpha2_sp * rho_samp) |>
    mutate(mu_b = beta0_samp + beta1_samp * p_g + beta2_sp * rho_samp) |>
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
    xj = rbind(1, data2$p_g, data2$rho_samp),
    # Xj1 = cbind(data2$rho_samp),
    # xk = cbind(rep(1, 3)),
    # xl = cbind(rep(1, 3)),
    y = data$fd,
    x = data$k,
    p_g = data$p_g,
    rho = data$rho_samp
  ) #|>
  # str()
}

clean_sap_data <- function(data, file) {
  d <- read_csv(data) |>
    janitor::clean_names() |>
    rename(species = species_name) |>
    mutate(species = ifelse(species == "Bauhinia  tenuiflor", "Bauhinia tenuiflor", species)) |>
    mutate(species = ifelse(species == "Dypsis  lutescen", "Dypsis lutescen", species)) |>
    mutate(genus_short = str_sub(species, 1, 1)) |>
    mutate(sp_only = str_split_fixed(species, " ", 2)[, 2]) |>
    mutate(sp_short = str_c(genus_short, ". ", sp_only))

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

  d <- d |>
    filter(!sample_id %in% str_c("Areca catechu_", c(5, 6))) |>
    filter(!sample_id %in% str_c("Chrysalidocarpus lutescens_", c(2, 3)))

  my_write_csv(d, file)

}

generate_sap_traits_stan_data <- function(data, remove_abnormal_values = FALSE, upper_pressure = FALSE, trait_set = "all") {
  # library(tidyverse)
  # d <- read_csv("data/fd_k_traits.csv")
  d <- read_csv(data)

  d <- d |>
    filter(!is.na(wood_density)) |>
    filter(!is.na(swc)) |>
    filter(!is.na(dh)) |>
    filter(!is.na(vaf)) |>
    filter(!is.na(vf)) |>
    filter(!is.na(ks))

  if (remove_abnormal_values) {
    d <- d |>
      filter(is.na(removed_k))
  }

  if (upper_pressure) {
   d <- d |>
    filter(p_g <= upper_pressure)
  }

  tmp <- d |>
    group_by(sample_id, species) |>
    nest() |>
    ungroup() |>
    arrange(sample_id)

  uj <- model.matrix(~ species, tmp)
  uj[apply(uj, 1, sum) == 2, 1] <- 0

  tmp0 <- d |>
    mutate(log_swc = log(swc)) |>
    mutate(log_dh = log(dh)) |>
    mutate(log_vaf = log(vaf)) |>
    mutate(log_vf = log(vf)) |>
    mutate(log_ks = log(ks)) |>
    group_by(sample_id) |>
    summarise_if(is.numeric, mean, na.rm = TRUE)

  if (trait_set == "all") {
    tmp <- tmp0 |>
      dplyr::select(wood_density, log_dh, log_vaf, log_vf, log_ks)
  } else if (trait_set == "nowd") {
    tmp <- tmp0 |>
      dplyr::select(log_dh, log_vaf, log_vf, log_ks)
  } else if (trait_set == "novf") {
    tmp <- tmp0 |>
      dplyr::select(log_dh, log_vaf, log_ks)
  } else if (trait_set == "novaf") {
    tmp <- tmp0 |>
      dplyr::select(log_dh, log_vf, log_ks)
  } else if (trait_set == "ks") {
    tmp <- tmp0 |>
      dplyr::select(log_ks)
  } else if (trait_set == "wd") {
    tmp <- tmp0 |>
      dplyr::select(wood_density)
  } else if (trait_set == "dh") {
    tmp <- tmp0 |>
      dplyr::select(log_dh)
  } else if (trait_set == "vaf") {
    tmp <- tmp0 |>
      dplyr::select(log_vaf)
  } else if (trait_set == "vf") {
    tmp <- tmp0 |>
      dplyr::select(log_vf)
  }

  tmp2 <- apply(tmp, 2, scale)
  #tmp2 <- na.omit(tmp) |> as.matrix()
  xj <- cbind(1, tmp2)

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
    y = log(d$fd),
    xj = xj,
    T = ncol(xj)
  )

  stan_data
}

generate_sap_stan_data <- function(data, remove_abnormal_values = FALSE, upper_pressure = FALSE, traits = FALSE) {
  # library(tidyverse)
  # d <- read_csv("data/fd_k_traits.csv")
  d <- read_csv(data)

  if (traits) {
    d <- d |>
      filter(!is.na(wood_density)) |>
      filter(!is.na(swc)) |>
      filter(!is.na(dh)) |>
      filter(!is.na(vaf)) |>
      filter(!is.na(vf)) |>
      filter(!is.na(ks))
  }

  if (remove_abnormal_values) {
    d <- d |>
      filter(is.na(removed_k))
  }

  if (upper_pressure) {
   d <- d |>
    filter(p_g <= upper_pressure)
  }

  tmp <- d |>
    group_by(sample_id, species) |>
    nest() |>
    ungroup() |>
    arrange(sample_id)

  uj <- model.matrix(~ species, tmp)
  uj[apply(uj, 1, sum) == 2, 1] <- 0

  tmp <- d |>
    group_by(sample_id) |>
    summarise_if(is.numeric, mean, na.rm = TRUE) |>
    dplyr::select(wood_density, swc, dh, vaf, vf, ks)

  tmp2 <- na.omit(tmp) |> as.matrix()
  xj <- cbind(1, tmp2)

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

  if (traits) {
    stan_data$xj <- xj
    stan_data$T <- ncol(xj)
  }

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
    filter(p_g <= upper_pressure)
  }

  nd <- d |>
    group_by(species) |>
    nest() |>
    ungroup() |>
    arrange(species)

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
    filter(p_g <= upper_pressure)
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

generate_ab_var_data <- function(path, draws) {
  # draws <- withr::with_dir(rprojroot::find_root('_targets.R'),
  #   targets::tar_read(fit_ab_draws_granier_without_traits_sap_all_clean_0.08)) |>
  #   janitor::clean_names()
  tmp <- draws |>
      janitor::clean_names()  |>
      dplyr::select(matches("tau|sigma")) |>
      mutate(var_a = tau_j_1^2 + tau_k_1^2 + tau_l_1^2) |>
      mutate(var_b = tau_j_2^2 + tau_k_2^2 + tau_l_2^2) |>
      mutate(var_a_segment = tau_j_1^2 / var_a) |>
      mutate(var_a_sp = tau_k_1^2 / var_a) |>
      mutate(var_a_xylem = tau_l_1^2 / var_a) |>
      mutate(var_b_segment = tau_j_2^2 / var_b) |>
      mutate(var_b_sp = tau_k_2^2 / var_b) |>
      mutate(var_b_xylem = tau_l_2^2 / var_b)

  tibble(
    variable = colnames(tmp),
    mean = apply(tmp * 100, 2, mean),
    q2.5 = apply(tmp * 100, 2, quantile, 0.025),
    q25 = apply(tmp * 100, 2, quantile, 0.25),
    q50 = apply(tmp * 100, 2, quantile, 0.5),
    q75 = apply(tmp * 100, 2, quantile, 0.75),
    q97.5 = apply(tmp * 100, 2, quantile, 0.975)
  ) |>
  my_write_csv(path)
}


coef_ab_vpart <- function(draws) {
   tmp <- draws |>
      janitor::clean_names()  |>
      dplyr::select(matches("tau|sigma")) |>
      mutate(var_a = tau_j_1^2 + tau_k_1^2 + tau_l_1^2) |>
      mutate(var_b = tau_j_2^2 + tau_k_2^2 + tau_l_2^2) |>
      mutate(var_a_segment = tau_j_1^2 / var_a) |>
      mutate(var_a_sp = tau_k_1^2 / var_a) |>
      mutate(var_a_xylem = tau_l_1^2 / var_a) |>
      mutate(var_b_segment = tau_j_2^2 / var_b) |>
      mutate(var_b_sp = tau_k_2^2 / var_b) |>
      mutate(var_b_xylem = tau_l_2^2 / var_b)

}

my_loo <- function(x) x$loo(cores = parallel::detectCores())

#' @title Check divergence from draws
div_check <- function(diags) {
  n1 <- diags |>
    filter(divergent__ == 1) |>
    nrow()
  # n2 <- diags |>
    nrow()
  print(paste(
    n1, "of", n2,
    "iterations ended with a divergence", n1 / n2 * 100, "%"
  ))
}



#' @title with trait csv
write_without_traits_csv <- function(stan_summary, output) {

  d <- read_csv(here("data/fd_k_traits.csv")) |>
  filter(is.na(removed_k))

  # s <- withr::with_dir(rprojroot::find_root('_targets.R'),
  #   targets::tar_read(stan_summary))
  s <- stan_summary

  sample_id <- tibble(sample_id = unique(d$sample_id) |> as.factor()) |>
    arrange(sample_id) |>
    mutate(sample_id_num = as.numeric(sample_id)) |>
    mutate(sample_id = as.character(sample_id))

  species <- tibble(species = unique(d$species) |> as.factor()) |>
    arrange(species) |>
    mutate(species_num = as.numeric(species)) |>
    mutate(species = as.character(species))

  xylem <- tibble(xylem = unique(d$xylem_type) |> as.factor()) |>
    arrange(xylem) |>
    mutate(xylem_num = as.numeric(xylem)) |>
    mutate(xylem = as.character(xylem))

  para_name <- c("a", "b")
  gamma <- s |>
    filter(str_detect(variable, "gamma")) |>
    mutate(para = variable)

  beta <- s |>
    filter(str_detect(variable, "beta")) |>
    mutate(para = variable) #|>

  for (i in 1:nrow(beta)) {
    beta <- beta |>
      mutate(para = case_when(
        str_detect(para, str_c("\\[", i)) ~ str_replace(para, str_c("\\[",i),
          str_c("_", as.character(para_name[i]))),
      TRUE ~ para)) |>
      mutate(para = case_when(
        str_detect(para, str_c(",", i, "\\]")) ~ str_replace(para, str_c("," ,i, "\\]"),
          str_c("_", as.character(xylem[i, 1]))),
      TRUE ~ para))
  }


  alpha <- s |>
    filter(str_detect(variable, "alpha")) |>
    mutate(para = variable) #|>
    # mutate(para = str_replace_all(para, "alpha", "trait"))
  for (i in 1:nrow(alpha)) {
    alpha <- alpha |>
      mutate(para = case_when(
        str_detect(para, str_c("\\[", i)) ~ str_replace(para, str_c("\\[",i),
          str_c("_", as.character(para_name[i]))),
      TRUE ~ para)) |>
      mutate(para = case_when(
        str_detect(para, str_c(",", i, "\\]")) ~ str_replace(para, str_c("," ,i, "\\]"),
          str_c("_", as.character(species[i, 1]))),
      TRUE ~ para))
  }

  A <- s |>
    filter(str_detect(variable, "A")) |>
    mutate(para = variable) #|>

  for (i in 1:nrow(A)) {
    A <- A |>
      # mutate(para = case_when(
      #   str_detect(para, str_c("\\[", i)) ~ str_replace(para, str_c("\\[",i),
      #     str_c("_", as.character(para_name[i]))),
      # TRUE ~ para)) |>
      mutate(para = case_when(
        str_detect(para, str_c(",", i, "\\]")) ~ str_replace(para, str_c("," ,i, "\\]"),
          str_c("_", as.character(sample_id[i, 1]))),
      TRUE ~ para)) |>
      mutate(para =case_when(
        str_detect(para, "A\\[1") ~ str_replace(para, "A\\[1", "a"),
        str_detect(para, "A\\[2") ~ str_replace(para, "A\\[2", "b"),
        TRUE ~ para
      ))
  }

  bind_rows(gamma, beta, alpha, A) |>
    my_write_csv(output)

}

#' @title with trait csv
write_with_traits_csv <- function(stan_summary, output) {
  s <- stan_summary
  # d <- withr::with_dir(rprojroot::find_root('_targets.R'),
  #   targets::tar_read(fd_k_traits_csv)) |>
  # d <- fd_k_traits_csv |>
  #   here() |>
  #   read_csv() |>
  d <- read_csv(here("data/fd_k_traits.csv")) |>
    filter(is.na(removed_k)) |>
    filter(!is.na(wood_density)) |>
    filter(!is.na(swc)) |>
    filter(!is.na(dh)) |>
    filter(!is.na(vaf)) |>
    filter(!is.na(vf)) |>
    filter(!is.na(ks))

sample_id <- tibble(sample_id = unique(d$sample_id) |> as.factor()) |>
  arrange(sample_id) |>
  mutate(sample_id_num = as.numeric(sample_id)) |>
  mutate(sample_id = as.character(sample_id))

species <- tibble(species = unique(d$species) |> as.factor()) |>
  arrange(species) |>
  mutate(species_num = as.numeric(species)) |>
  mutate(species = as.character(species))

xylem <- tibble(xylem = unique(d$xylem_type) |> as.factor()) |>
  arrange(xylem) |>
  mutate(xylem_num = as.numeric(xylem)) |>
  mutate(xylem = as.character(xylem))

traits <- tibble(trait = c("intercept", "wood_density", "log_dh", "log_vf"))

gamma <- s6 |>
  filter(str_detect(variable, "gamma")) |>
  mutate(para = variable) #|>
  # mutate(para = str_replace_all(para, "gamma", "trait"))
for (i in 1:4) {
  gamma <- gamma |>
    mutate(para = case_when(
      str_detect(para, str_c("\\[", i)) ~ str_replace(para, str_c("\\[",i),
        str_c("_", as.character(traits[i, 1]))),
    TRUE ~ para)) |>
    mutate(para = str_remove(para, ",1\\]"))
}

beta <- s6 |>
  filter(str_detect(variable, "beta")) |>
  mutate(para = variable) #|>
  # mutate(para = str_replace_all(para, "beta", "trait"))
for (i in 1:nrow(beta)) {
  beta <- beta |>
    mutate(para = case_when(
      str_detect(para, str_c("\\[", i)) ~ str_replace(para, str_c("\\[",i),
        str_c("_", as.character(traits[i, 1]))),
    TRUE ~ para)) |>
    mutate(para = case_when(
      str_detect(para, str_c(",", i, "\\]")) ~ str_replace(para, str_c("," ,i, "\\]"),
        str_c("_", as.character(xylem[i, 1]))),
    TRUE ~ para))
}

alpha <- s6 |>
  filter(str_detect(variable, "alpha")) |>
  mutate(para = variable) #|>
  # mutate(para = str_replace_all(para, "alpha", "trait"))
for (i in 1:nrow(alpha)) {
  alpha <- alpha |>
    mutate(para = case_when(
      str_detect(para, str_c("\\[", i)) ~ str_replace(para, str_c("\\[",i),
        str_c("_", as.character(traits[i, 1]))),
    TRUE ~ para)) |>
    mutate(para = case_when(
      str_detect(para, str_c(",", i, "\\]")) ~ str_replace(para, str_c("," ,i, "\\]"),
        str_c("_", as.character(species[i, 1]))),
    TRUE ~ para))
}


a_hat <- s6 |>
  filter(str_detect(variable, "a_hat")) |>
  mutate(para = variable) #|>
for (i in 1:nrow(a_hat)) {
  a_hat <- a_hat |>
    mutate(para = case_when(
      str_detect(para, str_c("\\[", i, "\\]")) ~ str_replace(para, str_c("\\[",i, "\\]"),
        str_c("_", as.character(sample_id[i, 1]))),
    TRUE ~ para))
}
b_hat <- s6 |>
  filter(str_detect(variable, "b_hat")) |>
  mutate(para = variable) #|>
for (i in 1:nrow(b_hat)) {
  b_hat <- b_hat |>
    mutate(para = case_when(
      str_detect(para, str_c("\\[", i, "\\]")) ~ str_replace(para, str_c("\\[",i, "\\]"),
        str_c("_", as.character(sample_id[i, 1]))),
    TRUE ~ para))
}

A_hat <- s6 |>
  filter(str_detect(variable, "A_hat")) |>
  mutate(para = variable) #|>
for (i in 1:nrow(A_hat)) {
  A_hat <- A_hat |>
    mutate(para = case_when(
      str_detect(para, str_c(",", i, "\\]")) ~ str_replace(para, str_c(",",i, "\\]"),
        str_c("_", as.character(sample_id[i, 1]))),
    TRUE ~ para)) |>
    mutate(para =case_when(
      str_detect(para, "A_hat\\[1") ~ str_replace(para, "A_hat\\[1", "a"),
      str_detect(para, "A_hat\\[2") ~ str_replace(para, "A_hat\\[2", "b"),
      TRUE ~ para
    ))
}
A_hat

bind_rows(gamma, beta, alpha, a_hat, b_hat, A_hat) |>
  my_write_csv(ouput)

}

get_tmp <- function(data, row, col, chr = FALSE) {
  out <- data |>
    filter(variable == {{row}}) |>
    pull({{col}})
  if (chr) {
    round(out, 2) |> format(nsmall = 2)
  } else {
    out
  }
}

get_tmp2 <- function(data, row, col, chr = FALSE) {
  out <- data |>
    filter(str_detect(variable, {{row}})) |>
    pull({{col}})
  if (chr) {
    round(out, 2) |> format(nsmall = 2)
  } else {
    out
  }
}

get_tmp2_chr <- function(data, col) {
  lwr <- get_tmp2(data, col, q2.5, chr = TRUE)
  mid <- get_tmp2(data, col, q50, chr = TRUE)
  upr <- get_tmp2(data, col, q97.5, chr = TRUE)
  str_c(mid, " [", lwr, ", ", upr, "]")
}

get_tmp_chr <- function(data, col) {
  lwr <- get_tmp(data$summary, col, q2.5, chr = TRUE)
  mid <- get_tmp(data$summary, col, q50, chr = TRUE)
  upr <- get_tmp(data$summary, col, q97.5, chr = TRUE)
  str_c(mid, " [", lwr, ", ", upr, "]")
}

write_ab_csv <- function(d, summary_full_pool, summary_full_segments, summary_sp,
                         out,
  with_traits = FALSE) {

  d <- read_csv(d) |>
    filter(is.na(removed_k))

  summary_sp2 <- summary_sp |>
    mutate(log_a_sp_pool = map_chr(fit_pool, get_tmp_chr, "log_a")) |>
    mutate(b_sp_pool = map_chr(fit_pool, get_tmp_chr, "b")) |>
    mutate(sigma_sp_pool = map_chr(fit_pool, get_tmp_chr, "sigma")) |>
    mutate(log_a_sp_segments = map_chr(fit_segments, get_tmp_chr, "log_a")) |>
    mutate(b_sp_segments = map_chr(fit_segments, get_tmp_chr, "b")) |>
    mutate(sigma_sp_segments = map_chr(fit_segments, get_tmp_chr, "sigma")) |>
    dplyr::select(species, log_a_sp_pool:sigma_sp_segments)

  d2 <- d |>
    group_by(species, sp_short, xylem_type) |> nest() |>
    ungroup() |>
    arrange(species) |>
    mutate(log_xx = map(data, \(x)seq(min(log(x$k)), max(log(x$k)), length = 80))) |>
    mutate(log_a_full_pool = get_tmp2_chr(summary_full_pool, "alpha\\[1")) |>
    mutate(b_full_pool = get_tmp2_chr(summary_full_pool, "alpha\\[2")) |>
    mutate(sigma_full_pool = get_tmp2_chr(summary_full_pool, "sigma")) |>
    mutate(log_a_full_segments = get_tmp2_chr(summary_full_segments, "alpha\\[1")) |>
    mutate(b_full_segments = get_tmp2_chr(summary_full_segments, "alpha\\[2")) |>
    mutate(sigma_full_segments = get_tmp2_chr(summary_full_segments, "sigma")) |>
    dplyr::select(xylem_type, species, log_a_full_pool:sigma_full_segments)

  full_join(summary_sp2, d2, by = "species") |>
    my_write_csv(out)

}


generate_trait_fig_data <- function(summary_data, draws, fd_k_traits_csv, trait_name) {
# summary_data <- tar_read(fit_abt_summary_granier_with_traits_sap_trait_clean_vaf)
# draws <- tar_read(fit_abt_draws_granier_with_traits_sap_trait_clean_vaf)
  a_mat <- summary_data |>
    filter(str_detect(variable, "^A\\[1"))
  b_mat <- summary_data |>
    filter(str_detect(variable, "^A\\[2"))

# summary_data <- tar_read(fit_abt_summary_granier_with_traits_sap_trait_clean_vaf)
# summary_data |>
#     filter(str_detect(variable, "gamma"))

  d <- read_csv(fd_k_traits_csv)
  d <- d |>
    filter(!is.na(wood_density)) |>
    filter(!is.na(swc)) |>
    filter(!is.na(dh)) |>
    filter(!is.na(vaf)) |>
    filter(!is.na(vf)) |>
    filter(!is.na(ks))

  d <- d |>
    filter(is.na(removed_k)) #|>

  tmp0 <- d |>
    mutate(log_swc = log(swc)) |>
    mutate(log_dh = log(dh)) |>
    mutate(log_vaf = log(vaf)) |>
    mutate(log_vf = log(vf)) |>
    mutate(log_ks = log(ks)) |>
    group_by(sample_id, xylem_type, species) |>
    summarise_if(is.numeric, mean, na.rm = TRUE) |>
    ungroup() |>
    dplyr::select(sample_id, xylem_type, species,
     wood_density, log_dh, log_vaf, log_vf, log_ks)

  draws <- janitor::clean_names(draws)

  coef_a <- draws |>
    dplyr::select(gamma_a_1_1, gamma_a_2_1) |>
    as.matrix()
  coef_b <- draws |>
    dplyr::select(gamma_b_1_1, gamma_b_2_1) |>
    as.matrix()


  trait <- tmp0 |> pull({{trait_name}})
  ts <- scale(trait) |> range()
  ts <- seq(ts[1], ts[2], length = 100)
  xx <- sd(trait) * ts + mean(trait)
  # # ts <- scale(tmp0$log_vaf) |> range()
  # ts <- scale(tmp0 |> pull(trait)) |> range()
  # ts <- seq(ts[1], ts[2], length = 100)
  # xx <- sd(tmp0$log_vaf) * ts + mean(tmp0$log_vaf)


  pred_a <- coef_a %*% t(cbind(1, ts))
  pred_a_m <- apply(pred_a, 2, median)
  pred_a_ll <- apply(pred_a, 2, quantile, 0.025)
  pred_a_l <- apply(pred_a, 2, quantile, 0.25)
  pred_a_h <- apply(pred_a, 2, quantile, 0.75)
  pred_a_hh <- apply(pred_a, 2, quantile, 0.975)

  pred_b <- coef_b %*% t(cbind(1, ts))
  pred_b_m <- apply(pred_b, 2, median)
  pred_b_ll <- apply(pred_b, 2, quantile, 0.025)
  pred_b_l <- apply(pred_b, 2, quantile, 0.25)
  pred_b_h <- apply(pred_b, 2, quantile, 0.75)
  pred_b_hh <- apply(pred_b, 2, quantile, 0.975)

  pred_line <- tibble(
    pred_a_m, pred_a_ll, pred_a_l, pred_a_h, pred_a_hh,
    pred_b_m, pred_b_ll, pred_b_l, pred_b_h, pred_b_hh, x = exp(xx))

  pred_points <- tmp0 |>
    mutate(a_mid = a_mat$q50) |>
    mutate(a_lwr = a_mat$q2.5) |>
    mutate(a_upr = a_mat$q97.5) |>
    mutate(b_mid = b_mat$q50) |>
    mutate(b_lwr = b_mat$q2.5) |>
    mutate(b_upr = b_mat$q97.5)

  list(pred_points = pred_points, pred_line = pred_line)
}


traits_points <- function(vaf_pred_data, ks_pred_data) {
# tar_load(vaf_pred_data)
# tar_load(ks_pred_data)
 fig_fun <- function(data, trait_name, coef_a = TRUE) {
   if (coef_a) {
     ggplot() +
        geom_line(data = data$pred_line, aes(x = x, y = pred_a_m)) +
        geom_ribbon(data = data$pred_line, aes(x = x, ymin = pred_a_ll, ymax = pred_a_hh), alpha = 0.3, fill = "grey") +
        geom_ribbon(data = data$pred_line, aes(x = x, ymin = pred_a_l, ymax = pred_a_h), alpha = 0.8, fill = "grey") +
        geom_point(data = data$pred_points, aes(y = a_mid, x = exp({{trait_name}}), col = xylem_type), alpha = 0.6) +
        geom_errorbar(data = data$pred_points, aes(ymin = a_lwr, ymax = a_upr, x = exp({{trait_name}}), col = xylem_type)) +
        ylab(expression(Coefficient~italic(a))) +
        scale_x_log10() +
        my_theme() +
        theme(legend.position = "none")
   } else {
     ggplot() +
        geom_line(data = data$pred_line, aes(x = x, y = pred_b_m)) +
        geom_ribbon(data = data$pred_line, aes(x = x, ymin = pred_b_ll, ymax = pred_b_hh), alpha = 0.2, fill = "grey") +
        geom_ribbon(data = data$pred_line, aes(x = x, ymin = pred_b_l, ymax = pred_b_h), alpha = 0.8, fill = "grey") +
        geom_point(data = data$pred_points, aes(y = b_mid, x = exp({{trait_name}}), col = xylem_type)) +
        geom_errorbar(data = data$pred_points, aes(ymin = b_lwr, ymax = b_upr, x = exp({{trait_name}}), col = xylem_type)) +
        ylab(expression(Coefficient~italic(b))) +
        scale_x_log10() +
        my_theme() +
        theme(legend.position = "none")
   }
  }

  p1 <- fig_fun(vaf_pred_data, log_vaf) +
    xlab("VAF (%)")
  p2 <- fig_fun(ks_pred_data, log_ks) +
    xlab(expression(K[s]~(kg~m^{-1}~s^{-1}~MPa^{-1}))) +
    theme(legend.position = c(0.2, 0.75))
  p3 <- fig_fun(vaf_pred_data, log_vaf, coef_a = FALSE) +
    xlab("VAF (%)")
  p4 <- fig_fun(ks_pred_data, log_ks, coef_a = FALSE) +
    xlab(expression(K[s]~(kg~m^{-1}~s^{-1}~MPa^{-1})))

  p1 + p2 + p3 + p4 +
    plot_annotation(tag_levels = "A") #3
}
