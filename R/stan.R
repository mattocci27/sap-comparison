
generate_anova_data <- function(data, log = FALSE) {
  data <- read_csv(data)
  list_data <- list(
    N = nrow(data),
    jj = data$species |>  as.factor() |> as.numeric(),
    kk = data$pressure |>  as.factor() |> as.numeric(),
    jk = paste(data$species, data$pressure, sep = "_")
      |> as.factor() |> as.numeric(),
    y = data$pres_calib - data$tens_calib
  )

  if (log) {
    list_data$y <- log(data$pres_calib / data$tens_calib)
  }

  list_data$J <- unique(data$species) |> length()
  list_data$K <- unique(data$pressure) |> length()
  list_data$JK <- unique(list_data$jk) |> length()
  list_data
}

#' @title Create summary table for posteriors
create_stan_tab <- function(draws) {
  tmp <- draws |>
    janitor::clean_names() |>
    dplyr::select(contains(c("alpha", "beta", "gamma", "pred", "mu_hat")))
  mean_ <- apply(tmp, 2, mean)
  q2_5 <- apply(tmp, 2, \(x)(quantile(x, 0.025)))
  q5 <- apply(tmp, 2, \(x)(quantile(x, 0.05)))
  q97_5 <- apply(tmp, 2, \(x)(quantile(x, 0.975)))
  q95 <- apply(tmp, 2, \(x)(quantile(x, 0.9)))
  tibble(para = names(mean_), mean_, q2_5, q5, q95, q97_5)
}
