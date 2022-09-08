
my_theme <- function(){
  theme_bw() %+replace%
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.y = unit(0.1, "cm"),
    legend.text.align = 0,
    # legend.key.height = unit(0.2, "cm"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9)
  )
}

# theme_set(my_theme)
sma_scatter <- function(data, log = FALSE) {
#  targets::tar_load(five_spp_csv)
#  data <- five_spp_csv
  data <- read_csv(data)
  my_col <- RColorBrewer::brewer.pal(3, "Set2")
  p <- data |>
    mutate(pres_fac = paste(pressure, "MPa") |> as.factor()) |>
    ggplot(aes(x = pres_calib, y = tens_calib)) +
    geom_abline(slope = 1, intercept = 0, lty = 2, col = "grey60") +
    geom_sma(se = TRUE) +
    geom_point(aes(col = pres_fac)) +
    ylab(expression("Flow rate under tension "~(g~s^{-1}))) +
    xlab(expression("Flow rate by pressure "~(g~s^{-1}))) +
    scale_colour_manual(
      values = my_col,
      name = "") +
    stat_cor(
        label.x.npc = 0.35,
        label.y.npc = 0.1,
        vjust = 1,
        p.accuracy = 0.0001,
        aes(label = paste(..rr.label.., ..p.label.. ,sep = "~`,`~"))
      ) +
    my_theme() +
    theme(
      #plot.margin = unit(c(0, 0, 0, 0), "npc"),
      legend.position = c(0.25, 0.8)
    )
   if (log) {
     p +
      annotate(geom = "text", x = 3.5, y = 6,
        label = "1:1 line", angle = 45, col = "grey60") +
      # scale_x_log10(limits = c(min(data$pres_calib, data$tens_calib), max(data$pres_calib, data$tens_calib))) +
      # scale_y_log10(limits = c(min(data$pres_calib, data$tens_calib), max(data$pres_calib, data$tens_calib))) +
      scale_x_log10(limits = c(0.1, max(data$pres_calib, data$tens_calib))) +
      scale_y_log10(limits = c(0.1, max(data$pres_calib, data$tens_calib))) +
      coord_fixed()
   } else {
     p +
      annotate(geom = "text", x = 8, y = 9,
        label = "1:1 line", angle = 45, col = "grey60") +
      coord_cartesian(
        xlim = c(0,max(c(data$tens_calib, data$pres_calib))),
        ylim = c(0,max(c(data$tens_calib, data$pres_calib)))) +
      coord_fixed()
   }
}


ks_bars <- function(data, pred_draws) {
  data2 <- read_csv(data) |>
    mutate(sp = as.factor(species)) |>
    mutate(pres = as.factor(pressure)) |>
    mutate(jk = paste(species, pressure, sep = "_")) |>
    mutate(jk_fct = as.factor(jk) |> as.numeric()) |>
    dplyr::select(sp, pres, jk, jk_fct) |>
    unique() |>
    arrange(jk_fct)
  pred_draws2 <- pred_draws |>
    filter(str_detect(para, "effect")) |>
    bind_cols(data2)
  ggplot(pred_draws2, aes(x = pres)) +
    geom_point(aes(y = mean_)) +
    geom_errorbar(aes(ymin = q2_5, ymax = q97_5), width = 0) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
    facet_grid(sp ~ .) +
    my_theme()
}

ks_box <- function(data) {
  my_col <- RColorBrewer::brewer.pal(3, "RdBu")
  data <- read_csv(data)
  data |>
    mutate(calib_type = ifelse(pres_type == "pres",
     "Pressure", "Tension")) |>
    mutate(species = factor(species,
     levels = c("HH", "VM", "HB", "TG", "AP"))) |>
    ggplot(aes(y = ks, x = as.factor(pressure), fill = calib_type)) +
    #geom_violin() +
    geom_boxplot() +
    geom_point(position = position_jitterdodge(jitter.width = 0), alpha = 0.5) +
    # facet_wrap( ~ species, scale = "free") +
    facet_grid(rows = vars(species), scale = "free") +
    scale_y_log10() +
    ylab(expression(K[s]~(kg~m^{-1}~s^{-1}~MPa^{-1}))) +
    xlab("Pressure (MPa)") +
    scale_fill_manual(
      name = "Calibration",
      values = my_col[-2],
      ) +
    my_theme() +
    theme(
      #plot.margin = unit(c(0, 0, 0, 0), "cm"),
      # legend.margin =  unit(c(-1, 0, 0, 0), "cm"),
      #plot.background = element_rect(fill='#e3fbff'),
      strip.text = element_text(size = 8),
      legend.position = "top")
}

ks_box2 <- function(data) {
  my_col <- RColorBrewer::brewer.pal(3, "RdBu")
  data <- read_csv(data)
  data |>
    mutate(calib_type = ifelse(pres_type == "pres",
     "Pressure", "Tension")) |>
    mutate(species = factor(species,
     levels = c("HH", "VM", "HB", "TG", "AP"))) |>
    ggplot(aes(y = ks, x = as.factor(pressure), fill = calib_type)) +
    #geom_violin() +
    geom_boxplot() +
    geom_point(position = position_jitterdodge(jitter.width = 0), alpha = 0.5) +
    facet_wrap( ~ species, scale = "free") +
   # facet_grid(rows = vars(species), scale = "free") +
    scale_y_log10() +
    ylab(expression(K[s]~(kg~m^{-1}~s^{-1}~MPa^{-1}))) +
    xlab("Pressure (MPa)") +
    scale_fill_manual(
      name = "Calibration",
      values = my_col[-2],
      ) +
    my_theme() +
    theme(
      #plot.margin = unit(c(0, 0, 0, 0), "cm"),
      # legend.margin =  unit(c(-1, 0, 0, 0), "cm"),
      #plot.background = element_rect(fill='#e3fbff'),
      strip.text = element_text(size = 8),
      legend.position = "top")
}

sma_ks <- function(p1, p2) {
  p1 + p2 +
    plot_layout(nrow = 1, width = c(1.8, 1)) +
    plot_annotation(tag_levels = "a") &
    theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12)
    )
}

sma_ks2 <- function(p1, p2) {
  p1 + p2 +
    plot_layout(nrow = 1, width = c(1, 2)) +
    plot_annotation(tag_levels = "a") &
    theme(
      # plot.margin = unit(c(0, 0, 0, 0), "cm"),
      axis.title = element_text(size = 9)
    )
}


coef_intervals_sd <- function(draws) {
  intervals_data <- mcmc_intervals_data(
    draws,
    pars = c("sigma"),
    point_est = "median",
    prob = 0.5,
    prob_outer = 0.95,
    regex_pars = "tau\\[[1-9]\\]") |>
    mutate(para = case_when(
      parameter == "sigma" ~ "Residuals",
      parameter == "tau[1]" ~ "Species",
      parameter == "tau[2]" ~ "Pressure",
    )) |>
    mutate(para = factor(para,
      levels = c("Species", "Pressure", "Residuals") |> rev()))

  ggplot(intervals_data, aes(y = para)) +
    geom_linerange(aes(xmin = ll, xmax = hh)) +
    geom_linerange(aes(xmin = l, xmax = h), size = 2) +
    geom_point(aes(x = m), size = 3) +
    ylab("") +
    xlab("Standard deviation") +
    my_theme()
}

coef_intervals_mean <- function(draws) {
  intervals_data <- mcmc_intervals_data(
    draws,
    pars = c("mu_hat"),
    point_est = "median",
    prob = 0.5,
    prob_outer = 0.95,
    regex_pars = "alpha\\[[1-9]\\]|beta\\[[1-9]\\]") |>
    mutate(para = case_when(
      parameter == "beta[1]" ~ "p_2",
      parameter == "beta[2]" ~ "p_5",
      parameter == "beta[3]" ~ "p_8",
      parameter == "alpha[1]" ~ "AP",
      parameter == "alpha[2]" ~ "HB",
      parameter == "alpha[3]" ~ "HH",
      parameter == "alpha[4]" ~ "TG",
      parameter == "alpha[5]" ~ "VM",
      parameter == "mu_hat" ~ "ks_ratio",
    )) |>
    mutate(para = factor(para,
      levels = c("ks_ratio",
      "AP", "HB", "HH", "TG", "VM",
      "p_2", "p_5", "p_8"
      ) |> rev()))

  ggplot(intervals_data, aes(y = para)) +
    geom_vline(xintercept = 0, lty = 2, col = "grey60") +
    geom_linerange(aes(xmin = ll, xmax = hh)) +
    geom_linerange(aes(xmin = l, xmax = h), size = 2) +
    geom_point(aes(x = m), size = 3) +
    ylab("") +
    xlab("Posterior estimates") +
    my_theme()
}

coef_intervals_diff <- function(draws) {
  tmp <- tibble(x = 1:5, sp_x = c("AP",
      "HB",
      "HH",
      "TG",
      "VM"))
  tmp2 <- tmp |>
    rename(y = x) |>
    rename(sp_y = sp_x)

  tmp3 <- expand_grid(x = 1:5, y = 1:5) |>
    full_join(tmp) |>
    full_join(tmp2) |>
    mutate(para = paste0("alpha_", x, y)) |>
    mutate(para_name = paste(sp_x, "-", sp_y))

  tmp4 <- tribble(~ para, ~ para_name,
    "beta_12", "p_2 - p_5",
    "beta_13", "p_2 - p_8",
    "beta_23", "p_5 - p_8"
  )

  tmp5 <- tmp3 |>
  dplyr::select(para, para_name) |>
    bind_rows(tmp4) |>
    mutate(para_name = factor(para_name,
    levels = c(tmp3$para_name, tmp4$para_name) |> rev()))

  intervals_data <- mcmc_intervals_data(
    draws,
    # pars = c("mu_hat"),
    point_est = "median",
    prob = 0.5,
    prob_outer = 0.95,
    regex_pars = "alpha_[1-5]|beta_[1-9]") |>
    left_join(tmp5, by = c("parameter" = "para"))

  ggplot(intervals_data, aes(y = para_name)) +
    geom_vline(xintercept = 0, lty = 2, col = "grey60") +
    geom_linerange(aes(xmin = ll, xmax = hh)) +
    geom_linerange(aes(xmin = l, xmax = h), size = 2) +
    geom_point(aes(x = m), size = 3) +
    ylab("") +
    xlab("Posterior estimates") +
    my_theme()
}

coef_density_sd <- function(draws) {
# library(tidyverse)
# library(targets)
# library(ggridges)
# tar_load(fit_anova_noint_err_log_draws_anova_noint_err)
# draws <- fit_anova_noint_err_log_draws_anova_noint_err

# draws_data <- draws |>
#     janitor::clean_names()  |>
#     dplyr::select(tau_1, tau_2, sigma) |>
#     rename(Species = tau_1) |>
#     rename(Pressure = tau_2) |>
#     rename(Residuals = sigma) |>
#     pivot_longer(1:3)
  mcmc_areas(
    draws,
    pars = c("tau[1]", "tau[2]", "sigma"),
    prob = 0.5
    ) +
    # xlim(c(0, 2)) +
    scale_y_discrete(labels = c("Species", "Pressure", "Residuals"))
}

my_ggsave <- function(filename, plot, units = c("in", "cm",
        "mm", "px"), height = NA, width = NA, dpi = 300, ...) {

  ggsave(
    filename = paste0(filename, ".png"),
    plot = plot,
    height = height,
    width = width,
    units = units,
    dpi = dpi,
    ...
  )

  ggsave(
    filename = paste0(filename, ".pdf"),
    plot = plot,
    height = height,
    width = width,
    units = units,
    dpi = dpi,
    ...
  )

  paste0(filename, c(".png", ".pdf"))
}


logit <- function(z) {
  1 / (1 + exp(-z)) * 100
}

plot_logistic_sp <- function(draws, data_file) {
  d <- read_csv(data_file) |>
   filter(!is.na(count))

  x_mean <- mean(d$pressure)
  x_sd <- sd(d$pressure)
  xx_raw <- seq(0.02, 0.08, length = 80)
  xx <- (xx_raw - x_mean) / x_sd
  xx_mat <- cbind(1, xx, xx^2)

  sp <- d$species |> unique() |> as.factor() |> sort()

  gamma_ <- draws |>
    janitor::clean_names() |>
    dplyr::select(contains("gamma")) |>
    as.matrix()

  beta_ <- draws |>
    janitor::clean_names() |>
    dplyr::select(contains("beta"))

  tmp <- list()
  for (i in 1:5) {
    tmp[[i]] <- beta_ |>
    dplyr::select(matches(str_c("_", i, "$"))) |>
    as.matrix()
  }

  tmp2 <- tibble(data = map(tmp, \(x) x %*% t(xx_mat))) |>
    mutate(pred_m = map(data, \(x)apply(x, 2, quantile, 0.5))) |>
    mutate(pred_l = map(data, \(x)apply(x, 2, quantile, 0.25))) |>
    mutate(pred_h = map(data, \(x)apply(x, 2, quantile, 0.75))) |>
    mutate(pred_ll = map(data, \(x)apply(x, 2, quantile, 0.025))) |>
    mutate(pred_hh = map(data, \(x)apply(x, 2, quantile, 0.975))) |>
    dplyr::select(-data) |>
    mutate(species = sp) |>
    unnest(cols = c(pred_m, pred_l, pred_h, pred_ll, pred_hh))

  fig_data <- tmp2 |>
    mutate_if(is.numeric, logit) |>
    mutate(xx = rep(xx_raw, 5))

  my_col <- RColorBrewer::brewer.pal(3, "Set2")

# "#66C2A5"
# "#FC8D62"
# "#8DA0CB"

  tag_data <- tibble(
    species = sp,
    label = str_c(LETTERS[1:5], ": ", sp),
    x = 0.03,
    y = 95
  )

  ggplot() +
    geom_point(data = d, aes(x = pressure, y = count/total * 100, size = total), alpha = 0.2, col = my_col[3]) +
    geom_line(data = fig_data, aes(x = xx, y = pred_m), size = 0.5) +
    geom_ribbon(data = fig_data, aes(x = xx, ymin = pred_l, ymax = pred_h), alpha = 0.4) +
    geom_ribbon(data = fig_data, aes(x = xx, ymin = pred_ll, ymax = pred_hh), alpha = 0.4) +
    ylab("Proportion of vessel with silicon (%)") +
    xlab("Pressure (MPa)") +
    geom_text(data = tag_data,
      aes(x = x, y = y, label = label),
      size = 3,
      col = "grey10"
      ) +
    coord_cartesian(xlim = c(0.015, 0.085)) +
    scale_x_continuous(breaks = c(0.02, 0.05, 0.08)) +
    facet_grid(rows = vars(species), scale = "fixed") +
    my_theme() +
    theme(
      legend.position = "none",
        strip.background = element_blank(),
        strip.text.y = element_blank()
        # plot.margin = unit(c(1,1,1,1) , units = "lines" )
    )

}

