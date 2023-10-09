
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

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
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
    plot_annotation(tag_levels = "A") &
    theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12)
    )
}

sma_ks2 <- function(p1, p2) {
  p1 + p2 +
    plot_layout(nrow = 1, width = c(1, 2)) +
    plot_annotation(tag_levels = "A") &
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


logistic <- function(z) {
  1 / (1 + exp(-z)) * 100
}

plot_logistic_sp <- function(draws, data_file, quad = TRUE) {
  draws <- tar_read(quad_logistic_draws_hierarchical_logistic)
  d <- read_csv(here::here("data/cond_count.csv")) |>
  # d <- read_csv(data_file) |>
   filter(!is.na(count)) |>
   mutate(xylem_type = case_when(
    species == "TG" ~ "RP",
    species == "AP" ~ "L",
    TRUE ~ "DP"
   )) |>
   mutate(xylem_long = case_when(
      xylem_type == "DP"  ~ "Diffuse-porous tree",
      xylem_type == "RP"  ~ "Ring-porous tree",
      xylem_type == "Pa"  ~ "Palm",
      xylem_type == "L"  ~ "Liana"
    )) |>
    mutate(xylem_long_fct = factor(xylem_long,
    levels = c("Diffuse-porous tree", "Ring-porous tree", "Palm", "Liana"))) |>
  mutate(species = factor(species, levels = c("HB", "HH", "VM", "TG", "AP")))

  x_mean <- mean(d$pressure)
  x_sd <- sd(d$pressure)
  xx_raw <- seq(0.02, 0.08, length = 80)
  xx <- (xx_raw - x_mean) / x_sd
  xx_mat <- cbind(1, xx, xx^2)
  if (!quad) {
    xx_mat <- cbind(1, xx)
  }

  sp <- d$species |> unique()

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
    mutate(species = as.character(sp) |> sort()) |> # alphabetical order in stan
    unnest(cols = c(pred_m, pred_l, pred_h, pred_ll, pred_hh))

  fig_data <- tmp2 |>
    mutate_if(is.numeric, logistic) |>
    mutate(xx = rep(xx_raw, 5)) |>
    mutate(species = factor(species, levels = c("HB", "HH", "VM", "TG", "AP")))

  tag_data <- tibble(
    species = sp,
    label = str_c(LETTERS[1:5], ": ", sp),
    x = 0.03,
    y = 95
  )

  my_cols <- gg_color_hue(4)

  ggplot() +
    geom_point(data = d, aes(x = pressure, y = count/total * 100,
      size = total, col = xylem_long_fct),
    alpha = 0.3) +
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
    scale_color_manual(values = my_cols[-3]) +
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


coef_intervals_logistic <- function(draws) {
   sp <- c("AP",
      "HB",
      "HH",
      "TG",
      "VM") |> as.factor() |> sort() |> as.character()

  intervals_data <- mcmc_intervals_data(
    draws,
    regex_pars = "beta\\[|gamma\\[",
    prob = 0.5,
    prob_outer = 0.95) |>
    mutate(para = case_when(
      str_detect(parameter, "\\[1") ~ "Intercept",
      str_detect(parameter, "\\[2") ~ "Linear",
      str_detect(parameter, "\\[3") ~ "Quadratic",
    )) |>
    mutate(sp = case_when(
      str_detect(parameter, "gamma") ~ "All",
      str_detect(parameter, "1\\]") ~ sp[1],
      str_detect(parameter, "2\\]") ~ sp[2],
      str_detect(parameter, "3\\]") ~ sp[3],
      str_detect(parameter, "4\\]") ~ sp[4],
      str_detect(parameter, "5\\]") ~ sp[5],
    )) |>
    # mutate(para = str_c(sp, "-", para0)) |>
    mutate(para = factor(para,
      levels = c("Intercept", "Linear", "Quadratic") |> rev()))

    ggplot(intervals_data, aes(y = para)) +
      geom_vline(xintercept = 0, lty = 2, col = "grey60") +
      geom_linerange(aes(xmin = ll, xmax = hh)) +
      geom_linerange(aes(xmin = l, xmax = h), size = 2) +
      geom_point(aes(x = m), size = 2, position = position_nudge(y = 0.3)) +
      facet_wrap(~sp) +
      ylab("") +
      xlab("Posterior estimates") +
      my_theme()
}

line_pg_multi <- function(data, xylem_lab, k_range, s_002, s_0025, s_003, s_0035, s_004, s_005, s_006, s_007, s_008) {

# s_0025 <- withr::with_dir(rprojroot::find_root('_targets.R'),
#   targets::tar_read(fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.025))
# s_0035 <- withr::with_dir(rprojroot::find_root('_targets.R'),
#   targets::tar_read(fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.035))
# s_002 <- withr::with_dir(rprojroot::find_root('_targets.R'),
#   targets::tar_read(fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.02))
# s_003 <- withr::with_dir(rprojroot::find_root('_targets.R'),
#   targets::tar_read(fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.03))
# s_004 <- withr::with_dir(rprojroot::find_root('_targets.R'),
#   targets::tar_read(fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.04))
# s_005 <- withr::with_dir(rprojroot::find_root('_targets.R'),
#   targets::tar_read(fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.05))
# s_006 <- withr::with_dir(rprojroot::find_root('_targets.R'),
#   targets::tar_read(fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.06))
# s_007 <- withr::with_dir(rprojroot::find_root('_targets.R'),
#   targets::tar_read(fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.07))
# s_008 <- withr::with_dir(rprojroot::find_root('_targets.R'),
#   targets::tar_read(fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.08))

  # data <- tar_read(fd_k_traits_csv)
  xylem_lab2 <- xylem_lab |>
    select(sp_short, sp_short_chr, xylem_long_fct)
  d <- read_csv(data) |>
    filter(is.na(removed_k)) |>
    rename(sp_short_chr = sp_short) |>
    left_join(xylem_lab2, by = "sp_short_chr")

  s_002 <- s_002 |> mutate(p_g_lim = "0.02")
  s_0025 <- s_0025 |> mutate(p_g_lim = "0.025")
  s_003 <- s_003 |> mutate(p_g_lim = "0.03")
  s_0035 <- s_0035 |> mutate(p_g_lim = "0.035")
  s_004 <- s_004 |> mutate(p_g_lim = "0.04")
  s_005 <- s_005 |> mutate(p_g_lim = "0.05")
  s_006 <- s_006 |> mutate(p_g_lim = "0.06")
  s_007 <- s_007 |> mutate(p_g_lim = "0.07")
  s_008 <- s_008 |> mutate(p_g_lim = "0.08")

  data <- bind_rows(s_002, s_0025, s_003, s_0035, s_004, s_005, s_006, s_007, s_008) |>
    janitor::clean_names()

  # tar_load(k_range)
  # tar_load(xylem_lab)

  xylem_lab <- xylem_lab |>
    mutate(a_chr = str_c("alpha[1,", sp_num,"]")) |>
    mutate(b_chr = str_c("alpha[2,", sp_num,"]"))

  k_range2 <- left_join(k_range, xylem_lab, by = "species")

  log_a <- data |>
    filter(str_detect(variable, "alpha\\[1")) |>
    dplyr::select(variable, log_a = q50, p_g_lim)

  b <- data |>
    filter(str_detect(variable, "alpha\\[2")) |>
    dplyr::select(variable, b = q50, p_g_lim)
  sig <- data |>
    filter(variable == "sigma") |>
    dplyr::select(sigma = q50, p_g_lim)

  log_a2 <- log_a |>
    mutate(b = b$b)

  log_a3 <- full_join(log_a2, sig)

  nd <- left_join(k_range2, log_a3, by = c("a_chr" = "variable", "p_g_lim" = "p_g_lim")) |>
    group_by(p_g_lim, species) |>
    nest() |>
    ungroup()

  # nd$data[[1]]

  pred_data <- nd |>
    mutate(log_xx = map(data, \(x)seq(log(x$k_lwr), log(x$k_upr), length = 80))) |>
    mutate(log_pred = pmap(list(log_xx, data), \(log_xx, data) {data$log_a + data$b * log_xx + data$sigma^2 / 2})) |>
    unnest(cols = c(data, log_xx, log_pred))

  my_cols <- gg_color_hue(4)

  ggplot() +
    geom_point(data = d |>
      filter(xylem_type == "DP"), aes(x = k, y = fd), col = my_cols[1]) +
    geom_point(data = d |>
      filter(xylem_type == "RP"), aes(x = k, y = fd), col = my_cols[2]) +
    geom_point(data = d |>
      filter(xylem_type == "Pa"), aes(x = k, y = fd), col = my_cols[3]) +
    geom_point(data = d |>
      filter(xylem_type == "L"), aes(x = k, y = fd), col = my_cols[4]) +
    # geom_point(data = d, aes(x = k, y = fd, col = xylem_long_fct)) +
    # geom_line(data = pred_data, aes(x = exp(log_xx), y = exp(log_pred), group = p_g_lim)) +
    geom_line(data = pred_data, aes(x = exp(log_xx), y = exp(log_pred), col = as.numeric(p_g_lim), group = p_g_lim)) +
    facet_wrap(vars(sp_short), ncol = 4, scale = "free") +
    ylab(expression("Sap flux density "(g~m^{-2}~s^{-1}))) +
    xlab(expression("K "((Delta~T[max]-Delta~T)/Delta~T))) +
    scale_color_viridis_c(option = "D") +
    labs(col =expression(Maximum~italic(P[g])~(MPa~m^{-1}))) +
    my_theme() +
    theme(
      strip.text = element_text(face = "italic", size = 8),
      legend.position = c(0.85, 0.05)
      )

}

generate_k_range <- function(data) {
  d <- read_csv(data) |>
    filter(is.na(removed_k))
  # d_0015 <- d |> filter(p_g <= 0.015) |> mutate(p_g_lim = "0.015")
  d_002 <- d |> filter(p_g <= 0.02) |> mutate(p_g_lim = "0.02")
  d_0025 <- d |> filter(p_g <= 0.025) |> mutate(p_g_lim = "0.025")
  d_003 <- d |> filter(p_g <= 0.03) |> mutate(p_g_lim = "0.03")
  d_0035 <- d |> filter(p_g <= 0.035) |> mutate(p_g_lim = "0.035")
  d_004 <- d |> filter(p_g <= 0.04) |> mutate(p_g_lim = "0.04")
  d_005 <- d |> filter(p_g <= 0.05) |> mutate(p_g_lim = "0.05")
  d_006 <- d |> filter(p_g <= 0.06) |> mutate(p_g_lim = "0.06")
  d_007 <- d |> filter(p_g <= 0.07) |> mutate(p_g_lim = "0.07")
  d_008 <- d |> filter(p_g <= 0.08) |> mutate(p_g_lim = "0.08")

  k_data <- bind_rows(d_002, d_0025, d_003, d_0035, d_004, d_005, d_006, d_007, d_008)
  k_data2 <- k_data |>
    group_by(species, p_g_lim) |>
    nest() |>
    ungroup() |>
    mutate(k_lwr = map_dbl(data, \(x) min(x$k))) |>
    mutate(k_upr = map_dbl(data, \(x) max(x$k))) |>
    mutate(n = map_dbl(data, nrow)) |>
    arrange(species) |>
    dplyr::select(-data)

  k_data2 |>
    filter(str_detect(species, "assa"))

  # k_range |>
  #   filter(str_detect(species, "assa"))
  # k_range |>
  #   filter(p_g_lim == "0.02")
  # k_range |>
  #   filter(p_g_lim == "0.025")

  tmp <- k_data2$n
  tmp2 <- as.numeric(tmp)
  for (i in 2:nrow(k_data2)) {
   if (tmp[i-1] == tmp[i]) {
    tmp2[i] <- 1
   } else tmp2[i] <- 0
  }
  tmp2[1] <- 0

  k_data2 |>
    mutate(check = tmp2) |>
    filter(check == 0) |>
    filter(n >= 10)


}

coef_cor <- function() {
#  install.packages("ellipse")
  library(ellipse)
  library(targets)
  library(tidyverse)
  s_008 <- withr::with_dir(rprojroot::find_root('_targets.R'),
  targets::tar_read(fit_ab_summary_granier_without_traits_sap_all_clean_0.08)) |>
  janitor::clean_names() |>
  mutate(pressure = "0.08")

  d_008 <- withr::with_dir(rprojroot::find_root('_targets.R'),
  targets::tar_read(fit_ab_draws_granier_without_traits_sap_all_clean_0.08)) |>
  janitor::clean_names() |>
  mutate(pressure = "0.08")

 s_008 |>
    filter(str_detect(variable, "alpha"))

#    dplyr::select(starts_with("rho_l"))

 sigma_l <- d_008 |>
    dplyr::select(starts_with("sigma_l"))
 gamma <- d_008 |>
    dplyr::select(starts_with("gamma"))
 sigma_l <- sigma_l  |>
    mutate(s = pmap(list(sigma_l_1_1, sigma_l_2_1, sigma_l_2_2), \(s1, rho, s2) matrix(c(s1, rho, rho, s2), nrow = 2))) |>
    mutate(s_cor = map(s, cov2cor))

 gamma <- gamma |>
  mutate(mu = map2(gamma_1_1, gamma_2_1, \(x, y) c(x, y)))

  bind_cols(sigma_l, gamma) |>
  mutate(data = pmap(list(mu, s), \(mu, s) ellipse(s, scale = c(s[1, 1], s[2, 2]), centre = mu, level = 0.9)))  |>
  mutate(log_a = map(data, \(x) x[, 1])) |>
  mutate(b = map(data, \(x) x[, 2]))

 tmp <- bind_cols(sigma_l, gamma) |>
  mutate(data = pmap(list(mu, s), \(mu, s) ellipse(s, scale = c(s[1, 1], s[2, 2]), centre = c(0, 0), level = 0.9))) |>
  mutate(log_a = map(data, \(x) x[, 1])) |>
  mutate(b = map(data, \(x) x[, 2])) |>
  mutate(.id = 1:nrow(sigma_l)) |>
  dplyr::select(log_a, b, .id) |>
  unnest(cols = c(log_a, b))

ggplot(tmp, aes(x = b, y = log_a, group = .id)) +
  geom_path(alpha = 0.2) +
  xlab("x") +
  ylab("y") +
  theme_bw()

library(ggthemes)

ggplot(tmp, aes(x = b, y = log_a, group = .id)) +
  geom_path(alpha = 0.1, col = "orange")+
  xlab("x") +
  ylab("y") +
  theme_solarized(light = FALSE)

  tmp |>
    rename("hoge" = "[,2]")
  dim(tmp)

  m <- c(.5, -.5)
  sigma <- matrix(c(1,.5,.5,1), nrow=2)

  gamma <- s_008 |>
    filter(str_detect(variable, "gamma"))
  beta <- s_008 |>
    filter(str_detect(variable, "beta"))
  Sigma <- s_008 |>
    filter(str_detect(variable, "Sigma_l")) |>
    pull(q50)
  tau <- s_008 |>
    filter(str_detect(variable, "tau_l"))

  m <- gamma$q50
  sigma <- matrix(Sigma, nrow = 2)


  alpha_levels <- seq(0.5, 0.95, by = 0.05) ## or whatever you want
  alpha_levels <- 0.9
  names(alpha_levels) <- alpha_levels ## to get id column in result
  contour_data <- plyr::ldply(alpha_levels,
        ellipse,
        x = sigma,
        scale = c(sigma[1, 1], sigma[2, 2]),  ## needed for positional matching
        centre = m)
  ggplot(contour_data,aes(x,y,group=.id))+geom_path()

  ellipse(sigma)


  s_008 <- tar_read(fit_ab_summary_granier_without_traits_sap_all_clean_0.08)

  s_008 |>
    filter(str_detect(variable, "rho"))

  s_008 |>
    filter(str_detect(variable, "Sigma"))
  s_008 |>
    filter(str_detect(variable, "Sigma"))

  s_008 |>
    filter(str_detect(variable, "tau"))

}

line_pool_multi <- function(d, xylem_lab, s_008, s2_008) {
  # s_008 <- tar_read(fit_ab_summary_granier_without_traits_sap_all_clean_0.08)
  # s2_008 <- tar_read(fit_ab_summary_granier_without_traits2_sap_all_clean_0.08)

  xylem_lab2 <- xylem_lab |>
    select(sp_short, sp_short_chr, xylem_long_fct)

  d <- read_csv(d) |>
    filter(is.na(removed_k)) |>
    rename(sp_short_chr = sp_short) |>
    left_join(xylem_lab2, by = "sp_short_chr")

  log_a <- s_008 |>
    filter(str_detect(variable, "alpha\\[1")) |>
    pull(q50)
  b <- s_008 |>
    filter(str_detect(variable, "alpha\\[2")) |>
    pull(q50)
  sig <- s_008 |>
    filter(variable == "sigma") |>
    pull(q50)

  log_a_pool <- s2_008 |>
    filter(str_detect(variable, "alpha\\[1")) |>
    pull(q50)
  b_pool <- s2_008 |>
    filter(str_detect(variable, "alpha\\[2")) |>
    pull(q50)
  sig_pool <- s2_008 |>
    filter(variable == "sigma") |>
    pull(q50)

  nd <- d |>
    group_by(species, sp_short, xylem_type) |>
    nest() |>
    ungroup() |>
    arrange(species) |>
    mutate(log_xx = map(data, \(x)seq(min(log(x$k)), max(log(x$k)), length = 80))) |>
    mutate(log_a = log_a) |>
    mutate(b = b) |>
    mutate(sig = sig) |>
    mutate(log_a_pool = log_a_pool) |>
    mutate(b_pool = b_pool) |>
    mutate(sig_pool = sig_pool) |>
    mutate(log_pred = pmap(list(log_xx, log_a, b, sig), \(log_xx, log_a, b, sig) log_a + b * log_xx + sig^2 / 2)) |>
    mutate(log_pred_pool = pmap(list(log_xx, log_a_pool, b_pool, sig_pool), \(log_xx, log_a, b, sig) log_a + b * log_xx + sig^2 / 2))

  pred_data <- nd |>
    dplyr::select(xylem_type, species, sp_short, log_xx, log_pred, log_pred_pool) |>
    unnest(c(log_xx, log_pred, log_pred_pool))

  d |>
    ggplot() +
    geom_point(aes(x = k, y = fd, col = xylem_long_fct)) +
    geom_line(data = pred_data, aes(x = exp(log_xx), y = exp(log_pred))) +
    geom_line(data = pred_data, aes(x = exp(log_xx), y = exp(log_pred_pool)), lty = 2) +
    facet_wrap(vars(sp_short), ncol = 4, scale = "free") +
    ylab(expression("Sap flux density "(g~m^{-2}~s^{-1}))) +
    xlab(expression("K "((Delta~T[max]-Delta~T)/Delta~T))) +
    labs(col = "") +
    my_theme() +
    theme(
      strip.text = element_text(face = "italic", size = 8),
      legend.position = c(0.85, 0.05)
      )

}


  # library(tidyverse)
  # library(here)
  # library(ggridges)
  # library(scales)
  # library(targets)
  # d <- read_csv(here("data/fd_k_traits.csv")) |>
  #   filter(is.na(removed_k))
  #  draws <- withr::with_dir(rprojroot::find_root('_targets.R'),
  #    targets::tar_read(fit_ab_draws_granier_without_traits_full_segments_sap_all_clean_0.08)) |>
  #    janitor::clean_names()
  # tar_load(xylem_lab)
coef_density <- function(xylem_lab, draws) {
  draws <- draws |>
    janitor::clean_names()

  n_iter <- nrow(draws)

  tmp <- xylem_lab |>
    dplyr::select(xylem_fct, xylem_long_fct) |> unique()

  xy_data_a <- draws |>
    dplyr::select(starts_with("beta_1")) |>
    pivot_longer(1:4)  |>
    arrange(name) |>
    mutate(xylem_fct = rep(unique(xylem_lab$xylem_fct) |> sort(), each = n_iter)) |>
    left_join(tmp, by = "xylem_fct") |>
    mutate(sp = xylem_long_fct) |>
    mutate(sp_short = xylem_long_fct) |>
    dplyr::select(sp, sp_short, xylem = xylem_long_fct, value)

  sp_data_a <- draws |>
    dplyr::select(starts_with("alpha_1")) |>
    pivot_longer(1:31) |>
    slice(str_order(name, numeric = TRUE)) |>
    mutate(sp = rep(str_sort(xylem_lab$species), each = n_iter)) |>
    left_join(xylem_lab, by = c("sp" = "species")) |>
    dplyr::select(sp, sp_short, xylem = xylem_long_fct, value)

  data_a <- bind_rows(xy_data_a, sp_data_a) |>
    mutate(sp_short = factor(sp_short,
     c(
      xy_data_a$xylem |> unique() |> sort() |> as.character(),
      sp_data_a$sp_short |> unique() |> str_sort()))) |>
    mutate(para = "Coefficient a") |>
    mutate(group = ifelse(str_detect(sp_short, "^*\\."), "Species", "Xylem")) |>
    mutate(group = factor(group, levels = c("Xylem", "Species")))

  xy_data_b <- draws |>
    dplyr::select(starts_with("beta_2")) |>
    pivot_longer(1:4)  |>
    arrange(name) |>
    mutate(xylem_fct = rep(unique(xylem_lab$xylem_fct) |> sort(), each = n_iter)) |>
    left_join(tmp, by = "xylem_fct") |>
    mutate(xylem_fct = rep(unique(xylem_lab$xylem_fct) |> sort(), each = n_iter)) |>
    mutate(sp = xylem_long_fct) |>
    mutate(sp_short = xylem_long_fct) |>
    dplyr::select(sp, sp_short, xylem = xylem_long_fct, value)

  sp_data_b <- draws |>
    dplyr::select(starts_with("alpha_2")) |>
    pivot_longer(1:31)  |>
    slice(str_order(name, numeric = TRUE)) |>
    mutate(sp = rep(sort(xylem_lab$species), each = n_iter)) |>
    left_join(xylem_lab, by = c("sp" = "species")) |>
    dplyr::select(sp, sp_short, xylem = xylem_long_fct, value)

  data_b <- bind_rows(xy_data_b, sp_data_b) |>
    mutate(sp_short = factor(sp_short,
     c(
      xy_data_a$xylem |> unique() |> sort() |> as.character(),
      sp_data_b$sp_short |> unique() |> str_sort()))) |>
    mutate(para = "Coefficient b") |>
    mutate(group = ifelse(str_detect(sp_short, "^*\\."), "xylem", "sp"))

  data <- bind_rows(data_a, data_b)

  tmp <- data_a  |>
    dplyr::select(sp_short, xylem) |>
    unique() |>
    group_by(xylem) |>
    nest() |>
    arrange(xylem) |>
    unnest(data) |>
    ungroup() |>
    mutate(sp_short_chr = as.character(sp_short)) |>
    pull(sp_short_chr)
  # tmp2 <- tmp[str_detect(tmp, "Palm|tree|Liana")]

  tmp2 <- xylem_lab |>
    pull(xylem_long_fct) |>
    unique() |>
    levels()

  tmp3 <- tmp[!str_detect(tmp, "Palm|tree|Liana")]

  data_a <- data_a |>
    mutate(sp_short2 = factor(sp_short, levels = c(tmp2, tmp3)))
  data_b <- data_b |>
    mutate(sp_short2 = factor(sp_short, levels = c(tmp2, tmp3)))

  y_lab_raw <- data_a$sp_short2  |> levels()
  y_lab <- as.expression(y_lab_raw)

  names(y_lab) <- y_lab_raw
  # R 4.2.2 can handle comparisons of expressions
  y_lab <- case_when(
    y_lab == "A. fraxinifolius" ~ expression(italic("A. fraxinifolius")),
    y_lab == "D. alatus" ~ expression(italic("D. alatus")),
    y_lab == "D. tonkinensis" ~ expression(italic("D. tonkinensis")),
    y_lab == "D. turbinatus" ~ expression(italic("D. turbinatus")),
    y_lab == "H. cordifolia" ~ expression(italic("H. cordifolia")),
    y_lab == "H. brasiliensis" ~ expression(italic("H. brasiliensis")),
    y_lab == "H. hongaychsis" ~ expression(italic("H. hongaychsis")),
    y_lab == "J. mimosifolia" ~ expression(italic("J. mimosifolia")),
    y_lab == "L. coromandelica" ~ expression(italic("L. coromandelica")),
    y_lab == "L. leucocephala" ~ expression(italic("L. leucocephala")),
    y_lab == "M. ferrea" ~ expression(italic("M. ferrea")),
    y_lab == "P. chinensis" ~ expression(italic("P. chinensis")),
    y_lab == "P. cerasoides" ~ expression(italic("P. cerasoides")),
    y_lab == "P. tomentosa" ~ expression(italic("P. tomentosa")),
    y_lab == "S. assamica" ~ expression(italic("S. assamica")),
    y_lab == "T. franchetii" ~ expression(italic("T. franchetii")),
    y_lab == "V. mangachapoi" ~ expression(italic("V. mangachapoi")),
    y_lab == "V. xishuangbannaensis" ~ expression(italic("V. xishuangbannaensis")),
    y_lab == "A. pennata" ~ expression(italic("A. pennata")),
    y_lab == "B. tenuiflor" ~ expression(italic("B. tenuiflor")),
    y_lab == "G. montanum" ~ expression(italic("G. montanum")),
    y_lab == "M. pachycarpa" ~ expression(italic("M. pachycarpa")),
    y_lab == "V. calyculata" ~ expression(italic("V. calyculata")),
    y_lab == "A. alexandrae" ~ expression(italic("A. alexandrae")),
    y_lab == "A. catechu" ~ expression(italic("A. catechu")),
    y_lab == "A. triandra" ~ expression(italic("A. triandra")),
    y_lab == "C. mitis" ~ expression(italic("C. mitis")),
    y_lab == "C. lutescens" ~ expression(italic("C. lutescens")),
    y_lab == "D. lutescen" ~ expression(italic("D. lutescen")),
    y_lab == "M. toosendan" ~ expression(italic("M. toosendan")),
    y_lab == "T. grandis" ~ expression(italic("T. grandis")),
    TRUE ~ y_lab
  )
  names(y_lab) <- y_lab_raw

  p1 <- ggplot(data_a, aes(x = exp(value), y = sp_short2, fill = xylem))  +
    facet_grid(group ~ ., scales = "free", space = "free") +
    scale_x_log10(
        breaks = c(10^2, 10^3, 10^4),
        labels = trans_format("log10", math_format(10^.x))) +
    theme_bw() +
    theme(
      legend.position = "none",
      # axis.text.y = element_text(face = "italic", size = 8),
      strip.background = element_blank(),
      strip.text = element_blank()
      ) +
    geom_vline(xintercept = 119, lty = 1, col = "grey40") +
    scale_y_discrete(labels = y_lab, limits = rev) +
    geom_density_ridges(col = "grey92") +
    xlab(expression(Coefficient~italic(a))) +
    ylab("")

  # my_col <- RColorBrewer::brewer.pal(4, "PuOr")
  p2 <- ggplot(data_b, aes(x = value, y = sp_short2, fill = xylem))  +
    facet_grid(group ~ ., scales = "free", space = "free") +
    theme_bw() +
    theme(
     legend.position = "none",
      axis.text.y = element_blank(),
      #axis.text.y = element_text(face = "italic", size = 8)
      strip.background = element_blank(),
      strip.text = element_blank()
      )  +
    geom_vline(xintercept = 1.23, lty = 1, col = "grey40") +
    geom_density_ridges(col = "grey92") +
    scale_y_discrete(limits = rev) +
    # scale_fill_manual(values = my_col) +
    xlab(expression(Coefficient~italic(b))) +
    ylab("")

  p1 + p2 +
    plot_annotation(tag_levels = "A")

}


generate_xylem_lab <- function(data, removed_k = TRUE) {
  d <- read_csv(data)
  if (removed_k) {
    d <- d |>  filter(is.na(removed_k))
  }
  xylem_lab <- d |>
    dplyr::select(xylem_type, species, sp_short) |>
    unique() |>
    mutate(xylem_fct = as.factor(xylem_type)) |>
    arrange(xylem_fct) |>
    mutate(xylem_long = case_when(
      xylem_type == "DP"  ~ "Diffuse-porous tree",
      xylem_type == "RP"  ~ "Ring-porous tree",
      xylem_type == "Pa"  ~ "Palm",
      xylem_type == "L"  ~ "Liana"
    )) |>
    mutate(xylem_long_fct = factor(xylem_long,
    levels = c("Diffuse-porous tree", "Ring-porous tree", "Palm", "Liana"))) |>
    mutate(sp_fct = as.factor(species)) |>
    mutate(sp_num = as.numeric(sp_fct)) |>
    mutate(sp_num1 = str_c("1_", sp_num))  |>
    mutate(sp_num2 = str_c("2_", sp_num))  |>
    arrange(xylem_long_fct, species) |>
    mutate(sp_short_chr = sp_short) |>
    mutate(sp_short = factor(sp_short_chr, levels = sp_short_chr))
  xylem_lab
}


# s_002 <- tar_read(fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.02)
# s_0025 <- tar_read(fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.025)
# s_003 <- tar_read(fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.03)
# s_0035 <- tar_read(fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.035)
# s_004 <- tar_read(fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.04)
# s_005 <- tar_read(fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.05)
# s_006 <- tar_read(fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.06)
# s_007 <- tar_read(fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.07)
# s_008 <- tar_read(fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.08)
# tar_load(xylem_lab)
# tar_load(k_range)
ab_pg_ribbon <- function(xylem_lab, k_range, s_002, s_0025, s_003, s_0035, s_004, s_005, s_006, s_007, s_008, coef_a = TRUE) {

  s_002 <- s_002 |> mutate(p_g_lim = "0.02")
  s_0025 <- s_0025 |> mutate(p_g_lim = "0.025")
  s_003 <- s_003 |> mutate(p_g_lim = "0.03")
  s_0035 <- s_0035 |> mutate(p_g_lim = "0.035")
  s_004 <- s_004 |> mutate(p_g_lim = "0.04")
  s_005 <- s_005 |> mutate(p_g_lim = "0.05")
  s_006 <- s_006 |> mutate(p_g_lim = "0.06")
  s_007 <- s_007 |> mutate(p_g_lim = "0.07")
  s_008 <- s_008 |> mutate(p_g_lim = "0.08")
  data <- bind_rows(s_002, s_0025, s_003, s_0035, s_004, s_005, s_006, s_007, s_008) |>
    janitor::clean_names()

  # xylem_lab$species |> unique()
  if (coef_a) {
    data <- data |>
      filter(str_detect(variable, "alpha\\[1"))
    fig_data <- xylem_lab |>
      mutate(variable = str_c("alpha[1,", sp_num,"]")) |>
      dplyr::select(sp_short, species, variable, xylem_long_fct) |>
      full_join(data, by = "variable") |>
      mutate(q2_5 = exp(q2_5)) |>
      mutate(q25 = exp(q25)) |>
      mutate(q50 = exp(q50)) |>
      mutate(q75 = exp(q75)) |>
      mutate(q97_5 = exp(q97_5))
  } else {
    data <- data |>
      filter(str_detect(variable, "alpha\\[2"))
    fig_data <- xylem_lab |>
      mutate(variable = str_c("alpha[2,", sp_num,"]")) |>
      dplyr::select(sp_short, species, variable, xylem_long_fct) |>
      full_join(data, by = "variable")
  }

  fig_data <- right_join(fig_data, k_range, by = c("species", "p_g_lim")) |>
    mutate(p_g_lim = as.numeric(p_g_lim))

  p <- ggplot(fig_data, aes(x = p_g_lim, y = q50, group = sp_short,
    fill = xylem_long_fct)) +
      facet_wrap(~sp_short, ncol = 4, scale = "free_y") +
      geom_ribbon(aes(ymin = q2_5, ymax = q97_5), alpha = 0.4) +
      geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.6) +
      geom_point() +
      geom_line() +
      scale_fill_discrete(name = "") +
      xlab(expression(Maximum~italic(P[g])~(MPa~m^{-1}))) +
      my_theme() +
      theme(
        strip.text = element_text(face = "italic", size = 8),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        legend.position = c(0.85, 0.05))

  if (coef_a) {
    p +
      scale_y_log10() +
      ylab(expression(Coefficient~italic(a)))

  } else {
    p +
      ylab(expression(Coefficient~italic(b)))
  }

}

#' @title predicitons for beta
#' @para mat_pred dataframe for MCMC draws of mean predictions (n.mcmc x 80)
#' @para df_pred dataframe for trait (80 x 2)
#' @ref https://www.tjmahr.com/visualizing-uncertainty-rstanarm/
tidy_predictions <- function(
  mat_pred,
  df_data,
  obs_name = "observation",
  prob_lwr = .025,
  prob_upr = .975
) {
  # Get dataframe with one row per fitted value per posterior sample
  df_pred <- mat_pred |>
    as_tibble() |>
    setNames(seq_len(ncol(mat_pred))) |>
    tibble::rownames_to_column("posterior_sample") |>
    tidyr::pivot_longer(
      cols = c(-posterior_sample),
      names_to = obs_name,
      values_to = "fitted"
    )

  # Helps with joining later
  class(df_pred[[obs_name]]) <- class(df_data[[obs_name]])

  # Summarise prediction interval for each observation
  df_pred |>
    group_by(.data[[obs_name]]) |>
    summarise(
      mean = mean(fitted),
      lower = quantile(fitted, prob_lwr),
      upper = quantile(fitted, prob_upr)
    ) |>
    left_join(df_data, by = obs_name)

}

#' @title generate beta plot data (takes time)
generate_beta_list <- function(fit_beta, fit_gamma, stan_data, draws, x, y, x_lab, y_lab) {
  tmp_beta <- fit_beta |>
    filter(pred_name == paste(y))
  gamma_slope <- fit_gamma |>
    filter(pred_name == paste(y)) |>
    filter(trait_name == paste(x))
  gamma_int <- fit_gamma |>
    filter(pred_name == paste(y)) |>
    filter(trait_name == "intercept")
  trait <- stan_data$u[paste(x), ]

  tmp_para <- gamma_slope |>
    pull(para) |>
    str_split_fixed("_", 3)

  beta_k <- tmp_beta$para |> str_split_fixed("_", 3)
  k <- beta_k[, 2] |> unique()

  y_lab_parse <- paste0("expression(", y_lab ,"(beta[paste(", k, ",',',", "j)]))")

  beta_trait <- bind_cols(tmp_beta, trait = trait)

  x_steps <- seq(min(beta_trait$trait), max(beta_trait$trait), length = 80)
  new_data <- tibble(
    observation = seq_along(x_steps) |> as.character(),
    x_lt = x_steps)

#  tic()
  tmp <- draws |>
    janitor::clean_names() |>
    dplyr::select(
      gamma_int |> pull(para),
      gamma_slope |> pull(para))
#  toc()

  colnames(tmp) <- c("gamma_int_draw", "gamma_slope_draw")

  tmp2 <- tmp |>
    mutate(rep = paste0("rep", 1:nrow(tmp))) |>
    nest(data = c(gamma_int_draw, gamma_slope_draw)) |>
    mutate(x = list(new_data$x_lt)) |>
    mutate(y = map2(x, data, \(x, data) {data$gamma_int_draw + data$gamma_slope_draw * x}))

  pred_lin <- matrix(unlist(tmp2$y), ncol = nrow(draws), byrow = FALSE) |> t()
  # dim(pred_lin)
  df_pred_lin <- tidy_predictions(pred_lin, new_data)

  list(
    beta_trait = beta_trait,
    df_pred_lin = df_pred_lin,
    x_lab  = x_lab,
    y_lab  = y_lab_parse
    )
}

ab_comp_points <- function(pool_csv, seg_csv, xylem_lab) {
  dp <- read_csv(pool_csv)
  dp2 <- dp |>
    filter(str_detect(variable, "alpha")) |>
    mutate(species = str_split_fixed(para, "_a_|_b_", 2)[, 2]) |>
    mutate(coef = ifelse(str_detect(para, "_a_"), "a", "b")) |>
    dplyr::select(species, coef,
      q50_pool = q50, q2.5_pool = q2.5 , q97.5_pool = q97.5)

  ds <- read_csv(seg_csv)
  ds2 <- ds |>
    filter(str_detect(variable, "alpha")) |>
    mutate(species = str_split_fixed(para, "_a_|_b_", 2)[, 2]) |>
    mutate(coef = ifelse(str_detect(para, "_a_"), "a", "b")) |>
    dplyr::select(species, coef,
      q50_seg = q50, q2.5_seg = q2.5 , q97.5_seg = q97.5)

  d <- full_join(ds2, dp2) |>
    full_join(xylem_lab)

  da <- d |>
    filter(coef == "a") |>
    mutate_if(is.numeric, exp)
  min_a <- min(da$q2.5_seg, da$q2.5_pool)
  max_a <- max(da$q97.5_seg, da$q97.5_pool)

  p1 <- da |>
    ggplot(aes(x = q50_pool, y = q50_seg, col = xylem_long_fct)) +
    geom_abline(slope = 1, intercept = 0, lty = 2) +
    geom_point() +
    geom_errorbar(aes(ymin = q2.5_seg, ymax = q97.5_seg), alpha = 0.5) +
    geom_errorbar(aes(xmin = q2.5_pool, xmax = q97.5_pool), alpha = 0.5) +
    scale_x_log10() +
    scale_y_log10() +
    xlab("a - traditional fitting") +
    ylab("a - multilevel model") +
    my_theme() +
    labs(col = "") +
    theme(
      # legend.position = c(0.75, 0.25),
      legend.position = "none"
      # legend.text = element_text(size = 7)
      )

  db <- d |>
    filter(coef == "b")
  min_b <- min(db$q2.5_seg, db$q2.5_pool)
  max_b <- max(db$q97.5_seg, db$q97.5_pool)

  p2 <- db |>
    ggplot(aes(x = q50_pool, y = q50_seg, col = xylem_long_fct)) +
    geom_abline(slope = 1, intercept = 0, lty = 2) +
    geom_point() +
    geom_errorbar(aes(ymin = q2.5_seg, ymax = q97.5_seg), alpha = 0.5) +
    geom_errorbar(aes(xmin = q2.5_pool, xmax = q97.5_pool), alpha = 0.5) +
    xlab("b - traditional fitting") +
    ylab("b - multilevel model") +
    coord_cartesian(xlim = c(min_b, max_b), ylim = c(min_b, max_b)) +
    labs(col = "") +
    my_theme() +
    theme(
      # legend.position = "none"
      legend.position = c(0.7, 0.25),
      legend.text = element_text(size = 7)
      )

  p1 + p2 +
      plot_annotation(tag_levels = "A")
}

ab_comp_four_models_points <- function(summary12, summary3, summary4, xylem_lab, rm_dip = TRUE) {

# Extract percentiles for model 1 and 2
  extract_percentiles <- function(fits, var_name) {
    sapply(fits, function(x) {
      x$summary %>%
        janitor::clean_names() %>%
        filter(variable == var_name) %>%
        dplyr::select(q2_5, q50, q97_5) %>%
        unlist()
    }) |> t() |> as_tibble()
  }

# Extract percentiles for model 3 and 4
  extract_alpha_percentiles <- function(data, index) {
    data %>%
      janitor::clean_names() %>%
      filter(str_detect(variable, paste0("alpha\\[", index))) %>%
      dplyr::select(q2_5, q50, q97_5)
  }

# Model 1 and 2
  # d <- targets::tar_read(fit_ab_each_sap_sp_clean_0.08)
  d <- summary12
  log_a <- extract_percentiles(d$fit_pool, "log_a") |> exp()
  b <- extract_percentiles(d$fit_pool, "b")
  log_a2 <- extract_percentiles(d$fit_segments, "log_a") |> exp()
  b2 <- extract_percentiles(d$fit_segments, "b")

# Model 3
  # d3 <- targets::tar_read(fit_ab_summary_granier_without_traits_full_pool_sap_all_clean_0.08)
  d3 <- summary3
  log_a3 <- extract_alpha_percentiles(d3, 1) |> exp()
  b3 <- extract_alpha_percentiles(d3, 2)

# Model 4
  # d4 <- targets::tar_read(fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.08)
  d4 <- summary4
  log_a4 <- extract_alpha_percentiles(d4, 1) |> exp()
  b4 <- extract_alpha_percentiles(d4, 2)

  d1 <- tibble(species = d$species) |>
    bind_cols(
      log_a |> rename(a1_q2_5 = q2_5, a1_q50 = q50, a1_q97_5 = q97_5),
      b |> rename(b1_q2_5 = q2_5, b1_q50 = q50, b1_q97_5 = q97_5)
    ) |> full_join(xylem_lab)
  d2 <- bind_cols(
      d1,
      log_a2 |> rename(a2_q2_5 = q2_5, a2_q50 = q50, a2_q97_5 = q97_5),
      b2 |> rename(b2_q2_5 = q2_5, b2_q50 = q50, b2_q97_5 = q97_5)
  )
  d3 <- bind_cols(
      d1,
      log_a3 |> rename(a3_q2_5 = q2_5, a3_q50 = q50, a3_q97_5 = q97_5),
      b3 |> rename(b3_q2_5 = q2_5, b3_q50 = q50, b3_q97_5 = q97_5)
  )
  d4 <- bind_cols(
      d1,
      log_a4 |> rename(a4_q2_5 = q2_5, a4_q50 = q50, a4_q97_5 = q97_5),
      b4 |> rename(b4_q2_5 = q2_5, b4_q50 = q50, b4_q97_5 = q97_5)
  )


# Common plotting function
  ab_each_tmp <- function(data, x_var, y_var, x_label, y_label, use_log_scale = TRUE, rm_dip = TRUE) {

    x_q50 <- sym(paste0(x_var, "_q50"))
    y_q50 <- sym(paste0(y_var, "_q50"))
    x_q2_5 <- sym(paste0(x_var, "_q2_5"))
    y_q2_5 <- sym(paste0(y_var, "_q2_5"))
    x_q97_5 <- sym(paste0(x_var, "_q97_5"))
    y_q97_5 <- sym(paste0(y_var, "_q97_5"))

    if (rm_dip) {
    data <- data |>
      filter(species != "Dipterocarpus tonkinensis")
    }

    p <- data |>
      ggplot(aes(!!x_q50, !!y_q50, col = xylem_long_fct)) +
      geom_abline(slope = 1, intercept = 0, lty = 2) +
      geom_point() +
      geom_errorbar(aes(ymin = !!y_q2_5, ymax = !!y_q97_5), alpha = 0.5) +
      geom_errorbar(aes(xmin = !!x_q2_5, xmax = !!x_q97_5), alpha = 0.5) +
      xlab(x_label) +
      ylab(y_label) +
      my_theme() +
      labs(col = "") +
      theme(legend.position = "none")

    if (use_log_scale) {
      p <- p + scale_x_log10() + scale_y_log10()
    } else {
      if (rm_dip) {
        p <- p + coord_cartesian(xlim = c(0.25, 1.8), ylim = c(0.25, 1.8))
      } else {
        p <- p + coord_cartesian(xlim = c(0.25, 3), ylim = c(0.25, 3))
      }
    }
    return(p)
  }

# Individual plots
  p_a2 <- ab_each_tmp(d2, "a1", "a2", "a - traditional fitting", "a - model 2", rm_dip = rm_dip)
  p_b2 <- ab_each_tmp(d2, "b1", "b2", "b - traditional fitting", "b - model 2", use_log_scale = FALSE, rm_dip = rm_dip)
  p_a3 <- ab_each_tmp(d3, "a1", "a3", "a - traditional fitting", "a - model 3", rm_dip = rm_dip)
  p_b3 <- ab_each_tmp(d3, "b1", "b3", "b - traditional fitting", "b - model 3", use_log_scale = FALSE, rm_dip = rm_dip)
  p_a4 <- ab_each_tmp(d4, "a1", "a4", "a - traditional fitting", "a - model 4", rm_dip = rm_dip)
  p_b4 <- ab_each_tmp(d4, "b1", "b4", "b - traditional fitting", "b - model 4", use_log_scale = FALSE, rm_dip = rm_dip)

# Combined plot
  combined_plot <- p_a2 + p_b2 +
    p_a3 + p_b3 +
    p_a4 + p_b4 +
    plot_layout(ncol = 2) +
    plot_annotation(tag_levels = "A") &
      theme(
        axis.ticks.length = unit(-0.1, "cm"),
        axis.text.x = element_text(size = 8, margin = margin(t = 0.5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 8, margin = margin(t = 0, r = 0.5, b = 0, l = 0)),
        # axis.title.x = element_text(size = 7.5, margin = margin(t = 0, r = 0, b = 0, l = 0)),
        # axis.title.y = element_text(size = 7.5, margin = margin(t = 0, r = 0, b = 0, l = 0)),
        # axis.text = element_text(size = 6),
        plot.margin=unit(c(0.1, 0.1, 0.1, 0.1), "lines")
        # plot.tag = element_text(size = 8)
    )

  combined_plot
}
# targets::tar_make(ab_points_four_models_plot)
# 800m2
tr_scaled_bars <- function(ab_uncertainty_full_df) {
  d <- ab_uncertainty_full_df |>
    filter(!is.na(tree)) |>
    mutate(s_total_granier = ifelse(pg == 0.02, s_total_granier, NA)) |>
    rename(s_total_granier_m = s_total_granier) |>
    group_by(tree, pg) |>
    summarize(
      across(
        .cols = where(is.numeric),
        .fns = sum
      ), .groups = "drop")

  d2 <- d |>
    pivot_longer(3:13,
      names_to = "variable",
      values_to = "value") |>
    mutate(model = str_split_fixed(variable, "_", 4)[,3]) |>
    mutate(quantile = str_split_fixed(variable, "_", 4)[,4])

  d3 <- d2 |>
    group_by(pg, model, quantile) |>
    summarize(tr = sum(value) / 800 * 1000) |>
    ungroup() |>
    mutate(pg = ifelse(model == "granier", "NA", pg)) |>
    mutate(pg = factor(pg, levels = c(
      "NA",
      "0.02", "0.025", "0.03", "0.035", "0.04", "0.05", "0.06", "0.07", "0.08"
    ))) |>
    filter(!is.na(tr))

  d4 <- d3 |>
    pivot_wider(names_from = quantile, values_from = tr) #|>
    # complete(pg, nesting(model), fill = list(tr = 0))

  ggplot(d4, aes(x = pg, y = m, fill = model, group = model)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = ll, ymax = hh), position = position_dodge(.9), width = 0.2) +
    scale_fill_viridis_d() +
    ylab(expression(paste("Annual sap flux (Kg m"^-2~" year"^-1~")"))) +
    xlab(expression(Applied~italic(P[g])~(MPa~m^{-1}))) +
    theme_bw()
}


generate_tr_scaled_bars_data <- function(ab_uncertainty_full_df, each = FALSE) {
  d <- ab_uncertainty_full_df |>
    filter(!is.na(tree)) |>
    mutate(s_total_granier = ifelse(pg == 0.02, s_total_granier, NA)) |>
    rename(s_total_granier_m = s_total_granier) |>
    group_by(tree, pg) |>
    summarize(
      across(
        .cols = where(is.numeric),
        .fns = sum
      ), .groups = "drop")

  if (each) tmp <- 12 else tmp <- 13
  d2 <- d |>
    pivot_longer(3:tmp,
      names_to = "variable",
      values_to = "value") |>
    mutate(model = str_split_fixed(variable, "_", 4)[,3]) |>
    mutate(quantile = str_split_fixed(variable, "_", 4)[,4])

  d3 <- d2 |>
    group_by(pg, model, quantile) |>
    summarize(tr = sum(value) / 800 * 1000) |>
    ungroup() |>
    mutate(pg = ifelse(model == "granier", "NA", pg)) |>
    mutate(pg = factor(pg, levels = c(
      "NA",
      "0.02", "0.025", "0.03", "0.035", "0.04", "0.05", "0.06", "0.07", "0.08"
    ))) |>
    filter(!is.na(tr))

  d3 |>
    pivot_wider(names_from = quantile, values_from = tr) #|>
}

tr_scaled_bars2 <- function(ab_uncertainty_full_df, ab_uncertainty_full_each_df) {
  tmp1 <- generate_tr_scaled_bars_data(ab_uncertainty_full_each_df, each = TRUE) |>
    mutate(model = paste("full", model, sep = "_"))
  tmp2 <- generate_tr_scaled_bars_data(ab_uncertainty_full_df) |>
  mutate(model = paste("sep", model, sep = "_"))
  tmp3 <- bind_rows(tmp1, tmp2)
  ggplot(tmp3, aes(x = pg, y = m, fill = model, group = model)) +
   geom_bar(stat = "identity", position = "dodge") +
   geom_errorbar(aes(ymin = ll, ymax = hh), position = position_dodge(.9), width = 0.2) +
   scale_y_continuous(breaks = c(0, 250, 500, 750, 1000)) +  # This line sets the y-axis breaks
   scale_fill_viridis_d() +
   geom_hline(yintercept = 750, lty = 2) +
   ylab(expression(paste("Annual sap flux (Kg m"^-2~" year"^-1~")"))) +
   xlab(expression(Applied~italic(P[g])~(MPa~m^{-1}))) +
   theme_bw()
}

dbh_points <- function(dbh_imp_df, girth_increment_csv) {
  girth <- read_csv(girth_increment_csv) |>
    janitor::clean_names() |>
    mutate(date = lubridate::dmy(date))

  tmp <- dbh_imp_df |>
    filter(date %in% c(as.Date("2014-12-24"), girth$date))

  tmp2 <- tmp |> filter(date == max(tmp$date))

  ggplot(dbh_imp_df, aes(date, dbh, group = tree)) +
    geom_line(lty = 2, width = 1) +
    geom_point(data = tmp) +
    geom_text_repel(data = tmp2,
      aes(label = tree),
      size = 3.5,
      nudge_y = 1, nudge_x = 110,
      arrow = grid::arrow(type = "open",
      length = unit(0.075, "inches"))) +
    ylab("DBH (cm)") +
    xlab("Date")  +
    # facet_wrap(~tree) +
    scale_x_date(
      breaks = seq(as.Date("2014-12-01"), as.Date("2017-01-01"), by="3 months"),
      date_labels = "%Y-%m",  # Format as "YYYY-MM"
      limits = c(min(dbh_imp_df$date), as.Date("2017-04-01"))) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 30, vjust = 0.8))
}

#fit_dbh_sapwood_draws_normal
sap_dbh_points <- function(sapwood_depth_csv, post_slen) {
  d <- read_csv(sapwood_depth_csv)

  calc_quantiles <- function(x, na.rm = TRUE) {
    q <- quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = na.rm)
    return(list(ll = q[1], l = q[2], m = q[3], h = q[4], hh = q[5]))
  }

  pred_sap <- function(x, para) {
    log_mu <- para$alpha + para$beta * log(x)
    exp(log_mu)
  }

  xx <- seq(10, 40, length = 80)

  tmp <- tibble(xx) |>
    mutate(post = list(post_slen)) |>
    mutate(pred = map2(xx, post, pred_sap)) |>
    mutate(dbh = map(pred, calc_quantiles)) |>
    unnest_wider(dbh)

  ggplot(tmp) +
    geom_point(data = d, aes(x = dbh, y = sapwood_depth)) +
    geom_hline(yintercept = 4, lty = 2) +
    geom_hline(yintercept = 6, lty = 2) +
    geom_line(aes(x = xx, y = m)) +
    geom_ribbon(aes(x = xx, ymin = ll, ymax = hh), alpha = 0.4) +
    scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
    xlab("DBH (cm)") +
    ylab("Sapwood length (cm)")  +
  theme_bw()

}

missForest_clean_keep_date <- function(csv, year = 2015, month = 1) {
  d <- read_csv(csv) |>
    janitor::clean_names() |>
    mutate(date = mdy_hm(date)) |>
    mutate(year = year(date)) |>
    mutate(month = month(date)) |>
    filter(month %in% {{month}}) |>
    mutate(day = day(date)) |>
    mutate(yday = yday(date)) |>
    mutate(time = hour(date) * 60 + minute(date)) |>
    mutate(cos_transformed_day = cos((yday - 1) / 365 * 2 * pi)) |>
    mutate(cos_transformed_time = cos((time / 1440) * 2 * pi)) |>
    dplyr::select(year, yday, date, time,
      vpd, par, t01_0_0:t15_0_0) |>
    pivot_longer(c(t01_0_0:t15_0_0), names_to = "id", values_to = "ks") |>
    mutate(tree = str_split_fixed(id, "_", 3)[, 1]) |>
    mutate(dir = str_split_fixed(id, "_", 3)[, 2]) |>
    mutate(dep = str_split_fixed(id, "_", 3)[, 3]) |>
    mutate(dir = case_when(
      dir == "0" ~ "S",
      dir == "1" ~ "E",
      dir == "2" ~ "N",
      dir == "3" ~ "W"
    )) |>
    mutate(dep = case_when(
      dep == "0" ~ 2,
      dep == "1" ~ 4,
      dep == "2" ~ 6,
    )) |>
    mutate(tree = as.factor(tree)) |>
    mutate(dir = as.factor(dir)) |>
    dplyr::select(-id)

   d |>
    filter(year == {{year}}) |>
    as.data.frame()
}

  clean_imputed_df <- function(imputed_df, tree_id, month, day) {
    imputed_df |>
      mutate(date = as.Date(yday - 1, origin = paste0(year, "-01-01"))) |>
      mutate(hour = time %/% 60) |>
      mutate(mins = time %% 60)  |>
      mutate(date_time = as.POSIXct(paste(date, sprintf("%02d:%02d:00", hour, mins)), format="%Y-%m-%d %H:%M:%S")) |>
      filter(yday >= 30 * (month - 1) + day) |>
      filter(yday <  30 * (month - 1) + day + 5) |>
      filter(tree == tree_id) |>
      dplyr::select(date_time, ks, dep, dir) |>
      mutate(model = "Imputed")
  }

  clean_raw_df <- function(rubber_raw_data_csv, year, month, day, tree_id) {
    missForest_clean_keep_date(
      csv = rubber_raw_data_csv,
      year = year,
      month = month) |>
      as_tibble() |>
      filter(yday >= 30 * (month - 1) + day) |>
      filter(yday <  30 * (month - 1) + day + 5) |>
      filter(tree == tree_id) |>
      dplyr::select(date, ks, dep, dir) |>
      mutate(date_time = as.POSIXct(date)) |>
      mutate(model = "Raw")
  }

imp_points <- function(imputed_df_1, rubber_raw_data_csv_1, year_1, month_1, day_1,
                       imputed_df_2, rubber_raw_data_csv_2, year_2, month_2, day_2) {

  # if (depth) tree_id <- "t11" else tree_id <- "t04"
  tree_id <- "t11"

  tmp <- clean_raw_df(rubber_raw_data_csv_1, year = year_1, month = month_1, day = day_1, tree_id = tree_id)
  tmp2 <- clean_imputed_df(imputed_df_1, tree_id, month_1, day_1)
  tmp3 <- clean_raw_df(rubber_raw_data_csv_2, year = year_2, month = month_2, day = day_2, tree_id = tree_id)
  tmp4 <- clean_imputed_df(imputed_df_2, tree_id, month_2, day_2)

  df1 <- bind_rows(tmp, tmp2) |>
    mutate(model = factor(model, levels = c("Raw", "Imputed")))
  df2 <- bind_rows(tmp3, tmp4) |>
    mutate(model = factor(model, levels = c("Raw", "Imputed")))

  dep_fun <- function(data) {
     ggplot(data, aes(x = date_time, y = ks, col = as.factor(dep))) +
      geom_line() +
      geom_point(size = 1) +
      theme_bw() +
      facet_grid(~ model) +
      guides(col = guide_legend(title="Depth (cm)")) +
      theme(legend.position = "bottom") +
      labs(x = "Date", y = "K")
  }

  p <- dep_fun(df1) +
    theme(
      legend.position = "none"
    )
  p2 <- dep_fun(df2)

  p / p2 +
    plot_annotation(tag_levels = "A")

  # if (depth) {
  #   p <- ggplot(tmp3, aes(x = date_time, y = ks, col = as.factor(dep))) +
  #     geom_line() +
  #     geom_point(size = 1) +
  #     theme_bw() +
  #     facet_grid(~ model) +
  #     guides(col = guide_legend(title="Depth (cm)")) +
  #     theme(legend.position = "bottom") +
  #     labs(x = "Date", y = "K")
  # } else {
  #   p <- ggplot(tmp3, aes(x = date_time, y = ks, col = as.factor(dir))) +
  #     geom_line() +
  #     geom_point(size = 1) +
  #     theme_bw() +
  #     facet_grid(~ model) +
  #     guides(col = guide_legend(title="Direction")) +
  #     theme(legend.position = "bottom") +
  #     labs(x = "Date", y = "K")
  # }
  # p
}
