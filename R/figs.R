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
    legend.title = element_text(size = 9),
    plot.tag = element_text(face = "bold")
  )
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# Color-blind safe palette for main groups
okabe <- c(
  "diffuse" = "#D55E00",   # Vermilion
  "ring"    = "#009E73",   # Bluish Green
  "palm"    = "#56B4E9",   # Sky Blue
  "liana"   = "#CC79A7"    # Reddish Purple
)

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
     levels = c("HB", "HH", "VM", "TG", "AP"))) |>
    filter(calib_type == "Pressure") |>
    ggplot(aes(y = ks, x = as.factor(pressure))) +
    geom_boxplot(width = 0.4) +
    geom_point(alpha = 0.6) +
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

sma_ks <- function(five_spp_csv, ks_trees_csv, log = FALSE)  {
  p1 <- sma_scatter(five_spp_csv, log = log)
  p2 <- ks_box(ks_trees_csv)
  p1 + p2 +
    plot_layout(nrow = 1, width = c(1.8, 1)) +
    plot_annotation(tag_levels = "a") &
    theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12)
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
        "mm", "px"), height = NA, width = NA, dpi = 600, ...) {

  ggsave(
    filename = paste0(filename, ".tiff"),
    plot = plot,
    height = height,
    width = width,
    units = units,
    dpi = dpi,
    compression = "lzw",
    ...
  )

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

  paste0(filename, c(".png", ".pdf", ".tiff"))
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

line_pg_multi <- function(data, xylem_lab, k_range, fit_summary_combined) {

  # data <- tar_read(fd_k_traits_csv)
  xylem_lab2 <- xylem_lab |>
    dplyr::select(sp_short, sp_short_chr, xylem_long_fct)
  d <- read_csv(data) |>
    filter(is.na(removed_k)) |>
    rename(sp_short_chr = sp_short) |>
    left_join(xylem_lab2, by = "sp_short_chr")

  data <- fit_summary_combined |>
    mutate(max_pg = as.character(max_pg))

  # tar_load(k_range)
  # tar_load(xylem_lab)

  xylem_lab <- xylem_lab |>
    mutate(a_chr = str_c("alpha[1,", sp_num,"]")) |>
    mutate(b_chr = str_c("alpha[2,", sp_num,"]"))

  k_range2 <- left_join(k_range, xylem_lab, by = "species")

  log_a <- data |>
    filter(str_detect(variable, "alpha\\[1")) |>
    dplyr::select(variable, log_a = q50, max_pg)

  b <- data |>
    filter(str_detect(variable, "alpha\\[2")) |>
    dplyr::select(variable, b = q50, max_pg)
  sig <- data |>
    filter(variable == "sigma") |>
    dplyr::select(sigma = q50, max_pg)

  log_a2 <- log_a |>
    mutate(b = b$b)

  log_a3 <- full_join(log_a2, sig)

  nd <- left_join(k_range2, log_a3, by = c("a_chr" = "variable", "max_pg" = "max_pg")) |>
    group_by(max_pg, species) |>
    nest() |>
    ungroup()

  # nd$data[[1]]

  pred_data <- nd |>
    mutate(log_xx = map(data, \(x)seq(log(x$k_lwr), log(x$k_upr), length = 80))) |>
    mutate(log_pred = pmap(list(log_xx, data), \(log_xx, data) {data$log_a + data$b * log_xx + data$sigma^2 / 2})) |>
    unnest(cols = c(data, log_xx, log_pred))

  my_cols <- gg_color_hue(4)
  # my_cols <- viridis_pal(option = "E")(4)
  # my_cols <- brewer.pal(n = 4, name = "Paired")

#   my_cols <- c(
#   "diffuse" = "#D55E00",   # Vermilion
#   "ring"    = "#009E73",   # Bluish Green
#   "palm"    = "#56B4E9",   # Sky Blue
#   "liana"   = "#CC79A7"    # Reddish Purple
# )


  ggplot() +
    geom_point(data = d |>
      filter(xylem_type == "DP"), aes(x = k, y = fd), col = okabe[1]) +
    geom_point(data = d |>
      filter(xylem_type == "RP"), aes(x = k, y = fd), col = okabe[2]) +
    geom_point(data = d |>
      filter(xylem_type == "Pa"), aes(x = k, y = fd), col = okabe[3]) +
    geom_point(data = d |>
      filter(xylem_type == "L"), aes(x = k, y = fd), col = okabe[4]) +
    # geom_point(data = d, aes(x = k, y = fd, col = xylem_long_fct)) +
    # geom_line(data = pred_data, aes(x = exp(log_xx), y = exp(log_pred), group = max_pg)) +
    geom_line(data = pred_data, aes(x = exp(log_xx), y = exp(log_pred), col = as.numeric(max_pg), group = max_pg)) +
    facet_wrap(vars(sp_short), ncol = 4, scale = "free") +
    ylab(expression("Sap flux density "(g~m^{-2}~s^{-1}))) +
    # xlab(expression("K "((Delta~T[max]-Delta~T)/Delta~T))) +
    xlab("K") +
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
  # d_0015 <- d |> filter(p_g <= 0.015) |> mutate(max_pg = "0.015")
  d_002 <- d |> filter(p_g <= 0.02) |> mutate(max_pg = "0.02")
  d_0025 <- d |> filter(p_g <= 0.025) |> mutate(max_pg = "0.025")
  d_003 <- d |> filter(p_g <= 0.03) |> mutate(max_pg = "0.03")
  d_0035 <- d |> filter(p_g <= 0.035) |> mutate(max_pg = "0.035")
  d_004 <- d |> filter(p_g <= 0.04) |> mutate(max_pg = "0.04")
  d_005 <- d |> filter(p_g <= 0.05) |> mutate(max_pg = "0.05")
  d_006 <- d |> filter(p_g <= 0.06) |> mutate(max_pg = "0.06")
  d_007 <- d |> filter(p_g <= 0.07) |> mutate(max_pg = "0.07")
  d_008 <- d |> filter(p_g <= 0.08) |> mutate(max_pg = "0.08")

  k_data <- bind_rows(d_002, d_0025, d_003, d_0035, d_004, d_005, d_006, d_007, d_008)
  k_data2 <- k_data |>
    group_by(species, max_pg) |>
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
  #   filter(max_pg == "0.02")
  # k_range |>
  #   filter(max_pg == "0.025")

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
    dplyr::select(sp_short, sp_short_chr, xylem_long_fct)

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
    xlab("K") +
    scale_color_manual(values = unname(okabe)) +
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

# draws <- tar_read(fit_draws_segments_xylem_0.08)
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

  y_lab_raw <- data_a$sp_short2 |> levels()
  y_lab <- y_lab_raw
  names(y_lab) <- y_lab_raw

  species_list <- y_lab[-1:-4]

  replace_with_italic <- function(label) {
    for (species in species_list) {
      if (label == species) {
        return(substitute(italic(SPECIES), list(SPECIES = as.name(species))))
      }
    }
    return(label) # If not replaced, return original label
  }

  y_lab <- sapply(y_lab, replace_with_italic)

  p1 <- ggplot(data_a, aes(x = exp(value), y = sp_short2, fill = xylem))  +
    facet_grid(group ~ ., scales = "free", space = "free") +
    scale_x_log10(
        breaks = c(10^2, 10^3, 10^4),
        labels = trans_format("log10", math_format(10^.x))) +
    scale_fill_manual(values = unname(okabe)) +
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
    scale_fill_manual(values = unname(okabe)) +
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
    plot_annotation(tag_levels = "a")

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
ab_pg_ribbon <- function(xylem_lab, k_range, fit_summary_combined, coef_a = TRUE) {

  data <- fit_summary_combined |>
    mutate(max_pg = as.character(max_pg))

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

  fig_data <- right_join(fig_data, k_range, by = c("species", "max_pg")) |>
    mutate(max_pg = as.numeric(max_pg))

  p <- ggplot(fig_data, aes(x = max_pg, y = q50, group = sp_short,
    fill = xylem_long_fct)) +
      facet_wrap(~sp_short, ncol = 4, scale = "free_y") +
      geom_ribbon(aes(ymin = q2_5, ymax = q97_5), alpha = 0.4) +
      geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.6) +
      geom_point() +
      geom_line() +
      scale_fill_manual(name = "", values = unname(okabe)) +
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
    geom_errorbar(aes(ymin = q2.5_seg, ymax = q97.5_seg), alpha = 0.6) +
    geom_errorbar(aes(xmin = q2.5_pool, xmax = q97.5_pool), alpha = 0.6) +
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
    geom_errorbar(aes(ymin = q2.5_seg, ymax = q97.5_seg), alpha = 0.6) +
    geom_errorbar(aes(xmin = q2.5_pool, xmax = q97.5_pool), alpha = 0.6) +
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
      plot_annotation(tag_levels = "a")
}

# Extract percentiles for model 3 and 4
extract_alpha_percentiles <- function(data, index) {
  data %>%
    janitor::clean_names() %>%
    filter(str_detect(variable, paste0("alpha\\[", index))) %>%
    dplyr::select(q2_5, q50, q97_5)
}
extract_A_percentiles <- function(data, index) {
  data %>%
    janitor::clean_names() %>%
    filter(str_detect(variable, paste0("A\\[", index))) %>%
    dplyr::select(q2_5, q50, q97_5)
}


# Helper function to prepare data frames
ab_points_model4_prepare_df <- function(data, a_values, b_values) {
  bind_cols(
    a_values |> rename(a_q2_5 = q2_5, a_q50 = q50, a_q97_5 = q97_5),
    b_values |> rename(b_q2_5 = q2_5, b_q50 = q50, b_q97_5 = q97_5)
  ) |>
  left_join(data, by = "sample_id")
}

# Helper function to create plots
ab_points_model4_create_plot <- function(df, pub_df, title, with_pub = FALSE) {
  base_plot <- ggplot() +
    scale_x_log10() +
    my_theme() +  # Adjusted theme for simplicity
    labs(col = "",
      x = expression(paste("Co-", italic(a))),
      y = expression(paste("Co-", italic(b))))  +
    theme(
      legend.background = element_rect(fill = "transparent", color = NA)
    )

  h7 <- scales::hue_pal()(7)
  # h4 <- scales::hue_pal()(4)
  h4 <- brewer.pal(n = 4, name = "Paired")

  h4 <- c(
  "diffuse" = "#D55E00",   # Vermilion
  "ring"    = "#009E73",   # Bluish Green
  "palm"    = "#56B4E9",   # Sky Blue
  "liana"   = "#CC79A7"    # Reddish Purple
  )

  my_col <- h7
  my_col <- rep("black", 7)
  my_col[1] <- h4[1]
  my_col[3] <- h4[2]
  my_col[4] <- h4[3]
  my_col[6] <- h4[4]

   # Define shapes
  shapes <- c(
    DP = 1,
    RP = 1,
    Pa = 1,
    Li = 1,
    Ba = 21,
    He = 22,
    NP = 23
  )

  df <- df |>
    mutate(xylem_type = ifelse(xylem_type == "L", "Li", xylem_type)) |>
    mutate(xylem_fct = factor(xylem_type, levels = c("DP", "RP", "Pa", "Li")))

  if (with_pub) {

    # Levels: Ba DP He Li NP Pa RP
    # Levels: DP L Pa RP
    base_plot +
      #geom_point(data = pub_df, aes(x = a, y = b, color = type), shape = 1) +
      geom_point(data = pub_df, aes(x = a, y = b, color = type)) +
      scale_shape_manual(values = shapes) +
      geom_point(data = df, aes(x = a_q50, y = b_q50, color = xylem_fct)) +
      geom_sma(data = pub_df, aes(x = a, y = b), method = "sma", se = TRUE, col = "grey40") +
      geom_sma(data = df, aes(x = a_q50, y = b_q50), method = "sma", se = FALSE) +
      scale_colour_manual(
        values = c(
          DP = my_col[1],
          RP = my_col[3],
          Pa = my_col[4],
          Li = my_col[6],
          # L = my_col[6],
          Ba = my_col[2],
          NP = my_col[5],
          He = my_col[7]
        )
      ) +
      # scale_shape_manual(values = shapes) +
      theme(legend.position = c(0.1, 0.7))
  } else {
    base_plot +
      geom_point(data = df, aes(x = a_q50, y = b_q50, col = xylem_fct)) +
      geom_errorbar(data = df, aes(x = a_q50, ymin = b_q2_5, ymax = b_q97_5, col = xylem_fct), alpha = 0.6, show.legend = FALSE) +
      geom_errorbar(data = df, aes(xmin = a_q2_5, xmax = a_q97_5, y = b_q50, col = xylem_fct), alpha = 0.6, show.legend = FALSE) +
      geom_sma(data = df, aes(x = a_q50, y = b_q50), method = "sma", se = TRUE) +
      scale_colour_manual(
        values = c(
          DP = my_col[1],
          RP = my_col[3],
          Pa = my_col[4],
          Li = my_col[6]
        )
      ) +
      theme(legend.position = c(0.15, 0.75))
  }
}

ab_points_model4_create_plot <- function(df, pub_df, title, with_pub = FALSE) {

  pub_df <- pub_df |> mutate(source = "pub") |>
    filter(references != "this study") |>
    mutate(
      type = factor(type, levels = c("DP", "RP", "Pa", "Li", "Ba", "He", "NP"))
    )

  df <- df |> mutate(source = "ours")

  base_plot <- ggplot() +
    scale_x_log10() +
    my_theme() +
    labs(
      col = "",
      x = expression(paste("Co-", italic(a))),
      y = expression(paste("Co-", italic(b)))
    ) +
    theme(
      legend.background = element_rect(fill = "transparent", color = NA)
    )

  # Color-blind safe palette for main groups
  my_col <- c(
    "diffuse" = "#D55E00",   # Vermilion
    "ring"    = "#009E73",   # Bluish Green
    "palm"    = "#56B4E9",   # Sky Blue
    "liana"   = "#CC79A7"    # Reddish Purple
  )

  df <- df |>
    mutate(xylem_type = ifelse(xylem_type == "L", "Li", xylem_type)) |>
    mutate(xylem_fct = factor(xylem_type, levels = c("DP", "RP", "Pa", "Li")))

  dp_df <- df |> filter(xylem_fct == "DP")
  rp_df <- df |> filter(xylem_fct == "RP")
  pa_df <- df |> filter(xylem_fct == "Pa")
  li_df <- df |> filter(xylem_fct == "Li")

  if (with_pub) {
    base_plot +
      geom_point(data = pub_df, aes(x = a, y = b, col = type, shape = type),
        size = case_when(
          pub_df$type == "NP" ~ 2,
          pub_df$type == "Ba" ~ 2,
          pub_df$type == "He" ~ 2,
          TRUE ~ 2.5
          )) +
      scale_shape_manual(values = c(21, 22, 23, 24, 3, 4, 8)) +
      scale_color_manual(values = c(unname(my_col), "black", "black", "black")) +
      scale_x_log10() +
      geom_point(data = dp_df, aes(x = a_q50, y = b_q50), size = 2.5, fill = my_col[1], shape = 21, alpha = 0.6) +
      geom_point(data = rp_df, aes(x = a_q50, y = b_q50), size = 2.5, fill = my_col[2], shape = 22, alpha = 0.6) +
      geom_point(data = pa_df, aes(x = a_q50, y = b_q50), size = 2.5, fill = my_col[3], shape = 23, alpha = 0.6) +
      geom_point(data = li_df, aes(x = a_q50, y = b_q50), size = 2.5, fill = my_col[4], shape = 24, alpha = 0.6) +
      geom_sma(data = pub_df, aes(x = a, y = b), method = "sma", se = TRUE, col = "grey40") +
      geom_sma(data = df, aes(x = a_q50, y = b_q50), method = "sma", se = FALSE) +
      guides(
        color = guide_legend(
          override.aes = list(
            shape = c(21, 22, 23, 24, 3, 4, 8),
            fill = c(my_col[1], my_col[2], my_col[3], my_col[4], "black", "black", "black"),
            color = rep("black", 7),
            size = c(rep(2.5, 4), rep(2, 3)),
            alpha = c(rep(.5, 4), rep(1, 3))
            )),
         shape = "none"
      ) +
      theme(legend.position = c(0.1, 0.7))
  } else {
    base_plot +
      geom_errorbar(data = df, aes(x = a_q50, ymin = b_q2_5, ymax = b_q97_5, col = xylem_fct), alpha = 0.6, show.legend = FALSE) +
      geom_errorbar(data = df, aes(xmin = a_q2_5, xmax = a_q97_5, y = b_q50, col = xylem_fct), alpha = 0.6, show.legend = FALSE) +
      geom_point(data = df, aes(x = a_q50, y = b_q50, fill = xylem_fct, shape = xylem_fct), color = "black", size = 2.5, alpha = 0.6) +
      geom_sma(data = df, aes(x = a_q50, y = b_q50), method = "sma", se = TRUE) +
      scale_fill_manual(values = unname(my_col)) +
      scale_colour_manual(values = unname(my_col)) +
      scale_shape_manual(values = c(21, 22, 23, 24)) +
      guides(
        fill = guide_legend(
          override.aes = list(
            shape = c(21, 22, 23, 24),
            fill = c(my_col[1], my_col[2], my_col[3], my_col[4]),
            color = rep("black", 4)
            )),
        shape = "none",
        color = "none"
      ) +
      theme(legend.position = c(0.15, 0.75), legend.title = element_blank())
  }
}

# post hoc
ab_points_model4_sma <- function(summary, fd_k_traits_csv, xylem_lab, pub_ab_path, rm_dip = TRUE) {
  a_sp <- exp(extract_alpha_percentiles(summary, 1))
  b_sp <- extract_alpha_percentiles(summary, 2)
  a_seg <- exp(extract_A_percentiles(summary, 1))
  b_seg <- extract_A_percentiles(summary, 2)

  sample_id <- read_csv(fd_k_traits_csv) |>
    filter(is.na(removed_k)) |>
    select(species, sample_id) |>
    distinct()

  sp_df <- xylem_lab |>
    arrange(species) |>
    bind_cols(
      a_sp |> rename(a_q2_5 = q2_5, a_q50 = q50, a_q97_5 = q97_5),
      b_sp |> rename(b_q2_5 = q2_5, b_q50 = q50, b_q97_5 = q97_5)
      )
  seg_df <- sample_id |>
    bind_cols(
      a_seg |> rename(a_q2_5 = q2_5, a_q50 = q50, a_q97_5 = q97_5),
      b_seg |> rename(b_q2_5 = q2_5, b_q50 = q50, b_q97_5 = q97_5)
      ) |>
    left_join(xylem_lab)

  pub_df <- read_csv(pub_ab_path) |>
    janitor::clean_names()

  # pub_df2 <- pub_df |>
  #   filter(b < 2.5) |>
  #   filter(b < 2 | !str_detect(references, "Dix"))

  # Linear models and equations
  fit_sp <- smatr::sma(b_q50 ~ log(a_q50), data = sp_df)
  fit_seg <- smatr::sma(b_q50 ~ log(a_q50), data = seg_df)
  fit_pub <- smatr::sma(b ~ log(a), data = pub_df)

  list(
    fit_sp,
    fit_seg,
    fit_pub
  )
}
#
ab_points_model4 <- function(summary, fd_k_traits_csv, xylem_lab, pub_ab_path, rm_dip = TRUE) {

  # summary <- tar_read(fit_summary_segments_xylem_0.08)
  # tar_load(fd_k_traits_csv)
  # tar_load(xylem_lab)
  # Prepare data
  a_sp <- exp(extract_alpha_percentiles(summary, 1))
  b_sp <- extract_alpha_percentiles(summary, 2)
  a_seg <- exp(extract_A_percentiles(summary, 1))
  b_seg <- extract_A_percentiles(summary, 2)

  sample_id <- read_csv(fd_k_traits_csv) |>
    filter(is.na(removed_k)) |>
    dplyr::select(species, sample_id) |>
    distinct()

  # Create data frames
  # sp_df <- ab_points_model4_prepare_df(xylem_lab, a_sp, b_sp)
  # seg_df <- ab_points_model4_prepare_df(sample_id, a_seg, b_seg)
  sp_df <- xylem_lab |>
    arrange(species) |>
    bind_cols(
      a_sp |> rename(a_q2_5 = q2_5, a_q50 = q50, a_q97_5 = q97_5),
      b_sp |> rename(b_q2_5 = q2_5, b_q50 = q50, b_q97_5 = q97_5)
      )
  seg_df <- sample_id |>
    bind_cols(
      a_seg |> rename(a_q2_5 = q2_5, a_q50 = q50, a_q97_5 = q97_5),
      b_seg |> rename(b_q2_5 = q2_5, b_q50 = q50, b_q97_5 = q97_5)
      ) |>
    left_join(xylem_lab)

  # Linear models and equations
  fit_sp <- lm(b_q50 ~ log(a_q50), data = sp_df)
  fit_seg <- lm(b_q50 ~ log(a_q50), data = seg_df)
  eq_sp <- sprintf("y=%.3f+%.3flnx", coef(fit_sp)[1], coef(fit_sp)[2])
  eq_seg <- sprintf("y=%.3f+%.3flnx", coef(fit_seg)[1], coef(fit_seg)[2])

  pub_df <- read_csv(pub_ab_path) |>
  # pub_df <- read_csv("data-raw/pub_ab.csv") |>
    janitor::clean_names()

  # Create plots
  p1 <- ab_points_model4_create_plot(sp_df, pub_df, with_pub = FALSE)  #+
      # annotate("point", x = 119, y = 1.231, colour = "red", size = 3)

  p2 <- ab_points_model4_create_plot(sp_df, pub_df, with_pub = TRUE)

  # Combine plots
  p1 + p2 + plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold"))
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
  d <- summary12
  log_a <- extract_percentiles(d$fit_species, "log_a") |> exp()
  b <- extract_percentiles(d$fit_species, "b")
  log_a2 <- extract_percentiles(d$fit_segments, "log_a") |> exp()
  b2 <- extract_percentiles(d$fit_segments, "b")

# Model 3
  d3 <- summary3
  log_a3 <- extract_alpha_percentiles(d3, 1) |> exp()
  b3 <- extract_alpha_percentiles(d3, 2)

# Model 4
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
      ggplot(aes(!!x_q50, !!y_q50, fill = xylem_long_fct, color = xylem_long_fct, shape = xylem_long_fct)) +
      geom_abline(slope = 1, intercept = 0, lty = 2) +
      geom_errorbar(aes(ymin = !!y_q2_5, ymax = !!y_q97_5), alpha = 0.6) +
      geom_errorbar(aes(xmin = !!x_q2_5, xmax = !!x_q97_5), alpha = 0.6) +
      geom_point(alpha = 0.6, color = "black", stroke = 0.25, size = 2.5) +
      scale_fill_manual(values = unname(okabe)) +
      scale_color_manual(values = unname(okabe)) +
      scale_shape_manual(values = c(21, 22, 23, 24)) +
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
    plot_annotation(tag_levels = "a") &
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
      limits = c(min(dbh_imp_df$date), as.Date("2017-05-01"))) +
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
    pivot_longer(c(t01_0_0:t15_0_0), names_to = "id", values_to = "k") |>
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

clean_imputed_df_tmp <- function(imputed_df, tree_id, month, day) {
  imputed_df |>
    mutate(date = as.Date(yday - 1, origin = paste0(year, "-01-01"))) |>
    mutate(hour = time %/% 60) |>
    mutate(mins = time %% 60)  |>
    mutate(date_time = as.POSIXct(paste(date, sprintf("%02d:%02d:00", hour, mins)), format="%Y-%m-%d %H:%M:%S")) |>
    filter(yday >= 30 * (month - 1) + day) |>
    filter(yday <  30 * (month - 1) + day + 5) |>
    filter(tree == tree_id) |>
    dplyr::select(date_time, k, dep, dir) |>
    mutate(model = "Imputed")
}

clean_raw_df_tmp <- function(rubber_raw_data_csv, year, month, day, tree_id) {
  missForest_clean_keep_date(
    csv = rubber_raw_data_csv,
    year = year,
    month = month) |>
    as_tibble() |>
    filter(yday >= 30 * (month - 1) + day) |>
    filter(yday <  30 * (month - 1) + day + 5) |>
    filter(tree == tree_id) |>
    dplyr::select(date, k, dep, dir) |>
    mutate(date_time = as.POSIXct(date)) |>
    mutate(model = "Raw")
}

prepare_imp2_df <- function(rubber_raw_data_csv, combined_imputed_k_mapped) {
  d <- read_csv(rubber_raw_data_csv) |>
    janitor::clean_names() |>
    mutate(date = mdy_hm(date)) |>
    mutate(year = year(date)) |>
    mutate(month = month(date)) |>
    filter(!(year == 2016 & (month == 12 | month == 10 | month == 9))) |>
    mutate(day = day(date)) |>
    mutate(yday = yday(date)) |>
    mutate(time = hour(date) * 60 + minute(date)) |>
    dplyr::select(year, month, day, time,
      vpd, par, t01_0_0:t16_0_0) |>
    pivot_longer(c(t01_0_0:t16_0_0), names_to = "id", values_to = "k") |>
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
    rename(k_ori = k) |>
    mutate(k_1st_imputed = combined_imputed_k_mapped$k_new_without_na) |>
    mutate(k_1st_imputed_with_na = combined_imputed_k_mapped$k_new_with_na) |>
    mutate(k_2nd_imputed = combined_imputed_k_mapped$k_imp)
}


imp_points2 <- function(imp2_df, year = 2015, month = 2, day = 1, days = 10, dep = 2, dir = "S", tree = "t11") {
  imp2_df_re <- imp2_df |>
    filter(year == {{year}}) |>
    filter(month == {{month}}) |>
    filter(tree == {{tree}}) |>
    filter(dir == {{dir}}) |>
    filter(dep == {{dep}}) |>
    filter(is.na(k_1st_imputed_with_na)) |>
    mutate(date = mdy_hm(paste(month, day, year, time %/% 60, time %% 60))) |>
    mutate(date_time = as.POSIXct(date)) |>
    pivot_longer(c(k_ori, k_2nd_imputed), names_to = "model", values_to = "k") |>
    mutate(model = case_when(
      model == "k_ori" ~ "Observed",
      model == "k_2nd_imputed" ~ "Re-imputed"
    )) |>
    filter(day <= {{days}} - {{day}} + 1)

  ggplot(imp2_df_re, aes(x = date_time, y = k, col = as.factor(model))) +
    geom_line() +
    geom_point(size = 1) +
    theme_bw() +
    guides(col = guide_legend(title=NULL)) +
    theme(legend.position = "bottom") +
    labs(x = "Feburary 2015", y = expression(italic(K))) # +
    # scale_x_datetime(date_labels = "%b %d\n%Y", date_breaks = "1 days")  # Use scale_x_datetime for POSIXct
}

imp_points <- function(imputed_df_1, rubber_raw_data_csv_1, year_1, month_1, day_1,
                       imputed_df_2, rubber_raw_data_csv_2, year_2, month_2, day_2) {

  # if (depth) tree_id <- "t11" else tree_id <- "t04"
  tree_id <- "t11"

  tmp <- clean_raw_df_tmp(rubber_raw_data_csv_1, year = year_1, month = month_1, day = day_1, tree_id = tree_id)
  tmp2 <- clean_imputed_df_tmp(imputed_df_1, tree_id, month_1, day_1)
  tmp3 <- clean_raw_df_tmp(rubber_raw_data_csv_2, year = year_2, month = month_2, day = day_2, tree_id = tree_id)
  tmp4 <- clean_imputed_df_tmp(imputed_df_2, tree_id, month_2, day_2)

  df1 <- bind_rows(tmp, tmp2) |>
    mutate(model = factor(model, levels = c("Raw", "Imputed")))
  df2 <- bind_rows(tmp3, tmp4) |>
    mutate(model = factor(model, levels = c("Raw", "Imputed")))

  dep_fun <- function(data) {
     ggplot(data, aes(x = date_time, y = k, col = as.factor(dep))) +
      geom_line() +
      geom_point(size = 1) +
      theme_bw() +
      facet_grid(~ model) +
      guides(col = guide_legend(title="Depth (cm)")) +
      theme(legend.position = "bottom") +
      labs(x = "Date", y = expression(italic(K)))
  }

  p <- dep_fun(df1) +
    theme(
      legend.position = "none"
    )
  p2 <- dep_fun(df2)

  p / p2 +
    plot_annotation(tag_levels = "a")
}


# Define a function for the repeated ggplot elements
tr_bar <- function(data, granier_df, x, y) {

  new_labels <- c(
    "dep_only"      = "Dep",
    # "dir_dep"  = expression(Dep %*% Dir),      # × rendered via matrix product
    "dir_dep"  = "Dep×Dir",      # × rendered via matrix product
    "dir_only"      = "Dir",
    "sapwood_aera"    = expression(A[s]),
    "total"    = "Total"
  )
  # col_pal <- viridis::viridis(n = 7)
  # colors <- c("dir_only" = col_pal[2], "dir_dep" = col_pal[3], "dep_only" = col_pal[5], "sarea" = col_pal[6], "total" = col_pal[7])
  granier_val <- granier_df  |>
    pull(tr_m)

  ggplot(data, aes(x = {{ x }}, y = {{ y }})) +
  # ggplot(data, aes(x = {{ x }}, y = {{ y }}, fill = {{fill}}, group = {{ group }})) +
    geom_bar(stat = "identity", position = "dodge", fill = "lightblue") +
    geom_errorbar(aes(ymin = tr_l, ymax = tr_h), position = position_dodge(.9), width = 0.2) +
    scale_y_continuous(breaks = c(0, 250, 500, 750, 1000)) +
    scale_x_discrete(labels = new_labels) +
    labs(fill = "B) Source") +
    geom_hline(yintercept = granier_val, lty = 2) +
    ylab(expression(paste("Transpiration (mm y"^-1~")"))) +
    xlab(expression(Applied~italic(P[g])~(MPa~m^{-1}))) +
    theme_bw() +
    theme(
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.4, "cm"),
      plot.tag = element_text(face = "bold")
    )
}

# Define a function for the repeated ggplot elements
tr_bar_ab <- function(data, granier_df, x, y, fill, group, model4 = TRUE) {
  col_pal <- viridis::mako(n = 5)
  granier_val <- granier_df |>
    pull(tr_m)

  if (model4) {
    ggplot(data, aes(x = {{ x }}, y = {{ y }})) +
      geom_bar(stat = "identity", position = "dodge", fill = col_pal[2]) +
      geom_errorbar(aes(ymin = tr_l, ymax = tr_h), position = position_dodge(.9), width = 0.2) +
      scale_y_continuous(breaks = c(0, 250, 500, 750, 1000)) +
      geom_hline(yintercept = granier_val, lty = 2) +
      ylab(expression(paste("Transpiration (mm y"^-1~")"))) +
      xlab(expression(Applied~italic(P[g])~(MPa~m^{-1}))) +
      theme_bw() +
      theme(
        legend.title = element_blank()
      )
  } else {
    ggplot(data, aes(x = {{ x }}, y = {{ y }}, fill = {{ fill }}, group = {{ group }})) +
      geom_bar(stat = "identity", position = "dodge") +
      geom_errorbar(aes(ymin = tr_l, ymax = tr_h), position = position_dodge(.9), width = 0.2) +
      scale_y_continuous(breaks = c(0, 250, 500, 750, 1000)) +
      scale_fill_manual(values = rev(col_pal[-1])) +
      geom_hline(yintercept = granier_val, lty = 2) +
      ylab(expression(paste("Transpiration (mm y"^-1~")"))) +
      xlab(expression(Applied~italic(P[g])~(MPa~m^{-1}))) +
      # xlab(expression(
      #    atop("Statistical models", italic(P[g]) == 0.08 ~ "(" * "MPa m"^{-1} * ")")
      #    )) +
      theme_bw() +
      theme(
        # legend.title = element_blank()
      )

  }

}

generate_tr_bar_ab_df <- function(ab_uncertainty_combined_df, ab_granier_uncertainty_combined_df) {
# Read and prepare granier_df
    granier_df <- ab_granier_uncertainty_combined_df |>
      mutate(
        id = "granier",
        across(c(tr_l, tr_h), ~ NA_real_) # if you have multiple such columns, use across
      )

# Read ab_uncertainty_combined_df and bind with granier_df
    tmp <- ab_uncertainty_combined_df |>
      bind_rows(granier_df) |>
      mutate(
        pg = as.factor(str_extract(id, "\\d+\\.\\d+")),
        model = str_extract(id, "species_xylem|segments_xylem|segments_inclusive|species_only|granier"),
        model = case_when(
          model == "species_only" ~ "Model 1",
          model == "segments_inclusive" ~ "Model 2",
          model == "species_xylem" ~ "Model 3",
          model == "segments_xylem" ~ "Model 4",
          model == "granier" ~ "Granier",
        ),
        model = as.factor(model)
      )
}

generate_tr_bar_all_df <- function(total_uncertainty_combined_df, sarea_uncertainty_combined_df, dir_dep_uncertainty_combined_df) {

  total_df <- total_uncertainty_combined_df |>
  mutate(id = "total")

  sarea_df <- sarea_uncertainty_combined_df |>
    mutate(id = "sapwood_aera")

  dir_dep_df <- dir_dep_uncertainty_combined_df |>
    mutate(id = str_extract(id, "dir_only|dep_only|dir_dep$"))

  bind_rows(total_df, sarea_df, dir_dep_df)
}


generate_rel_cont_df <- function(total_uncertainty_combined_df, ab_uncertainty_df, sarea_uncertainty_combined_df, dir_dep_uncertainty_combined_df) {
  total_df <- total_uncertainty_combined_df |>
    mutate(id = "total") |>
    mutate(var = tr_sd^2) |>
    dplyr::select(id, everything())
  sarea_df <- sarea_uncertainty_combined_df |>
    mutate(id = "sarea") |>
    mutate(var = tr_sd^2) |>
    dplyr::select(id, everything())
  ab_df <- ab_uncertainty_df |>
    mutate(id = "ab") |>
    mutate(var = tr_sd^2) |>
    dplyr::select(id, everything())

# Read in your data frame
  dir_dep_df <- dir_dep_uncertainty_combined_df |>
    mutate(
      id = str_extract(id, "dir_only|dep_only|dir_dep$"),
      var = tr_sd^2
    )

# Compute the covariance directly
  cov_dir_dep <- (dir_dep_df$var[dir_dep_df$id == "dir_dep"] -
                  (dir_dep_df$var[dir_dep_df$id == "dir_only"] +
                   dir_dep_df$var[dir_dep_df$id == "dep_only"]))

# Add the covariance as a new row
  dir_dep_df <- dir_dep_df %>%
    add_row(id = "dir_dep_cov", tr_m = NA, tr_l = NA, tr_h = NA,
            tr_mean = NA, tr_sd = NA, var = cov_dir_dep)

  total_var <- total_df |>
    pull(var)

  bind_rows(sarea_df, ab_df, dir_dep_df) |>
    filter(id != "dir_dep") |>
    mutate(rel = var / total_var * 100) |>
    mutate(rel_adj = (rel / sum(rel)) * 100) |>
    mutate(id = factor(id, levels = c("sarea", "ab", "dir_only",  "dir_dep_cov", "dep_only")))
}


rel_bar <- function(rel_bar_df) {
  new_labels <- c(
    "sarea" = "Sapwood area",
    "ab" = "Coefficient a-b",
    "dir_only" = "Direction",
    "dep_only" = "Depth",
    "dir_dep_cov" = "Cov(Direction, Depth)")

  col_pal <- viridis::viridis(n = 7)
  colors <- c("dir_only" = col_pal[2], "dir_dep_cov" = col_pal[4], "dep_only" = col_pal[5], "sarea" = col_pal[6], "ab" = col_pal[1])

  ggplot(rel_bar_df, aes(x = factor(1), y = rel_adj, fill = id)) +
    geom_bar(stat = "identity") +
    # scale_fill_viridis_d(labels = new_labels) +
    scale_fill_manual(labels = new_labels, values = colors) +
    # scale_fill_viridis_d() +
    ylab("Relative contribution (%)") +
    xlab("") +
    labs(fill = "C) Source") +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
}

rel_bar_wide <- function(rel_bar_df) {
  new_labels <- c(
    "sarea" = "Sapwood area",
    "ab" = "Coefficient a-b",
    "dir_only" = "Direction",
    "dep_only" = "Depth",
    "dir_dep_cov" = "Direction and depth")

  rel_bar_df <- rel_bar_df |>
    mutate(id =  reorder(id, -rel_adj))

  ggplot(rel_bar_df, aes(x = id, y = rel_adj, fill = id)) +
    geom_bar(stat = "identity") +
    scale_fill_viridis_d(labels = new_labels) +
    ylab("Relative contribution (%)") +
    xlab("") +
    labs(fill = "Source") +
    theme_bw() +
    theme(
      # axis.ticks.x = element_blank(),
      # legend.title = element_blank(),
      axis.text.x = element_blank(),
      legend.position = c(0.7, 0.7),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.4, "cm")
    )
}


generate_ab_pg_data <- function(s, d, xylem_lab, a = TRUE) {
  if (a) {
    alpha_vars <- str_c("alpha_1_", 1:31)
    alpha_s_vars <- str_c("alpha[1,", 1:31, "]")
  } else {
    alpha_vars <- str_c("alpha_2_", 1:31)
    alpha_s_vars <- str_c("alpha[2,", 1:31, "]")
  }

# Initialize a list to store results
  results <- tibble()

# Loop over each alpha variable
  for (i in seq_along(alpha_vars)) {
    # Selecting relevant columns from d for the current alpha
    d2 <- d |>
      dplyr::select(matches(paste0(alpha_vars[i], "$")), max_pg)

    # Filtering s for the current alpha summary
    s2 <- s |>
      filter(variable == alpha_s_vars[i])

    # Find max_pg that maximizes and minimizes q50
    max_pg_maximize_q50 <- s2 |>
      arrange(desc(q50)) |>
      slice(1) |>
      pull(max_pg)

    max_pg_minimize_q50 <- s2 |>
      arrange(q50) |>
      slice(1) |>
      pull(max_pg)

    # Process the data
    d3 <- d2 |>
      filter(max_pg %in% c(max_pg_maximize_q50, max_pg_minimize_q50)) |>
      mutate(max_pg = ifelse(max_pg == max_pg_maximize_q50, "max", "min"))

    d4 <- d3 |>
      pivot_wider(names_from = max_pg, values_from = all_of(alpha_vars[i])) |>
      unnest(cols = c(min, max)) |>
      mutate(diff = max - min)

    summary_stats <- d4 |>
      summarize(m = median(diff), ll = quantile(diff, 0.025), hh = quantile(diff, 0.975))

    # Store results
    results <- bind_rows(results, summary_stats)
    # results[[alpha_vars[i]]] <- summary_stats
  }

  results <- results |>
    mutate(sp = 1:31) |>
    mutate(sp = str_c("sp", sp))

  xylem_lab2 <-  xylem_lab |>
    mutate(sp = str_c("sp", sp_num))

  full_join(results, xylem_lab2)
}

ab_pg_summary_bars <- function(s, d, xylem_lab) {
  a_df <- generate_ab_pg_data(s, d, xylem_lab, a = TRUE) |>
    mutate(across(where(is.numeric), exp))
  b_df <- generate_ab_pg_data(s, d, xylem_lab, a = FALSE)

  p1 <- a_df |>
    mutate(sp_short = factor(sp_short, levels = sp_short[order(m)])) |>
    ggplot(aes(y = sp_short, x = m, fill = xylem_long_fct, color = xylem_long_fct, shape = xylem_long_fct)) +
    geom_errorbar(aes(xmin = ll, xmax = hh), width = 0.2) +
    geom_point() +
    scale_x_continuous(breaks = c(0, 1, 5, 10)) +
    xlab(expression(italic(a) ~ " range" ~ "(=" ~ italic(a[max]) / italic(a[min]) ~ ")")) +
    geom_vline(xintercept = 1, lty = 2, col = "grey30") +
    my_theme() +
    theme(
       legend.position = "none")

  p2 <- b_df |>
    mutate(sp_short = factor(sp_short, levels = sp_short[order(m)])) |>
    ggplot(aes(y = sp_short, x = m, fill = xylem_long_fct, color = xylem_long_fct, shape = xylem_long_fct)) +
    geom_errorbar(aes(xmin = ll, xmax = hh), width = 0.2) +
    geom_point() +
    xlab(expression(italic(b) ~ " range" ~ " = " * "(" * italic(b[max]) - italic(b[min]) * ")")) +
    geom_vline(xintercept = 0, lty = 2, col = "grey30") +
    my_theme() +
    theme(
       legend.position = "none")

  p3 <- a_df |>
    mutate(sp_short = factor(sp_short, levels = sp_short[order(m)])) |>
    ggplot(aes(y = sp_short, x = m, fill = xylem_long_fct, color = xylem_long_fct, shape = xylem_long_fct)) +
    geom_point() +
    geom_errorbar(aes(xmin = ll, xmax = hh), width = 0.2) +
    guides(color = guide_legend(ncol = 4, title = "")) #+
    my_theme() +
    theme(
      plot.background = element_blank(),
      legend.background = element_blank(),
      plot.margin = margin(0, 0, 0, 0),
      legend.margin = margin(0, 0, 0, 0),
      legend.box = element_blank()
      )
  # extracted_legend <- g_legend(p3)

  extracted_legend <- cowplot::get_legend(
    p1 +
    scale_color_manual(values = unname(okabe)) +
    guides(
        fill = guide_legend(
          ncol = 4, title = "",
          override.aes = list(
            shape = c(21, 22, 23, 24),
            fill = unname(okabe),
            # color = rep("black", 4)
            color = unname(okabe)
            )),
        shape = "none",
        color = "none"
      ) +
    #  guides(color = guide_legend(
    #   ncol = 4, title = "")) +
      theme(legend.text = element_text(size = 9),
        legend.position = "bottom"))

  p <- p1 + p2 +
    plot_annotation(tag_levels = "a") &
    scale_fill_manual(values = unname(okabe)) &
    scale_shape_manual(values = c(21, 22, 23, 24)) &
    theme(
      plot.margin = margin(2, 2, 0, 2),
      axis.text.y = element_text(face = "italic", size = 8),
      axis.title.y = element_blank(),
      strip.background = element_blank(),
      strip.text.y = element_blank()
    )
  cowplot::plot_grid(p, extracted_legend, ncol = 1, rel_heights = c(1, .1))
}

generate_tr_example_list <- function(post_dir_dep_mid, dir_dep_imp_df, sarea_df, segments_xylem_post_ab_fit_draws_segments_xylem_0.02, segments_xylem_post_ab_fit_draws_segments_xylem_0.08) {
  post_dir_dep_mid_df <- tibble(dir_dep = names(post_dir_dep_mid),
    dir_dep_effects = as.numeric(post_dir_dep_mid))
  dir_dep_imp_df_re <- prepare_dir_dep_imp_df(dir_dep_imp_df, post_dir_dep_mid_df)
  full_df_processed <- process_full_df(dir_dep_imp_df_re, sarea_df)

  # fd vs k
  tmp1 <- segments_xylem_post_ab_fit_draws_segments_xylem_0.02
  tmp2 <- segments_xylem_post_ab_fit_draws_segments_xylem_0.08
  post_ab_df1 <- tmp1 |>
    summarize(log_a = median(log_a), b = median(b))
  post_ab_df2 <- tmp2 |>
    summarize(log_a = median(log_a), b = median(b))

  log_k <- seq(0.001, 0.5, length = 100)  |> log()
  fd_fun1 <- function(x) exp(post_ab_df1$log_a + post_ab_df1$b * x)
  fd_fun2 <- function(x) exp(post_ab_df2$log_a + post_ab_df2$b * x)
  y1 <- fd_fun1(log_k)
  y2 <- fd_fun2(log_k)

  # fd
  calc_summary2 <- function(row, s_df) {
  s_df2 <- s_df |>
    mutate(log_fd = row$log_a + row$b * log_k) |>
    # mutate(fd = exp(log_fd)) # scale by sapwood area
    mutate(fd = exp(log_fd) * s) # scale by sapwood area
    return(s_df2)
  }

  summary_stats1 <- calc_summary2(post_ab_df1, full_df_processed)
  summary_stats2 <- calc_summary2(post_ab_df2, full_df_processed)

  tmp1 <- summary_stats1 |>
    # filter(date == "2015-01-01") |>
    filter(tree == "t01")
  tmp2 <- summary_stats2 |>
    # filter(date == "2015-01-01") |>
    filter(tree == "t01")
  log_fd2 <- tmp2 |> pull(log_fd)
  fd2 <- tmp2 |> pull(fd)

  list(
    df1 = full_df_processed |> dplyr::select(log_k),
    df2 = tibble(y1, y2, log_k),
    df3 = tmp1 |> mutate(fd2)
   )
}

tr_example_panel <- function(tr_example_list) {
  p1 <- ggplot(tr_example_list$df1, aes(x = exp(log_k))) +
    geom_histogram() +
    xlab(expression(italic(K))) +
    theme_bw()

  p2 <- tr_example_list$df2 |>
    pivot_longer(-log_k) |>
    mutate(pg = ifelse(name == "y1", "0.02", "0.08")) |>
    ggplot(aes(x = exp(log_k), y = value, col = pg)) +
    geom_line() +
    xlab(expression(italic(K))) +
    ylab(expression(italic(F[d])~(g~m^{-2}~s^{-1}))) +
    guides(col = guide_legend(title = expression(italic(P[g])~(MPa~m^{-1})))) +
    theme_bw() +
    theme(
      legend.position = c(0.3, 0.72),
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 7),
      legend.margin = margin(1, 1, 1, 1)
    )

  df3 <- tr_example_list$df3 #|>
    # filter(str_detect(date, "2015-07|2015-08"))

  p3 <- ggplot(df3, aes(x = fd * 600 * 1e-7, y = fd2 * 600 * 1e-7)) +
    geom_point() +
    # geom_pointdensity() +
    geom_abline(intercept = 0, slope = 1, lty = 2) +
    scale_color_viridis() +
    xlab("Sap flux; Pg == 0.02") +
    ylab("Sap flux; Pg = 0.08") +
    # coord_cartesian(xlim = c(0, 5000), ylim = c(0, 5000)) +
    theme_bw()

  p1 + p2 + p3
}

imp_r2_scatter2 <- function(data = combined_imputed_k_mapped, draws, granier = FALSE) {
  # library(targets)
  # library(tidyverse)
  # data <- tar_read(combined_imputed_k_mapped)
  # draws <- tar_read(segments_xylem_post_ab_fit_draws_segments_xylem_0.08)

  post_ab <- draws |> apply(2, median)

  if (granier) {
    post_ab[1] <- log(119)
    post_ab[2] <- 1.231
  }

  df <- data |>
    filter(is.na(k_new_with_na)) |>
    filter(!is.na(k_ori)) |>
    mutate(log_fd_ori = post_ab["log_a"] + post_ab["b"] * log(k_ori)) |>
    mutate(log_fd_reimp = post_ab["log_a"] + post_ab["b"] * log(k_imp)) |>
    mutate(flux_obs = exp(log_fd_ori)) |>
    mutate(flux_reimp = exp(log_fd_reimp))

  dens <- with(df, kde2d(flux_obs, flux_reimp, n = 500))
  ix <- findInterval(df$flux_obs, dens$x)
  iy <- findInterval(df$flux_reimp, dens$y)
  df$density <- dens$z[cbind(ix, iy)]

  p <- ggplot(df, aes(x = flux_obs, y = flux_reimp, col = density)) +
    geom_point(alpha = 0.01) +
    scale_color_viridis_c(option = "D") +
    geom_abline(intercept = 0, slope = 1, linetype = 2) +
    theme_bw()
  p
}

imp_r2_scatter <- function(data = combined_imputed_k_mapped) {
  df <- data |>
    filter(is.na(k_new_with_na)) |>
    filter(!is.na(k_ori))
    # filter(k_imp < 0.1 & k_new_without_na < 0.01)
  # fit <- lm(k_imp ~ k_new_without_na, data = df)
  r2 <- round(cor(df$k_imp, df$k_new_without_na)^2, 2)
  n <- nrow(df)

  p1 <- ggplot(df, aes(x = k_new_without_na, y = k_imp)) +
    # geom_point(alpha = 0.05) +
    geom_bin2d(bins = 100) +
    # scale_fill_viridis(option = "D", direction = 1) +
    scale_fill_viridis_c(
      option = "D",
      direction = 1,
      breaks = c(1, 100000),  # Only two breaks, similar to the map legend
      labels = c("1", "100,000"),
      name = "Number of data",  # Change the label
      guide = guide_colorbar(
        barwidth = 8,
        barheight = 0.5,
        title.position = "top",
        label.position = "top",
        title.hjust = 0.5,
        label.vjust = 2,  # Adjust vertical position to avoid overlap
        label.theme = element_text(size = 8, angle = 0)  # Fine-tune text size and angle
     )
    ) +
   annotate(
      "text",
      label = paste("italic(R)^2 == ", r2, "*','~italic(N) == ", n),
      x = 0, y = 1.2,
      size = 3,
      hjust = 0,
      vjust = 0,
      parse = TRUE
    ) +
    geom_abline(intercept = 0, slope = 1, linetype = 2) +
    # geom_smooth(method = "lm", se = FALSE) +
    # stat_cor(
    #   aes(label = paste(..rr.label.., ..n.label.., sep = "~`,`~")),
    #   show.legend = FALSE
    # ) +
    labs(
      x = expression(Observed~italic(K)),
      y = expression(Re*"-"*imputed~italic(K))) +
    my_theme() +
    theme(
      legend.position = "top",
      legend.margin = margin(t = -10, unit = "pt"),  # Reduce space between legend and plot
      plot.margin = margin(t = 0, r = 5, b = 5, l = 5, unit = "pt")  # Adjust plot margins as needed
    )

  p2 <- ggplot(df, aes(x = k_new_without_na, y = k_imp)) +
    geom_bin2d(bins = 100) +
    scale_fill_viridis_c(
      option = "D",
      direction = 1,
      name = "Number of data"  # Change the label
    ) +
    geom_abline(intercept = 0, slope = 1, linetype = 2) +
    # geom_smooth(method = "lm", se = FALSE) +
    labs(
      x = expression(Observed~italic(K)),
      y = expression(Re*"-"*imputed~italic(K))) +
    my_theme() +
    coord_cartesian(xlim = c(0, 0.1), ylim = c(0, 0.1)) +
    theme(
      legend.position = "none",
      axis.text.y = element_text(size = 6),
      axis.text.x = element_text(size = 6, angle = 60, margin = margin(t = 5)),
      axis.title = element_blank())

  p1 +
    annotation_custom(ggplotGrob(p2),
                    xmin = 0.8, xmax = 1.4,
                    ymin = -0.1, ymax = 0.5)
}

# p2 <- p1 +
#   annotation_custom(ggplotGrob(p_slope),
#                     xmin = 0.3, xmax = 1.05,
#                     ymin = 0.4, ymax = 1.05)

generate_reimp_bin_df <- function(combined_imputed_k_mapped, draws, granier = FALSE) {
  df <- combined_imputed_k_mapped |>
    filter(is.na(k_new_with_na)) |>
    filter(!is.na(k_ori))

  post_ab <- draws |> apply(2, median)

  if (granier) {
    post_ab[1] <- log(119)
    post_ab[2] <- 1.231
  }

  df <- df  |>
    mutate(log_fd_ori = post_ab["log_a"] + post_ab["b"] * log(k_ori)) |>
    mutate(log_fd_reimp = post_ab["log_a"] + post_ab["b"] * log(k_imp)) |>
    mutate(flux_obs = exp(log_fd_ori)) |>
    mutate(flux_reimp = exp(log_fd_reimp))

  # Define intervals (bins)
  bins <- seq(min(df$flux_obs), max(df$flux_obs), length.out = 21)

  # Calculate total flux in each interval
  df |>
    mutate(interval = cut(flux_obs, breaks = bins, include.lowest = TRUE)) |>
    group_by(interval) |>
    summarise(
      total_flux_obs = sum(flux_obs),
      total_flux_reimp = sum(flux_reimp),
      avg_flux_obs = mean(flux_obs),
      n = n()
    ) |>
    # mutate(difference = total_flux_obs - total_flux_reimp)
    mutate(relative_difference = total_flux_reimp / total_flux_obs)
}

generate_reimp_bin_ci_list <- function(reimp_bin_df, combined_imputed_k_mapped, draws) {
  df <- combined_imputed_k_mapped |>
    filter(is.na(k_new_with_na)) |>
    filter(!is.na(k_ori))
  post_ab <- as.data.frame(draws)

  calculate_relative_difference <- function(i) {
    df <- df |>
      mutate(
        log_fd_ori = post_ab[i, "log_a"] + post_ab[i, "b"] * log(k_ori),
        log_fd_reimp = post_ab[i, "log_a"] + post_ab[i, "b"] * log(k_imp),
        flux_obs = exp(log_fd_ori),
        flux_reimp = exp(log_fd_reimp)
      )

    bins <- seq(min(df$flux_obs, na.rm = TRUE), max(df$flux_obs, na.rm = TRUE), length.out = 21)

    # Calculate total flux in each interval
    df |>
      mutate(interval = cut(flux_obs, breaks = bins, include.lowest = TRUE)) |>
      group_by(interval) |>
      summarise(
        total_flux_obs = sum(flux_obs, na.rm = TRUE),
        total_flux_reimp = sum(flux_reimp, na.rm = TRUE),
        avg_flux_obs = mean(flux_obs, na.rm = TRUE),
        n = n(),
        .groups = 'drop'
      ) |>
      mutate(relative_difference = total_flux_reimp / total_flux_obs) #|>
      # pull(relative_difference)
  }

  # Use map to iterate over the posterior samples and calculate relative differences

  nd <- tibble(data = map(1:nrow(post_ab), calculate_relative_difference)) |>
    mutate(rel_diff = map(data, pull, relative_difference)) |>
    mutate(total_changes = map_dbl(data, \(x) {
      total_flux_reimp <- x |> pull(total_flux_reimp) |> sum()
      total_flux_obs <- x |> pull(total_flux_obs) |> sum()
      tmp <- 1 - (total_flux_reimp / total_flux_obs)
      tmp * 100
    }))


  rel_diff_mat <- nd$rel_diff %>%
    map(~ as.numeric(.)) %>%  # Ensure each element is numeric
    reduce(cbind)

  total_changes <- nd |>
    pull(total_changes) |>
    quantile(c(0.025, 0.5, 0.975))

  # Calculate quantiles for each interval
  tmp <- apply(rel_diff_mat, 1, quantile, c(0.025, 0.5, 0.975)) |> t() |> as_tibble() |> janitor::clean_names()
  list(reimp_bin_ci_df = bind_cols(reimp_bin_df, tmp), total_changes = total_changes)

}


reimp_bin_bar <- function(reimp_bin_df, error = TRUE) {
  # Plot the differences
  p <- ggplot(reimp_bin_df, aes(x = avg_flux_obs, y = relative_difference)) +
    geom_bar(stat = "identity") +
    # geom_errorbar(aes(ymin = x2_5_percent, ymax = x97_5_percent), linewidth = 0.25, width = 0.2) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    coord_cartesian(ylim = c(0, 1.35)) +
    labs(
      # x = "Midpoint of average observed flux bins",
      x = expression("Sap flux density bins "(g~m^{-2}~s^{-1})),
      y = "Relative difference\n(re-imputed / observed)") +
    theme_bw()

  if (error) {
    p + geom_errorbar(aes(ymin = x2_5_percent, ymax = x97_5_percent), linewidth = 0.25, width = 0.2)
  } else {
    p
  }
}

# fig. 6

ab_violin <- function(pub_ab_csv) {

  d <- read_csv(pub_ab_csv) |>
    janitor::clean_names() |>
    filter(!is.na(a) & !is.na(b)) |>
    mutate(type = factor(type, levels = c("DP", "RP", "NP", "Pa","Li", "Ba", "He")))

  nn <- d |>
    group_by(type) |>
    summarise(
      n = n()
    )

  # p_df <- d |>
  #   group_by(type) |>
  #   nest() |>
  #   filter(type != "He") |>
  #   mutate(p_a = map_dbl(data, \(x) t.test(log(x$a), mu = log(119))$p.val)) |>
  #   mutate(p_b = map_dbl(data, \(x) t.test(x$b, mu = 1.231, alternaive = "less")$p.val))

  # p_df |>
  #   dplyr::select(type, p_a, p_b) |>
  #   mutate(p_a = round(p_a, 3)) |>
  #   mutate(p_b = round(p_b, 3)) |>
  #   arrange(type) |>
  #   as.data.frame()
  #     type   p_a   p_b
  # 1   DP 0.000 0.753
  # 2   RP 0.001 0.238
  # 3   NP 0.141 0.432
  # 4   Pa 0.797 0.026
  # 5   Li 0.000 0.001
  # 6   Ba 0.105 0.197


  p_a <- tibble(
    type = c("DP", "RP", "NP", "Pa", "Li", "Ba", "He"),
    p_label = c("italic(P)<0.001", "italic(P)<0.001", "italic(P)==0.141",
                "italic(P)==0.797", "italic(P)<0.001", "italic(P)==0.105", "italic(Na)"),
    y_pos = c(1e6, 1e6, 1e6, 1e6, 1e6, 1e6, 1e6)  # adjust based on your plot scale
  )

  p_b <- tibble(
    type = c("DP", "RP", "NP", "Pa", "Li", "Ba", "He"),
    p_label = c("italic(P)==0.753", "italic(P)==0.238", "italic(P)==0.432",
                "italic(P)==0.026", "italic(P)<0.001", "italic(P)==0.197", "italic(Na)"),
    y_pos = rep(3.2, 7)  # adjust based on your plot scale
  )

  n_b <- tibble(
    type = c("DP", "RP", "NP", "Pa", "Li", "Ba", "He"),
    n_label = c("italic(N)==61", "italic(N)==18", "italic(N)==18",
                "italic(N)==9", "italic(N)==8", "italic(N)==4", "italic(N)==1"),
    y_pos = rep(0, 7)  # adjust based on your plot scale
  )

  p1 <- ggplot(d, aes(x = type, y = a)) +
    geom_violin(aes(fill = type, col = type), trim = FALSE, alpha = 0.6) +
    geom_jitter(aes(fill = type), width = 0.1, size = 1, alpha = 0.6, shape = 21) +
    geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.8) + # nolint
    geom_text(data = p_a, aes(x = type, y = y_pos, label = p_label),
      parse = TRUE, size = 3.5, vjust = 0) +
    scale_fill_manual(values = c(unname(okabe)[1:2], "grey70", unname(okabe)[3:4], "grey70", "grey70")) +
    scale_color_manual(values = c(unname(okabe)[1:2], "grey70", unname(okabe)[3:4], "grey70", "grey70")) +
    geom_hline(yintercept = 119, lty = 2, col = "grey40") +
    scale_y_log10(
      breaks = c(10, 1e2, 1e3, 1e4, 1e5, 1e6),
      labels = c(expression(10^1), expression(10^2),
        expression(10^3), expression(10^4),
        expression(10^5), expression(10^6))
    ) +
    ylab(expression(Coefficient~italic(a))) +
    my_theme() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      legend.position = "none"
      )

  p2 <- ggplot(d, aes(x = type, y = b)) +
    geom_violin(aes(fill = type, col = type), trim = FALSE, alpha = 0.6) +
    geom_jitter(aes(fill = type), width = 0.1, size = 1, alpha = 0.6, shape = 21) +
    geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.8) + # nolint
    geom_text(data = p_b, aes(x = type, y = y_pos, label = p_label),
      parse = TRUE, size = 3.5, vjust = 0) +
    geom_text(data = n_b, aes(x = type, y = y_pos, label = n_label),
      parse = TRUE, size = 3.5, vjust = 0) +
    scale_fill_manual(values = c(unname(okabe)[1:2], "grey70", unname(okabe)[3:4], "grey70", "grey70")) +
    scale_color_manual(values = c(unname(okabe)[1:2], "grey70", unname(okabe)[3:4], "grey70", "grey70")) +
    geom_hline(yintercept = 1.231, lty = 2, col = "grey40") +
    ylab(expression(Coefficient~italic(b))) +
    my_theme() +
    theme(
      axis.title.x = element_blank(),
      legend.position = "none"
      )

  p1 / p2 +
    plot_annotation(tag_levels = "a") &
    theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12)
    )

}



