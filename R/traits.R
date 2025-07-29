# Create a function to prepare the data
prepare_plot_data <- function(data, value_prefix) {
  data |>
    dplyr::select(-contains(setdiff(c("a_", "b_"), value_prefix))) |>
    rename_with(~ sub(value_prefix, "", .x)) |>
    mutate(ab = substr(value_prefix, 1, 1))  |>
    mutate(
      mid = ifelse(ab == "a", exp(mid), mid),
      lwr = ifelse(ab == "a", exp(lwr), lwr),
      upr = ifelse(ab == "a", exp(upr), upr)
    )
}

# Create a function to prepare the data
prepare_line_data <- function(data, ab = "a") {
  if (ab == "a") {
    data <- data |>
      mutate(
        ymin = exp(pred_a_ll),
        ymax = exp(pred_a_hh),
        y = exp(pred_a_m)
      )
  } else {
    data <- data |>
      mutate(
        ymin = pred_b_ll,
        ymax = pred_b_hh,
        y = pred_b_m
      )
  }
}

# Create a function for the repeated plotting code
plot_data <- function(data, r2_list, ab_value, plot_title = NULL, inner_tag = "A", sp = FALSE) {

  data <- data |>
    filter(ab == ab_value)

  if (!is.null(r2_list)) {
    r2 <- lapply(r2_list, "[[", ifelse(ab_value == "a", 1 , 2)) |>
      sapply(function(x) format(round(x[[2]], 2), nsmall = 2))
    r2_fmt <- paste0("italic(R^2) == \"", r2, "\"")
  } else {
    r2_fmt <- rep(NA, 4)
  }

  r2_data <- data |>
    group_by(trait) |>
    mutate(r2_x = max(exp(val)), tag_x = min(exp(val)), r2_y = min(mid)) |>
    ungroup() |>
    dplyr::select(trait, tag_x, r2_x, r2_y) |>
    distinct() |>
    mutate(r2 = r2_fmt) |>
    mutate(r2_y = ifelse(ab_value == "a", 10, 0)) |>
    mutate(tag_y = ifelse(ab_value == "a", 1e+4, 4)) |>
    arrange(trait) |>
    mutate(inner_tag = inner_tag)

  if (sp) {
    r2_data <- r2_data |>
      mutate(tag_y = ifelse(ab_value == "a", 8000, 2))
  }

 p <- data |>
    ggplot() +
      geom_point(aes(y = mid, x = exp(val), color = xylem_long_fct), alpha = 0.5) +
      geom_errorbar(aes(x = exp(val), ymin = lwr, ymax = upr, color = xylem_long_fct), linewidth = 0.25) +
      geom_text(data = r2_data, aes(x = r2_x, y = r2_y, label = r2),
        parse = TRUE, size = 3, hjust = 1, vjust = 0) +
      geom_text(data = r2_data, aes(x = tag_x, y = tag_y, label = inner_tag),
        parse = TRUE, size = 3, hjust = 0, vjust = -0.25) +
      facet_grid(. ~ trait, scales = "free") +
      scale_x_log10() +
      scale_color_manual(values = unname(okabe)) +
      my_theme() +
      theme(
        axis.title.x = element_blank(),
        strip.background.x = element_blank(),
        strip.text.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 12),
        plot.margin = margin(t = 0, r = 2, b = 2, l = 0)
      )

  if (!is.null(plot_title)) {
    p <- p + ggtitle(plot_title)
  }

  if (ab_value == "a") {
    p <- p +
      ylab(expression(Coefficient~italic(a))) +
      scale_y_log10(
        breaks = c(10, 100, 1000, 10000),
        labels = c(expression(10^1), expression(10^2), expression(10^3), expression(10^4))
      ) +
      coord_cartesian(ylim = c(10, 20000)) +
      theme(
        strip.text.x = element_blank(),
        axis.ticks.x = element_blank()
        )
  } else {
    p <- p +
      coord_cartesian(ylim = c(0, 4.3)) +
      ylab(expression(Coefficient~italic(b))) +
      theme(strip.text.x = element_blank())
  }

  p
}

add_lines <- function(p, data) {
  p +
    geom_ribbon(data = data, aes(x = x, ymin = ymin, ymax = ymax), alpha = 0.6, fill = "grey") +
    geom_line(data = data, aes(x = x, y = y), linewidth = 0.4, col = "grey40")
}

create_line <- function(pred_data, ymin_trans, ymax_trans, y_trans) {
  bind_rows(
    pred_data[[1]]$pred_line |> mutate(trait = "log_vaf"),
    pred_data[[2]]$pred_line |> mutate(trait = "log_ks")
  ) |>
    mutate(
      ymin = {{ymin_trans}},
      ymax = {{ymax_trans}},
      y = {{y_trans}}
    ) |>
    mutate(trait = fct_relevel(trait, "log_vaf", "log_ks"))
}

create_line2 <- function(pred_data, ymin_trans, ymax_trans, y_trans) {
  bind_rows(
    pred_data[[5]]$pred_line |> mutate(trait = "log_dh")
    # pred_data[[6]]$pred_line |> mutate(trait = "log_vf")
  ) |>
    mutate(
      ymin = {{ymin_trans}},
      ymax = {{ymax_trans}},
      y = {{y_trans}}
    ) |>
    mutate(trait = fct_relevel(trait, "log_dh", "log_vf"))
}

# Function to extract legend
g_legend <- function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Add the extracted legend to the bottom of the combined plot
prepare_x_fake_lab <- function(fig_data) {
  tmp <- fig_data |>
    dplyr::select(mid, xylem_long_fct, trait, val)

  tmp2 <- tmp |>
    group_by(trait) |>
    summarize(text_x = median(val), mid = median(mid)) |>
    ungroup() |>
    mutate(mid = mean(mid)) |>
    mutate(trait_lab = case_when(
      trait == "log_vaf" ~ "VAF~(`%`)",
      trait == "log_ks" ~ "K[S]~(kg~m^{-1}~s^{-1}~MPa^{-1})",
      trait == "log_vf" ~ "VF~(no.~mm^{-2})",
      trait == "log_dh" ~ "D[h]~(µm)",
      trait == "log_swc" ~ "SWC~(`%`)",
      trait == "wood_density" ~ "rho~(g~cm^{-3})"))

  ggplot(tmp2, aes(x = val, y = mid)) +
    geom_text(aes(x = text_x, y = mid, label = trait_lab), parse = TRUE, size = 3, hjust = 0.5, vjust = 0.4) +
    facet_grid(. ~ trait, scales = "free") +
    theme_void() +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      strip.text.x = element_blank(),
      panel.background = element_blank(),
      plot.margin = margin(t = 0.5, r = 0, b = 0, l = 0)
    )
}

traits_sp_points_main <- function(pred_data_seg, pred_data_sp, vaf_r2, ks_r2, vaf_sp_r2, ks_sp_r2) {

  r2_list <- list(vaf_r2, ks_r2)
  r2_sp_list <- list(vaf_sp_r2, ks_sp_r2)

#Apply the data preparation function to segment data
  tmp_seg_vaf <- prepare_plot_data(pred_data_seg[[1]]$pred_points, "a_") |> mutate(trait = "log_vaf")
  tmp_seg_vaf_b <- prepare_plot_data(pred_data_seg[[1]]$pred_points, "b_") |> mutate(trait = "log_vaf")
  tmp_seg_ks <- prepare_plot_data(pred_data_seg[[2]]$pred_points, "a_") |> mutate(trait = "log_ks")
  tmp_seg_ks_b <- prepare_plot_data(pred_data_seg[[2]]$pred_points, "b_") |> mutate(trait = "log_ks")

  fig_data_seg <- bind_rows(
    tmp_seg_vaf,
    tmp_seg_vaf_b,
    tmp_seg_ks,
    tmp_seg_ks_b) |>
    mutate(trait = fct_relevel(trait, "log_vaf", "log_ks")) |>
    mutate(val = ifelse(trait == "log_vaf", log_vaf, log_ks))

  tmp_sp_vaf <- prepare_plot_data(pred_data_sp[[1]]$pred_points, "a_") |> mutate(trait = "log_vaf")
  tmp_sp_vaf_b <- prepare_plot_data(pred_data_sp[[1]]$pred_points, "b_") |> mutate(trait = "log_vaf")
  tmp_sp_ks <- prepare_plot_data(pred_data_sp[[2]]$pred_points, "a_") |> mutate(trait = "log_ks")
  tmp_sp_ks_b <- prepare_plot_data(pred_data_sp[[2]]$pred_points, "b_") |> mutate(trait = "log_ks")

  fig_data_sp <- bind_rows(
    tmp_sp_vaf,
    tmp_sp_vaf_b,
    tmp_sp_ks,
    tmp_sp_ks_b) |>
    mutate(trait = fct_relevel(trait, "log_vaf", "log_ks")) |>
    mutate(val = ifelse(trait == "log_vaf", log_vaf, log_ks))

  line_a_sp <- create_line(pred_data_sp, exp(pred_a_ll), exp(pred_a_hh), exp(pred_a_m))
  line_b_sp <- create_line(pred_data_sp, pred_b_ll, pred_b_hh, pred_b_m)
  line_a_seg <- create_line(pred_data_seg, exp(pred_a_ll), exp(pred_a_hh), exp(pred_a_m))
  line_b_seg <- create_line(pred_data_seg, pred_b_ll, pred_b_hh, pred_b_m)

  p1 <- plot_data(fig_data_seg, r2_list, "a", "Segments", inner_tag = c("A", "B"))
  p1 <- add_lines(p1, data = line_a_seg) +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 7)
    )

  p2 <- plot_data(fig_data_seg, r2_list, "b", inner_tag = c("E", "F"))
  p2 <- add_lines(p2, data = line_b_seg) +
    theme(
      # strip.text.x = element_blank()
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 7)
    )

  p3 <- plot_data(fig_data_sp, r2_sp_list, "a", "Species", inner_tag = c("C", "D"))
  p3 <- add_lines(p3, data = line_a_sp) +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.text = element_text(size = 7),
    )

  p4 <- plot_data(fig_data_sp, r2_sp_list, "b", inner_tag = c("G", "H"))
  p4 <- add_lines(p4, data = line_b_sp) +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 7),
      axis.ticks.y = element_blank(),
      legend.position = "none",
      legend.title = element_blank(),
      legend.box.margin = margin(t = 0, b = 0, unit = "pt"),  # Adjust top and bottom margin of the legend box
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0)  # Adjust the space around the individual legend items
    )

# Create a dummy plot which will be used only to extract the legend
  p5 <- ggplot(fig_data_sp, aes(x = val, y = mid, color = xylem_long_fct)) +
    geom_point() +
    scale_color_manual(values = unname(okabe)) +
    theme_void() +
    theme(
      legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.spacing.x = unit(-0.2, "lines"),
      legend.title = element_blank())

# Extract the legend from p5
  legend <- g_legend(p5)

  x_fake_lab_sp <- prepare_x_fake_lab(fig_data_sp)
  x_fake_lab_seg <- prepare_x_fake_lab(fig_data_seg)

# Combine the plots
  combined_plot <- p1 + p3 + p2 + p4 + plot_layout(nrow = 2, ncol = 2)
  combined_plot / (x_fake_lab_seg + x_fake_lab_sp) /  legend +
    plot_layout(heights = c(1, 0.05, 0.05))

}


traits_seg_points_si <- function(pred_data_seg, pred_data_sp) {
#Apply the data preparation function to segment data
  tmp_seg_wd <- prepare_plot_data(pred_data_seg[[3]]$pred_points, "a_") |> mutate(trait = "wood_density")
  tmp_seg_wd_b <- prepare_plot_data(pred_data_seg[[3]]$pred_points, "b_") |> mutate(trait = "wood_density")
  tmp_seg_swc <- prepare_plot_data(pred_data_seg[[4]]$pred_points, "a_") |> mutate(trait = "log_swc")
  tmp_seg_swc_b <- prepare_plot_data(pred_data_seg[[4]]$pred_points, "b_") |> mutate(trait = "log_swc")
  tmp_seg_dh <- prepare_plot_data(pred_data_seg[[5]]$pred_points, "a_") |> mutate(trait = "log_dh")
  tmp_seg_dh_b <- prepare_plot_data(pred_data_seg[[5]]$pred_points, "b_") |> mutate(trait = "log_dh")
  tmp_seg_vf <- prepare_plot_data(pred_data_seg[[6]]$pred_points, "a_") |> mutate(trait = "log_vf")
  tmp_seg_vf_b <- prepare_plot_data(pred_data_seg[[6]]$pred_points, "b_") |> mutate(trait = "log_vf")

  fig_data_seg <- bind_rows(
    tmp_seg_wd,
    tmp_seg_wd_b,
    tmp_seg_swc,
    tmp_seg_swc_b,
    tmp_seg_vf,
    tmp_seg_vf_b,
    tmp_seg_dh,
    tmp_seg_dh_b) |>
    mutate(trait = fct_relevel(trait, "wood_density", "log_swc", "log_dh", "log_vf")) |>
    mutate(val = case_when(
      trait == "wood_density" ~ log(wood_density),
      trait == "log_swc" ~ log_swc,
      trait == "log_dh" ~ log_dh,
      trait == "log_vf" ~ log_vf,
    ))


  p1 <- plot_data(fig_data_seg, r2_list = NULL, "a", "Segments", inner_tag = LETTERS[1:4]) +
    theme(
      axis.text.x = element_blank(),
      strip.text.x = element_blank()
    )
  p2 <- plot_data(fig_data_seg, r2_list = NULL, "b", inner_tag = LETTERS[5:8]) +
    theme(
      strip.text.x = element_blank()
    )

# Create a dummy plot which will be used only to extract the legend
  p3 <- ggplot(fig_data_seg, aes(x = val, y = mid, color = xylem_long_fct)) +
    geom_point() +
    theme_void() +
    theme(
      legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.spacing.x = unit(0, "lines"),
      legend.title = element_blank())

# Extract the legend from p5
  legend <- g_legend(p3)
  x_fake_lab_seg <- prepare_x_fake_lab(fig_data_seg)

  p1 + p2 + x_fake_lab_seg + legend + plot_layout(nrow = 4, ncol = 1) +
    plot_layout(heights = c(0.5, 0.5, 0.06, 0.05))
}

traits_points_si <- function(pred_data, r2_list = NULL, title = "Segments", sp = FALSE) {
# r2_sp_list <- list(vaf_sp_r2, ks_sp_r2)
#Apply the data preparation function to spment data

  # pred_data <- tar_read(trait_pred_data_noxylem_sp_combined)

  tmp_wd <- prepare_plot_data(pred_data[[3]]$pred_points, "a_") |> mutate(trait = "wood_density")
  tmp_wd_b <- prepare_plot_data(pred_data[[3]]$pred_points, "b_") |> mutate(trait = "wood_density")
  tmp_swc <- prepare_plot_data(pred_data[[4]]$pred_points, "a_") |> mutate(trait = "log_swc")
  tmp_swc_b <- prepare_plot_data(pred_data[[4]]$pred_points, "b_") |> mutate(trait = "log_swc")
  tmp_dh <- prepare_plot_data(pred_data[[5]]$pred_points, "a_") |> mutate(trait = "log_dh")
  tmp_dh_b <- prepare_plot_data(pred_data[[5]]$pred_points, "b_") |> mutate(trait = "log_dh")
  tmp_vf <- prepare_plot_data(pred_data[[6]]$pred_points, "a_") |> mutate(trait = "log_vf")
  tmp_vf_b <- prepare_plot_data(pred_data[[6]]$pred_points, "b_") |> mutate(trait = "log_vf")

  fig_data <- bind_rows(
    tmp_wd,
    tmp_wd_b,
    tmp_swc,
    tmp_swc_b,
    tmp_dh,
    tmp_dh_b,
    tmp_vf,
    tmp_vf_b) |>
    mutate(trait = fct_relevel(trait, "wood_density", "log_swc", "log_dh", "log_vf")) |>
    mutate(val = case_when(
      trait == "wood_density" ~ log(wood_density),
      trait == "log_swc" ~ log_swc,
      trait == "log_dh" ~ log_dh,
      trait == "log_vf" ~ log_vf,
    ))

  line_a_sp <- create_line2(pred_data, exp(pred_a_ll), exp(pred_a_hh), exp(pred_a_m))
  line_b_sp <- create_line2(pred_data, pred_b_ll, pred_b_hh, pred_b_m) |>
    filter(trait == "log_vf")

  p1 <- plot_data(fig_data, r2_list = r2_list, "a", title, inner_tag = LETTERS[1:4], sp = sp) +
    theme(
      axis.text.x = element_blank(),
      strip.text.x = element_blank()
    )
  p2 <- plot_data(fig_data, r2_list = r2_list, "b", inner_tag = LETTERS[5:8], sp = sp) +
    theme(
      strip.text.x = element_blank()
    )

  if (sp) {
    p1 <- add_lines(p1, data = line_a_sp) +
      # coord_cartesian(ylim = c(10, 15667))
      coord_cartesian(ylim = c(10, 20000))
    p2 <- add_lines(p2, data = line_b_sp) +
      coord_cartesian(ylim = c(0, 2.15))
      # coord_cartesian(ylim = c(0, 4.3))
  }

# Create a dummy plot which will be used only to extract the legend
  p3 <- ggplot(fig_data, aes(x = val, y = mid, color = xylem_long_fct)) +
    geom_point() +
    theme_void() +
    theme(
      legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.spacing.x = unit(0, "lines"),
      legend.title = element_blank())

# Extract the legend from p5
  legend <- g_legend(p3)
  x_fake_lab <- prepare_x_fake_lab(fig_data)

  p1 + p2 + x_fake_lab + legend + plot_layout(nrow = 4, ncol = 1) +
    plot_layout(heights = c(0.5, 0.5, 0.06, 0.05))
}
