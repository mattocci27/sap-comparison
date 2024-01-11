knn_impute <- function(rubber_raw_data_csv) {
  d <- read_csv(rubber_raw_data_csv) |>
    janitor::clean_names() |>
    mutate(date = mdy_hm(date))
  kNN(d)
}

calc_summary2 <- function(row, s_df) {
  s_df2 <- s_df |>
    mutate(log_fd = row$log_a + row$b * log_k) |>
    # mutate(fd = exp(log_fd)) # scale by sapwood area
    mutate(fd = exp(log_fd) * s) # scale by sapwood area
  return(s_df2)
}

process_full_ab_df <- function(sarea_df, post_dir_dep_mid, dir_dep_imp_df) {
  post_dir_dep_mid_df <- tibble(dir_dep = names(post_dir_dep_mid),
    dir_dep_effects = as.numeric(post_dir_dep_mid))
  dir_dep_imp_df_re <- prepare_dir_dep_imp_df(dir_dep_imp_df, post_dir_dep_mid_df)
  process_full_df(dir_dep_imp_df_re, sarea_df)
}

prepare_summary_stats <- function(full_df_processed, post1, post2) {
  post_ab_df1 <- post1 |>
    summarize(log_a = median(log_a), b = median(b))
  post_ab_df2 <- post2 |>
    summarize(log_a = median(log_a), b = median(b))

  summary_stats1 <- calc_summary2(post_ab_df1, full_df_processed)
  summary_stats2 <- calc_summary2(post_ab_df2, full_df_processed)

  summary_stats1 <- summary_stats1 |>
    group_by(tree, date, time) |>
    summarize(fd = sum(fd, na.rm = TRUE), k = mean(exp(log_k), na.rm = TRUE)) |>
    ungroup()
  summary_stats2 <- summary_stats2 |>
    group_by(tree, date, time) |>
    summarize(fd = sum(fd, na.rm = TRUE), k = mean(exp(log_k), na.rm = TRUE)) |>
    ungroup()


  tmp1 <- summary_stats1 |>
    filter(date == "2015-08-01") #|>
  tmp2 <- summary_stats2 |>
    filter(date == "2015-08-01") #|>
  # log_fd2 <- tmp2 |> pull(log_fd)
  fd2 <- tmp2 |> pull(fd)

  x_max <- max(tmp1$fd, na.rm = T)
  y_max <- max(fd2, na.rm = T)
  max_lab <- max(x_max, y_max)

  list(
    df = tmp1 |> mutate(fd2),
    max_lab = max_lab)
}

ab_example <- function(full_df_processed, post1, post2, summary_stats_list, ab_summarized_df, segments_post) {
  post_ab_df1 <- post1 |>
    summarize(log_a = median(log_a), b = median(b))
  post_ab_df2 <- post2 |>
    summarize(log_a = median(log_a), b = median(b))

  p1 <- ggplot(full_df_processed, aes(x = exp(log_k))) +
    geom_histogram() +
    xlab("K") +
    theme_bw()

  log_k <- seq(0.001, 0.5, length = 100) |> log()
  fd_fun1 <- function(x) exp(post_ab_df1$log_a + post_ab_df1$b * x)
  fd_fun2 <- function(x) exp(post_ab_df2$log_a + post_ab_df2$b * x)
  y1 <- fd_fun1(log_k)
  y2 <- fd_fun2(log_k)

  p2 <- tibble(y1, y2, log_k) |>
    pivot_longer(-log_k) |>
    mutate(pg = ifelse(name == "y1", "0.02", "0.08")) |>
    ggplot(aes(x = exp(log_k), y = value, col = pg)) +
    geom_line() +
    xlab("K") +
    ylab(expression("Sap flux density "(g~m^{-2}~s^{-1}))) +
    labs(color = expression(italic(P[g])~(MPa~m^{-1}))) +
    theme_bw() +
    theme(
      legend.position = "right",
      legend.justification = c(0, 0.5),
      legend.key.width = unit(0.2, "cm"),
      legend.key.height = unit(0.4, "cm")
    )

  max_lab <- summary_stats_list$max_lab
  p3 <- summary_stats_list$df |>
    # mutate(tree_dir = paste0(tree, dir, sep = "_")) |>
    # sample_n(5000) |>
    ggplot(aes(x = fd, y = fd2)) +
    # geom_line(aes( group = tree_dir)) +
    # geom_point(alpha = 0.2) +
    geom_point(alpha = 0.2, aes(color = k)) +
    scale_color_viridis_c(name = "K", option = "C") +
    geom_abline(intercept = 0, slope = 1, lty = 2) +
    ylab(expression(Sap~flux~(g~s^{-1}) * ";"~italic(P[g])==0.08~(MPa ~m^{-1}))) +
    xlab(expression(Sap~flux~(g~s^{-1}) * ";"~italic(P[g])==0.02~(MPa ~m^{-1}))) +
    coord_cartesian(xlim = c(0, max_lab), ylim = c(0, max_lab)) +
    theme_bw() +
    theme(
      legend.position = c(0.85, 0.38),
      legend.key.width = unit(0.2, "cm"),
      legend.key.height = unit(0.4, "cm")
    )

  segments_post2 <- segments_post %>%
    mutate(id = 1:nrow(.)) |>
    mutate(id = as.character(id))
  tmp <- full_join(ab_summarized_df, segments_post2)
  p4 <- ggplot(tmp |> sample_n(1000), aes(x = exp(log_a), y = b, col = tr)) +
    geom_point(alpha = 0.9) +
    scale_color_viridis_c(name = expression("Transpiration (mm"~y^-1*")"), option = "D") +
    scale_x_log10() +
    xlab(expression("Coefficient"~italic(a))) +
    ylab(expression("Coefficient"~italic(b))) +
    theme_bw() +
    theme(
      legend.position = "right",
      legend.key.width = unit(0.2, "cm"),
      legend.key.height = unit(0.4, "cm")
    )

  p1 + p2 + p3 + p4 + plot_layout(ncol = 2, nrow = 2) +
     plot_annotation(tag_levels = "A") &
     theme(
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 7),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      plot.tag = element_text(size = 10)
     )
}
