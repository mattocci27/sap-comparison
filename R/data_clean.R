#' @title write_csv for targets
#' @inheritParams readr::write_csv
my_write_csv <- function(x, path, append = FALSE, col_names = !append) {
    write_csv(x, path, append = FALSE, col_names = !append)
    paste(path)
}

clean_ks_trees <- function(data) {
  #d2 <- read_csv("data-raw/ks_pres_tens_trees_raw.csv") |>
  d2 <- read_csv(data) |>
    pivot_longer(p_2:p_8, names_to = "pressure", values_to = "ks") |>
    mutate(pressure = case_when(
      pressure == "p_2" ~ 0.02,
      pressure == "p_5"  ~ 0.05,
      pressure == "p_8"  ~ 0.08
    )) |>
  write_csv("data/ks_pres_tens_trees.csv")
  paste("data/ks_pres_tens_trees.csv")
}

#' @para data  e.g., ks_trees_csv
clean_ks_trees_err <- function(data) {
#  d <- read_csv("data/ks_pres_tens_trees.csv") |>
  d <- read_csv(data) |>
    group_by(species, tree_id, pres_type, pressure) |>
    summarize(ks_mean = mean(ks, na.rm = TRUE),
     ks_sd = sd(ks, na.rm = TRUE)) |>
    ungroup()

    d_tens <- d |>
      filter(pres_type == "tens") |>
      rename(tens_calib_mean = "ks_mean") |>
      rename(tens_calib_sd = "ks_sd") |>
      dplyr::select(!pres_type)

    d_pres <- d |>
      filter(pres_type == "pres") |>
      rename(pres_calib_mean = "ks_mean") |>
      rename(pres_calib_sd = "ks_sd") |>
      dplyr::select(!pres_type)

    full_join(d_pres, d_tens) |>
      write_csv("data/ks_pres_tens_spp_err.csv")
    paste("data/ks_pres_tens_spp_err.csv")
}


clean_cond_count <- function(data, file) {
#  d <- read_csv("data-raw/cond_count.csv")
  d <- read_csv(data)
  d |>
    pivot_longer(c(count_2, count_5, count_8), names_to = "count_p", values_to = "count")  |>
    pivot_longer(c(total_2, total_5, total_8), names_to = "total_p", values_to = "total") |>
    filter(
      (count_p == "count_2" & total_p == "total_2") |
      (count_p == "count_5" & total_p == "total_5") |
      (count_p == "count_8" & total_p == "total_8")) |>
    mutate(pressure = case_when(
      count_p == "count_2" ~ 0.02,
      count_p == "count_5" ~ 0.05,
      count_p == "count_8" ~ 0.08,
    )) |>
    dplyr::select(species, count, total, pressure) |>
    write_csv(file)
  paste(file)
}

rubber_impute <- function(csv, y2015 = TRUE) {
  d <- read_csv(csv) |>
    janitor::clean_names() |>
    mutate(date = mdy_hm(date))
  d2015 <- d |>
    filter(str_detect(date, "2015")) |>
    dplyr::select(-t16_0_0)
  d2016 <- d |>
    filter(str_detect(date, "2016"))
  if (y2015) {
    d2015 |>
      kNN()
  } else {
    d2016 |>
      kNN()
  }
}

rubber_impute <- function(csv, y2015 = TRUE) {
  d <- read_csv(csv) |>
    janitor::clean_names() |>
    mutate(date = mdy_hm(date))
  d2015 <- d |>
    filter(str_detect(date, "2015")) |>
    dplyr::select(-t16_0_0)
  d2016 <- d |>
    filter(str_detect(date, "2016"))
  if (y2015) {
    d2015 |>
      kNN()
  } else {
    d2016 |>
      kNN()
  }
}

missForest_all <- function(csv) {
  d <- read_csv(csv) |>
    janitor::clean_names() |>
    mutate(date = mdy_hm(date)) |>
    dplyr::filter(date <= ymd("2016-01-31")) |>
    dplyr::filter(date >= ymd("2016-01-01")) |>
    dplyr::select(-t11_0_2, -t16_0_0) |>
    dplyr::select(-date)

  d |>
    as.data.frame() |>
    missForest(parallelize = "variables")
}

# library(targets)
# library(tidyverse)
# tar_read(imputed_long_2015)

missForest_long <- function(csv, year = 2015, month = 1) {
  d <- read_csv(csv) |>
    janitor::clean_names() |>
    mutate(date = mdy_hm(date)) |>
    mutate(year = year(date)) |>
    mutate(month = month(date)) |>
    filter(month == {{month}}) |>
    mutate(day = day(date)) |>
    mutate(yday = yday(date)) |>
    mutate(time = hour(date) * 60 + minute(date)) |>
    mutate(cos_transformed_day = cos((yday - 1) / 365 * 2 * pi)) |>
    mutate(cos_transformed_time = cos((time / 1440) * 2 * pi)) |>
    dplyr::select(year, cos_transformed_day, cos_transformed_time,
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
    as.data.frame() |>
    # as_tibble()
    missForest::missForest(parallelize = "forests")
}

missForest_clean <- function(csv, year = 2015, month = 1) {
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
    dplyr::select(year, yday, cos_transformed_day, cos_transformed_time,
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

backtransform_date <- function(csv, year = 2015, month = 1, imputed_df) {
  df_time <- read_csv(csv) |>
    janitor::clean_names() |>
    mutate(date = mdy_hm(date)) |>
    mutate(year = year(date)) |>
    mutate(month = month(date)) |>
    filter(month %in% {{month}}) |>
    mutate(day = day(date)) |>
    mutate(yday = yday(date)) |>
    mutate(time = hour(date) * 60 + minute(date)) |>
    dplyr::select(year, yday, time, t01_0_0:t15_0_0) |>
    pivot_longer(c(t01_0_0:t15_0_0), names_to = "id", values_to = "ks") |>
    filter(year == {{year}})

  imputed_df |>
    rename(yday = cos_transformed_day) |>
    rename(time = cos_transformed_time) |>
    mutate(yday = df_time$yday) |>
    mutate(time = df_time$time)

}

missForest_each <- function(csv, tree) {
  d <- read_csv(csv) |>
    janitor::clean_names()
  vpd <- d |>
    dplyr::select(vpd)
  missing_df <- d |>
    dplyr::select(matches(tree))
  bind_cols(vpd, missing_df) |>
    as.data.frame() |>
    missForest::missForest()
}

clean_imputed_df <- function(imputed_rest) {
  imputed_rest$ximp |>
    as_tibble() |>
    dplyr::select(-yday) |>
    filter(year == 2016)
}

make_long_nonimputed_df <- function(csv) {
  read_csv(csv) |>
  janitor::clean_names() |>
  mutate(date = mdy_hm(date)) |>
  mutate(year = year(date)) |>
  mutate(yday = yday(date)) |>
  mutate(time = hour(date) * 60 + minute(date)) |>
  mutate(h = time %/% 60) |>
  mutate(m = time %% 60) |>
  mutate(time = sprintf("%02d:%02d:%02d", h, m, 0)) |>
  mutate(date = ymd(paste(year, "01", "01", sep= "-")) + days(yday - 1)) |>
  dplyr::select(year, date, time,
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
  dplyr::select(-id) |>
  arrange(date) |>
  arrange(dir) |>
  arrange(dep) |>
  arrange(tree)
}

missForest_comb <- function(csv, combined_imputed_mapped) {
  d <- read_csv(csv) |>
    janitor::clean_names()
  vpd <- d |>
    dplyr::select(vpd)
  missing_df <- d |>
    dplyr::select(matches("1[2-6]"))
  imputed_df <- combined_imputed_mapped |>
    dplyr::select(matches("0_0"))
  bind_cols(vpd, missing_df, imputed_df) |>
    as.data.frame() |>
    missForest(parallelize = "variables")
}

generate_imputed_df <- function(csv, impute_data_full, combined_imputed_mapped) {
  date_df <- read_csv(here(csv)) |>
    janitor::clean_names() |>
    mutate(date = mdy_hm(date)) |>
    dplyr::select(date:par)

  impute_data_full_df <- impute_data_full$ximp |> as_tibble() |>
    dplyr::select(-vpd)

  names_ <- c(names(impute_data_full_df), names(combined_imputed_mapped))
  dup_names <- names_[duplicated(names_)]

  imputed_cleaned_df  <- impute_data_full_df |>
    dplyr::select(-dup_names)

  bind_cols(date_df, imputed_cleaned_df, combined_imputed_mapped) %>%
    dplyr::select(date:vpd, sort(names(.[7:ncol(.)])))

}

generate_dbh_sap_stan_data <- function(sapwood_depth_csv) {
  d <- read_csv(sapwood_depth_csv)
  list(
    y = log(d$sapwood_depth),
    x = log(d$dbh),
    N = nrow(d)
  )
}

generate_full_date_dbh <- function(girth_increment_csv, initial_dbh_csv) {

}

generate_dir_dep_stan_data <- function(imputed_full_df, time_res = c("daily", "hourly", "10mins"), fct = c("dir", "dep")) {
  # Create a reference data frame
  ref_df <- imputed_full_df |>
    filter(dir == "S" & dep == 2) |>
    select(year, date, time, tree, ks) |>
    filter(ks > 1e-6) |>
    rename(ks_ref = ks)

  if (time_res == "daily") {
    imp_df <- imputed_full_df |>
      left_join(ref_df, by = c("year", "date", "time", "tree")) |>
      group_by(date, tree, dir, dep) |>
      summarise(ks = mean(ks, na.rm = TRUE),
                ks_ref = mean(ks_ref, na.rm = TRUE),
                .groups = "drop") |>
      filter(dir != "S" | dep != 2) |>
      mutate(dir_fct = factor(dir, levels = c("S", "N", "E", "W"))) |>
      mutate(dep_fct = as.factor(dep)) |>
      mutate(time2 = date |> as.factor() |> as.numeric()) |>
      filter(ks > 1e-6) |>
      filter(!is.na(ks_ref))
  } else if (time_res == "hourly") {
    imp_df <- imputed_full_df |>
      left_join(ref_df, by = c("year", "date", "time", "tree")) |>
      mutate(hour = as.integer(substr(time, 1, 2))) |>
      group_by(year, date, hour, tree, dir, dep) |>
      summarise(ks = mean(ks, na.rm = TRUE),
                ks_ref = mean(ks_ref, na.rm = TRUE),
                .groups = "drop") |>
      filter(dir != "S" | dep != 2) |>
      mutate(dir_fct = factor(dir, levels = c("S", "N", "E", "W"))) |>
      mutate(dep_fct = as.factor(dep)) |>
      mutate(time2 = paste(date, hour) |> as.factor() |> as.numeric()) |>
      filter(ks > 1e-6) |>
      filter(!is.na(ks_ref))
  } else {
    imp_df <- imputed_full_df |>
      left_join(ref_df, by = c("year", "date", "time", "tree")) |>
      filter(dir != "S" | dep != 2) |>
      mutate(dir_fct = factor(dir, levels = c("S", "N", "E", "W"))) |>
      mutate(dep_fct = as.factor(dep)) |>
      mutate(time2 = paste(date, time) |> as.factor() |> as.numeric()) |>
      filter(ks > 1e-6) |>
      filter(!is.na(ks_ref))
  }

  if (fct == "dep") {
    imp_df <- imp_df |>
      filter(dir == "S")
    xd <- model.matrix(log(ks) ~  dep_fct, data = imp_df)
  } else {
    imp_df <- imp_df |>
      filter(dep == 2)
    xd <- model.matrix(log(ks) ~  dir_fct, data = imp_df)
  }

  # no intercept
  xd <- xd[, -1]

  list(
    N = nrow(imp_df),
    K = ncol(xd),
    M = imp_df |> pull(tree) |> unique() |> length(),
    log_ks = log(imp_df$ks),
    log_ks_ref = log(imp_df$ks_ref),
    tree = imp_df |> pull(tree) |> as.character() |> as.factor() |> as.numeric(),
    x = xd
  )
}

generate_dir_dep_imp_data <- function(imputed_full_df, post_dir, post_dep) {

  post_dir <- post_dir |>
    janitor::clean_names()
  post_dep <- post_dep |>
    janitor::clean_names()

  beta_dir_post <- post_dir |>
    dplyr::select(starts_with("beta"))
  beta_dir_post_2 <- bind_cols(1, beta_dir_post)
  colnames(beta_dir_post_2) <- c("S", "N", "E", "W")

  beta_dep_post <- post_dep |>
    dplyr::select(starts_with("beta"))
  beta_dep_post2 <- bind_cols(1, beta_dep_post)
  colnames(beta_dep_post2) <- c(2, 4, 6)

  beta_dir_post_4 <- beta_dir_post_2 + beta_dep_post2 |> pull(`4`)
  beta_dir_post_4 <- as_tibble(beta_dir_post_4)
  beta_dir_post_6 <- beta_dir_post_2 + beta_dep_post2 |> pull(`6`)
  beta_dir_post_6 <- as_tibble(beta_dir_post_6)

  beta_df <- bind_rows(
      beta_dir_post_2 |>
        mutate(dep = 2),
      beta_dir_post_4 |>
        mutate(dep = 4),
      beta_dir_post_6 |>
        mutate(dep = 6)
    ) |>
    pivot_longer(-dep, names_to = "dir") |>
    mutate(dir = factor(dir, levels = c("S", "N", "E", "W"))) |>
    group_by(dep, dir) |>
    nest() |>
    rename(beta = data)

  date_vec <- imputed_full_df |>
    select(date) |>
    distinct() |>
    pull()

  time_vec <- imputed_full_df |>
    select(time) |>
    distinct() |>
    pull()

  tmp <- expand_grid(date = date_vec, time = time_vec,
    dir = factor(c("S", "N", "E", "W"), levels = c("S", "N", "E", "W")),
    dep = c(2, 4, 6), tree = sprintf("t%02d", 1:15)) |>
    mutate(year = str_sub(date, 1, 4) |> as.numeric())

  ref_df <- imputed_full_df |>
    filter(dir == "S" & dep == 2) |>
    select(year, date, time, tree, ks) |>
    rename(ks_ref = ks)

  imp_df <- imputed_full_df |>
    full_join(ref_df, by = c("year", "date", "time", "tree")) |>
    mutate(dir_fct = factor(dir, levels = c("S", "N", "E", "W"))) |>
    mutate(dep_fct = as.factor(dep)) |>
    mutate(time2 = paste(date, time) |> as.factor() |> as.numeric())  |>
    dplyr::select(year, date, time, ks, ks_ref, tree, dir, dep) |>
    mutate(tree = as.character(tree))

  imp_long <- full_join(tmp, imp_df, by = c("year", "date", "time", "tree", "dir", "dep"))

  imp_nd <- imp_long |>
    group_by(dep, dir) |>
    nest()

  imp_nd2 <- full_join(imp_nd, beta_df)

  imp_nd3 <- imp_nd2 |>
    mutate(beta_mid = map_dbl(beta, \(x)unlist(x) |> median())) |>
    mutate(beta_var = map_dbl(beta, var)) |>
    mutate(pred_sd = map(data, \(x) sqrt(log(x$ks_ref)^2 * beta_var))) |>
    mutate(pred_m = map2(data, pred_sd, \(x, y)log(x$ks_ref) + beta_mid)) |>
    mutate(pred_ll = map2(data, pred_sd, \(x, y)log(x$ks_ref) + beta_mid - 1.96 * beta_var)) |>
    mutate(pred_uu = map2(data, pred_sd, \(x, y)log(x$ks_ref) + beta_mid + 1.96 * beta_var))

  imp_new <- imp_nd3 |>
    dplyr::select(dir, dep, data, pred_m, pred_ll, pred_uu) |>
    unnest(cols = c(data, pred_m, pred_ll, pred_uu)) |>
    ungroup() |>
    mutate(across(c(pred_m, pred_ll, pred_uu),
                  ~ case_when(is.infinite(.x) ~ NA,
                              # is.na(.x) ~ 0,
                              TRUE ~ .x)))
  imp_new
  # imp_nd3
}
