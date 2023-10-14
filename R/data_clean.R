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

generate_dir_dep_stan_data2 <- function(imputed_full_df){
  # Create a reference data frame
  ref_df <- imputed_full_df |>
    filter(dir == "S" & dep == 2) |>
    select(year, date, time, tree, ks) |>
    filter(ks > 1e-6) |>
    rename(ks_ref = ks)

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

  xd <- model.matrix(log(ks) ~  dep_fct + dir_fct, data = imp_df)

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

generate_post_dir <- function(draws_dir) {
  # post_dir <- tar_read(fit_draws_dir_dep_dir)
  post_dir <- draws_dir |>
    janitor::clean_names()

  beta_dir_post <- post_dir |>
    dplyr::select(starts_with("beta"))
  beta_dir_post_2 <- bind_cols(0, beta_dir_post)
  colnames(beta_dir_post_2) <- c("S", "N", "E", "W")
  beta_dir_post_2
}

generate_post_dep <- function(draws_dep) {
  # post_dep <- tar_read(fit_draws_dir_dep_dep)
  post_dep <- draws_dep |>
    janitor::clean_names()

  beta_dep_post <- post_dep |>
    dplyr::select(starts_with("beta"))
  beta_dep_post_2 <- bind_cols(0, beta_dep_post)
  colnames(beta_dep_post_2) <- c("2", "4", "6")
  beta_dep_post_2
}


# tar_load(fit2_draws_dir_dep)
# generate_post_dir_dep(fit2_draws_dir_dep, "dir_only")
# generate_post_dir_dep(fit2_draws_dir_dep, "dep_only")
# generate_post_dir_dep(fit2_draws_dir_dep, "dir_dep")

generate_post_dir_dep <- function(draws_dir_dep, data_type = c("dir_only", "dep_only", "dir_dep")) {
  post <- draws_dir_dep |>
    janitor::clean_names()

  beta_post <- post |>
    dplyr::select(starts_with("beta")) |>
    mutate(s_2 = 0) |>
    rename(
      "4" = beta_1,
      "6" = beta_2,
      "n" = beta_3,
      "e" = beta_4,
      "w" = beta_5)

  if (data_type == "dir_only") {
    beta_post <- beta_post |>
      mutate(across(c(`4`, `6`), ~ median(.))) |>
      mutate(s_4 = `4`) |>
      mutate(s_6 = `6`)
  } else if (data_type == "dep_only") {
    beta_post <- beta_post |>
      mutate(across(c(n, e, w), ~ median(.))) |>
      mutate(s_4 = `4`) |>
      mutate(s_6 = `6`)
  }

   beta_post <- beta_post |>
    mutate(s_4 = `4`,
           s_6 = `6`,
           n_2 = n,
           n_4 = `4` + n,
           n_6 = `6` + n,
           e_2 = e,
           e_4 = `4` + e,
           e_6 = `6` + e,
           w_2 = w,
           w_4 = `4` + w,
           w_6 = `6` + w) |>
    select(s_2, s_4, s_6, n_2, n_4, n_6, e_2, e_4, e_6, w_2, w_4, w_6)
  beta_post
}

# 15: Hevea brasiliensis
generate_post_ab <- function(draws) {
  draws |>
    janitor::clean_names() |>
    dplyr::select(alpha_1_15, alpha_2_15, sigma) |>
    rename(log_a = alpha_1_15, b = alpha_2_15)
}

generate_post_ab_each <- function(draws, pool = TRUE) {
  s1 <- draws |>
    filter(species == "Hevea brasiliensis")
  if (pool) {
    s1 <- s1$fit_pool[[1]]$draws
  } else {
    s1 <- s1$fit_segments[[1]]$draws
  }
   s1 |>
    posterior::as_draws_df() |>
    dplyr::select(log_a, b)
}

generate_dir_dep_imp_data <- function(imputed_full_df) {

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
    full_join(tmp, by = c("year", "date", "time", "tree", "dir", "dep")) |>
    full_join(ref_df, by = c("year", "date", "time", "tree")) |>
    mutate(dir_fct = factor(dir, levels = c("S", "N", "E", "W"))) |>
    mutate(dep_fct = as.factor(dep)) |>
    mutate(time2 = paste(date, time) |> as.factor() |> as.numeric())  |>
    dplyr::select(year, date, time, ks, ks_ref, tree, dir, dep) |>
    mutate(tree = as.character(tree))

  imp_df
}

pred_sap <- function(x, para) {
  log_mu <- para$alpha + para$beta * log(x)
  exp(log_mu)
}

calc_s <- function(dbh, depth, bark = 0.77) {
  r <- dbh / 2 - bark
  b <- depth - 2
  pi * (2 * r - b - depth) * 2
}

calc_s_plus <- function(dbh, sap, bark = 0.77) {
  r1 <- dbh / 2 - 6 - bark
  r2 <- dbh / 2 - sap - bark
  pi * (r1^2 - r2^2)
}

calc_s_cut <- function(dbh, sap, bark = 0.77) {
  r1 <- dbh / 2 - 4 - bark
  r2 <- dbh / 2 - sap - bark
  pi * (r1^2 - r2^2)
}

generate_dbh_imp_data2 <- function(girth_increment_csv, initial_dbh_csv) {
  d1 <- read_csv(girth_increment_csv)  |>
    janitor::clean_names() |>
    mutate(across(where(is.numeric), \(x) x * 0.1)) # mm to cm
  d2 <- read_csv(initial_dbh_csv) |>
    rename(tree = tree_id) |>
    mutate(tree = sprintf("t%02d", tree))
  d1[is.na(d1)] <- 0

  d1_re <- data.frame(date = c("24/12/2014"), matrix(rep(0, 16), nrow = 1, ncol = 16))  |>
    janitor::clean_names() |>
    bind_rows(d1)

  names(d1_re)[-1] <- 1:16

# transform the girth increment data frame to a long format of dbh increment
  d1_long <- d1_re |>
    rename_with(~ paste0("tree_", .x), -date) |>
    pivot_longer(-date, names_to = "tree", values_to = "girth_increment") |>
    mutate(
      date = lubridate::dmy(date),
      tree = as.numeric(str_replace(tree, "tree_", ""))) |>
    mutate(tree = sprintf("t%02d", tree))

# Join the transformed data frames by tree
  dbh_data <- d1_long |>
    left_join(d2, by = "tree") |>
    group_by(tree) |>
    arrange(date) |>
    mutate(girth = pi * initial_dbh + girth_increment) |>
    mutate(dbh = girth / pi) |>
    ungroup() |>
    arrange(tree, date)  # Arrange data by tree and date

  # dbh_data |>
  #   mutate(initial_dbh = round(initial_dbh, 2)) |>
  #   mutate(dbh = round(dbh, 2)) |>
  #   mutate(girth = round(girth, 2))  |>
  #   DT::datatable)

# Now create a new data frame where each tree has a row for each date
  all_dates <- seq(min(dbh_data$date), max(dbh_data$date), by = "day")

  dbh_data_interpolated <- dbh_data %>%
    group_by(tree) %>%
    do(data.frame(date = all_dates,
                  dbh = approx(x = .$date, y = .$dbh, xout = all_dates)$y,
                  stringsAsFactors = FALSE)) %>%
    ungroup() |>
    filter(date < "2017-03-01")

  dbh_data_interpolated
}


generate_dbh_imp_data <- function(girth_increment_csv, initial_dbh_csv) {
  d1 <- read_csv(girth_increment_csv)  |>
    janitor::clean_names() |>
    mutate(across(where(is.numeric), \(x) x * 0.1)) # mm to cm
  d2 <- read_csv(initial_dbh_csv) |>
    rename(tree = tree_id) |>
    mutate(tree = sprintf("t%02d", tree))
  d1[is.na(d1)] <- 0

  d1_re <- data.frame(date = c("24/12/2014"), matrix(rep(0, 16), nrow = 1, ncol = 16))  |>
    janitor::clean_names() |>
    bind_rows(d1)

  names(d1_re)[-1] <- 1:16

# transform the girth increment data frame to a long format of dbh increment
  d1_long <- d1_re |>
    rename_with(~ paste0("tree_", .x), -date) |>
    pivot_longer(-date, names_to = "tree", values_to = "girth_increment") |>
    mutate(
      date = lubridate::dmy(date),
      tree = as.numeric(str_replace(tree, "tree_", ""))) |>
    mutate(tree = sprintf("t%02d", tree))

# Join the transformed data frames by tree
  dbh_data <- d1_long |>
    left_join(d2, by = "tree") |>
    group_by(tree) |>
    arrange(date) |>
    mutate(girth = pi * initial_dbh + girth_increment) |>
    mutate(dbh = girth / pi) |>
    ungroup() |>
    arrange(tree, date)  # Arrange data by tree and date

  # dbh_data |>
  #   mutate(initial_dbh = round(initial_dbh, 2)) |>
  #   mutate(dbh = round(dbh, 2)) |>
  #   mutate(girth = round(girth, 2))  |>
  #   DT::datatable)

# Now create a new data frame where each tree has a row for each date
  all_dates <- seq(min(dbh_data$date), max(dbh_data$date), by = "day")

  dbh_data_interpolated <- dbh_data %>%
    group_by(tree) %>%
    do(data.frame(date = all_dates,
                  dbh = approx(x = .$date, y = .$dbh, xout = all_dates)$y,
                  stringsAsFactors = FALSE)) %>%
    ungroup() |>
    filter(date < "2017-01-01") |>
    filter(date >= "2015-01-01")

  dbh_data_interpolated
}

calc_fd <- function(log_ks, post) {
  mu <- post$log_a + post$b * log_ks
  exp(mu)
}

calc_fd_granier_a <- function(log_ks, post) {
  mu <- log(119) + post$b * log_ks
  exp(mu)
}

calc_fd_granier_b <- function(log_ks, post) {
  mu <- post$log_a + 1.23 * log_ks
  exp(mu)
}

#' @title Modified createFolds function
create_single_fold <- function(data, k, i) {
  fold_sizes <- floor(nrow(data) / k)
  folds <- split(data, rep(1:k, each = fold_sizes, length.out = nrow(data)))
  tmp <- folds[[i]]
  tmp |>
    dplyr::select(colnames(data))
}

calc_quantiles <- function(x, na.rm = TRUE) {
  q <- quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = na.rm)
  return(list(ll = q[1], l = q[2], m = q[3], h = q[4], hh = q[5]))
}

generate_sarea_df <- function(dbh_imp_df, post_slen_m) {
  dbh_df <- dbh_imp_df |>
    mutate(slen = pred_sap(dbh, post_slen_m))  |>
    mutate(s_0_2 = calc_s(dbh, 2)) |>
    mutate(s_2_4 = calc_s(dbh, 4)) |>
    mutate(s_4_6 = ifelse(slen < 6, calc_s_cut(dbh, slen), calc_s(dbh, 6))) |>
    mutate(s_6_c = ifelse(slen >= 6, calc_s_plus(dbh, slen), 0))
  dbh_df
}

generate_ab_uncertainty <- function(dir_dep_imp_df, dbh_imp_df,
  post_ab_pool_mc, post_ab_segments_mc, post_slen, post_dir_dep,
  k = 50, i = 1) {

  post_slen_m  <- post_slen |>
    summarize(across(everything(), median))

  post_dir_dep_m <- post_dir_dep |>
    mutate(dep_dir_mid = map_dbl(beta, median)) |>
    dplyr::select(-beta) |>
    unnest(cols = c()) |>
    ungroup()

  # add folds here
  dir_dep_imp_df <- create_single_fold(dir_dep_imp_df, k, i)
  tmp <- full_join(dir_dep_imp_df, post_dir_dep_m, by = c("dir", "dep")) |>
    mutate(log_ks = ifelse(is.na(ks), log(ks_ref) + dep_dir_mid, log(ks))) |>
    mutate(log_ks = ifelse(is.infinite(log_ks), NA, log_ks))

  dbh_df <- dbh_imp_df |>
    mutate(slen = pred_sap(dbh, post_slen_m))  |>
    mutate(s_0_2 = calc_s(dbh, 2)) |>
    mutate(s_2_4 = calc_s(dbh, 4)) |>
    mutate(s_4_6 = ifelse(slen < 6, calc_s_cut(dbh, slen), calc_s(dbh, 6))) |>
    mutate(s_6_c = ifelse(slen >= 6, calc_s_plus(dbh, slen), 0))

  tmp2 <- full_join(tmp, dbh_df, by = c("date", "tree"))

  s_df <- tmp2 |>
    mutate(s = case_when(
      dep == 2 ~ s_0_2,
      dep == 4 ~ s_2_4,
      dep == 6 ~ s_4_6
    )) |>
    dplyr::select(-s_0_2, -s_2_4, -s_4_6, -s_6_c)

  # deeper than 6cm
  s6_df <- tmp2 |>
    filter(dep == 6) |>
    filter(s_6_c > 0) |>
    mutate(dep = 7) |>
    mutate(s = s_6_c) |>
    dplyr::select(-s_0_2, -s_2_4, -s_4_6, -s_6_c) |>
    mutate(log_ks = log_ks - log(2))

  s_df <- bind_rows(s_df, s6_df)

  rm(tmp)
  rm(tmp2)
  rm(dbh_df)
  rm(dir_dep_imp_df)
  rm(dbh_imp_df)
  rm(s6_df)

  s_df |>
    setDT() |>
    mutate(post_ab_pool = list(post_ab_pool_mc)) |>
    mutate(post_ab_segments = list(post_ab_segments_mc)) |>
    mutate(fd_10min_pool = map2(log_ks, post_ab_pool, calc_fd)) |>
    mutate(fd_10min_segments = map2(log_ks, post_ab_segments, calc_fd)) |>
    mutate(fd_pool = map(fd_10min_pool, calc_quantiles)) |>
    mutate(fd_segments = map(fd_10min_segments, calc_quantiles)) |>
    mutate(fd_granier = 119 * exp(log_ks)^1.23) |>
    dplyr::select(-post_ab_pool, -post_ab_segments, -fd_10min_pool, -fd_10min_segments) |>
    unnest_wider(fd_pool, names_sep = "_") |>
    unnest_wider(fd_segments, names_sep = "_") |>
    ab_scaling()
}

generate_ab_uncertainty_granier <- function(dir_dep_imp_df, dbh_imp_df,
  post_ab_pool_mc, post_ab_segments_mc, post_slen, post_dir_dep,
  k = 50, i = 1) {

  post_slen_m  <- post_slen |>
    summarize(across(everything(), median))

  post_dir_dep_m <- post_dir_dep |>
    mutate(dep_dir_mid = map_dbl(beta, median)) |>
    dplyr::select(-beta) |>
    unnest(cols = c()) |>
    ungroup()

  # add folds here
  dir_dep_imp_df <- create_single_fold(dir_dep_imp_df, k, i)
  tmp <- full_join(dir_dep_imp_df, post_dir_dep_m, by = c("dir", "dep")) |>
    mutate(log_ks = ifelse(is.na(ks), log(ks_ref) + dep_dir_mid, log(ks))) |>
    mutate(log_ks = ifelse(is.infinite(log_ks), NA, log_ks))

  dbh_df <- dbh_imp_df |>
    mutate(slen = pred_sap(dbh, post_slen_m))  |>
    mutate(s_0_2 = calc_s(dbh, 2)) |>
    mutate(s_2_4 = calc_s(dbh, 4)) |>
    mutate(s_4_6 = ifelse(slen < 6, calc_s_cut(dbh, slen), calc_s(dbh, 6))) |>
    mutate(s_6_c = ifelse(slen >= 6, calc_s_plus(dbh, slen), 0))

  tmp2 <- full_join(tmp, dbh_df, by = c("date", "tree"))

  s_df <- tmp2 |>
    mutate(s = case_when(
      dep == 2 ~ s_0_2,
      dep == 4 ~ s_2_4,
      dep == 6 ~ s_4_6
    )) |>
    dplyr::select(-s_0_2, -s_2_4, -s_4_6, -s_6_c)

  # deeper than 6cm
  s6_df <- tmp2 |>
    filter(dep == 6) |>
    filter(s_6_c > 0) |>
    mutate(dep = 7) |>
    mutate(s = s_6_c) |>
    dplyr::select(-s_0_2, -s_2_4, -s_4_6, -s_6_c) |>
    mutate(log_ks = log_ks - log(2))

  s_df <- bind_rows(s_df, s6_df)

  rm(tmp)
  rm(tmp2)
  rm(dbh_df)
  rm(dir_dep_imp_df)
  rm(dbh_imp_df)
  rm(s6_df)

  s_df |>
    setDT() |>
    mutate(post_ab_pool = list(post_ab_pool_mc)) |>
    mutate(post_ab_segments = list(post_ab_segments_mc)) |>
    mutate(fd_10min_pool_a = map2(log_ks, post_ab_pool, calc_fd_granier_a)) |>
    mutate(fd_10min_segments_a = map2(log_ks, post_ab_segments, calc_fd_granier_a)) |>
    mutate(fd_10min_pool_b = map2(log_ks, post_ab_pool, calc_fd_granier_b)) |>
    mutate(fd_10min_segments_b = map2(log_ks, post_ab_segments, calc_fd_granier_b)) |>
    mutate(fd_pool_a = map(fd_10min_pool_a, calc_quantiles)) |>
    mutate(fd_segments_a = map(fd_10min_segments_a, calc_quantiles)) |>
    mutate(fd_pool_b = map(fd_10min_pool_b, calc_quantiles)) |>
    mutate(fd_segments_b = map(fd_10min_segments_b, calc_quantiles)) |>
    dplyr::select(-post_ab_pool, -post_ab_segments, -fd_10min_pool_a, -fd_10min_segments_a, -fd_10min_pool_b, -fd_10min_segments_b) |>
    unnest_wider(fd_pool_a, names_sep = "_") |>
    unnest_wider(fd_segments_a, names_sep = "_") |>
    unnest_wider(fd_pool_b, names_sep = "_") |>
    unnest_wider(fd_segments_b, names_sep = "_") |>
    ab_scaling()
}

# dividing by 4 is to get the quarter of the area (each direction)
# 600 is the number of seconds in 10 min
# dividing by 10000 to convert from cm2 to m2
ab_scaling <- function(ab_uncertainty_full_df) {
  ab_uncertainty_full_df |>
    mutate(
      across(
        .cols = starts_with("fd"),
        .fns = ~ .x * s / 4 * 600 / 10000,
        .names = "s_10m_{.col}"
      )
    ) |>
  group_by(tree) |>
  summarize(
    across(
      .cols = starts_with("s_10m_fd_"),
      .fns = ~ sum(.x, na.rm = TRUE) / 2 * 1e-6,
      .names = "s_total_{.col}"
    )
  ) |>
  rename_with(~ str_replace(., "s_total_s_10m_fd_", "s_total_"), starts_with("s_total_s_10m_fd_"))
}


generate_t16_df <- function(rubber_raw_data_csv) {
  t16_df <- read_csv(rubber_raw_data_csv) |>
    janitor::clean_names() |>
    dplyr::select(date, t16_0_0)

  t16_df_re <- t16_df |>
    rename(ks = t16_0_0) |>
    mutate(dir = factor("S", levels = c("S", "N", "E", "W"))) |>
    mutate(dep = 2) |>
    mutate(date = mdy_hm(date)) |>
    mutate(year = year(date)) |>
    mutate(time = format(date, "%H:%M:%S")) |>
    mutate(date = as_date(date))
  t16_df_re
}


add_t16 <- function(rubber_raw_data_csv, dir_dep_imp_df, post_dir_dep) {
  t16_df <- read_csv(rubber_raw_data_csv) |>
    janitor::clean_names() |>
    dplyr::select(date, t16_0_0)

  t16_df_re <- t16_df |>
    rename(ks_t16 = t16_0_0) |>
    mutate(dir = factor("S", levels = c("S", "N", "E", "W"))) |>
    mutate(dep = 2) |>
    mutate(date = mdy_hm(date)) |>
    mutate(year = year(date)) |>
    mutate(time = format(date, "%H:%M:%S")) |>
    mutate(date = as_date(date))

  post_dir_dep_m <- post_dir_dep |>
    mutate(dep_dir_mid = map_dbl(beta, median)) |>
    dplyr::select(-beta) |>
    unnest(cols = c()) |>
    ungroup()

  tmp <- full_join(dir_dep_imp_df, post_dir_dep_m, by = c("dir", "dep")) |>
    mutate(log_ks = ifelse(is.na(ks), log(ks_ref) + dep_dir_mid, log(ks))) |>
    mutate(log_ks = ifelse(is.infinite(log_ks), NA, log_ks)) |>
    mutate(ks = exp(log_ks))

  mean_ks_df <- tmp |>
    group_by(year, date, time, dir, dep) |>
    summarize(ks = mean(ks, na.rm = TRUE), .groups = "drop")

  t16_df_full <- full_join(mean_ks_df, t16_df_re, by = c("year", "date", "time", "dir", "dep")) |>
    mutate(ks = ifelse(is.na(ks_t16), ks, ks_t16)) |>
    dplyr::select(-ks_t16) |>
    mutate(tree = "t16")

  full_join(dir_dep_imp_df, t16_df_full, by = c("year", "date", "time", "ks", "tree", "dir", "dep"))

}


generate_k_data <- function(dir_dep_imp_full_df, post_dir_dep_mid, sarea_df){
  tmp <- full_join(dir_dep_imp_full_df, post_dir_dep_mid, by = c("dir", "dep")) |>
      mutate(log_ks = ifelse(is.na(ks), log(ks_ref) + dep_dir_mid, log(ks))) |>
      mutate(log_ks = ifelse(is.infinite(log_ks), NA, log_ks))
  tmp2 <- full_join(tmp, sarea_df, by = c("date", "tree"))
  s_df <- tmp2 |>
    mutate(s = case_when(
      dep == 2 ~ s_0_2,
      dep == 4 ~ s_2_4,
      dep == 6 ~ s_4_6
    )) |>
    dplyr::select(-s_0_2, -s_2_4, -s_4_6, -s_6_c)

# deeper than 6cm
  s6_df <- tmp2 |>
    filter(dep == 6) |>
    filter(s_6_c > 0) |>
    mutate(dep = 7) |>
    mutate(s = s_6_c) |>
    dplyr::select(-s_0_2, -s_2_4, -s_4_6, -s_6_c) |>
    mutate(log_ks = log_ks - log(2))

  s_df <- bind_rows(s_df, s6_df)
}

test_fun <- function(csv, year = 2015, month = 1) {
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
      vpd, par, t01_0_0:t16_0_0) |>
    pivot_longer(c(t01_0_0:t16_0_0), names_to = "id", values_to = "ks") |>
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
