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

# missForest_long <- function(clean_df, year = 2015, month = 1) {
#    d |>
#     filter(year == {{year}}) |>
#     as.data.frame() |>
#     # as_tibble()
#     missForest::missForest(parallelize = "forests")
# }

clean_for_missForest <- function(csv, year = 2015, month = 1) {
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
    dplyr::select(year, cos_transformed_day, cos_transformed_time,
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
    dplyr::select(year, yday, time, t01_0_0:t16_0_0) |>
    pivot_longer(c(t01_0_0:t16_0_0), names_to = "id", values_to = "k") |>
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
  tmp <- imputed_rest$ximp
  start_row <- nrow(tmp) / 2
  tmp2 <- tmp[(start_row + 1):nrow(tmp), ]
  as_tibble(tmp2) |>
  mutate(year = as.integer(2016)) |>
  dplyr::select(year, everything())
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

generate_dir_dep_stan_data <- function(imputed_full_df){
  # Create a reference data frame
  ref_df <- imputed_full_df |>
    filter(dir == "S" & dep == 2) |>
    select(year, date, time, tree, k) |>
    filter(k > 1e-6) |>
    rename(k_ref = k)

  imp_df <- imputed_full_df |>
    left_join(ref_df, by = c("year", "date", "time", "tree")) |>
    mutate(hour = as.integer(substr(time, 1, 2))) |>
    group_by(year, date, hour, tree, dir, dep) |>
    summarise(k = mean(k, na.rm = TRUE),
              k_ref = mean(k_ref, na.rm = TRUE),
              .groups = "drop") |>
    filter(dir != "S" | dep != 2) |>
    mutate(dir_fct = factor(dir, levels = c("S", "N", "E", "W"))) |>
    mutate(dep_fct = as.factor(dep)) |>
    mutate(time2 = paste(date, hour) |> as.factor() |> as.numeric()) |>
    filter(k > 1e-6) |>
    filter(!is.na(k_ref))

  xd <- model.matrix(log(k) ~  dep_fct + dir_fct, data = imp_df)

  # no intercept
  xd <- xd[, -1]

  list(
    N = nrow(imp_df),
    K = ncol(xd),
    M = imp_df |> pull(tree) |> unique() |> length(),
    log_k = log(imp_df$k),
    log_k_ref = log(imp_df$k_ref),
    tree = imp_df |> pull(tree) |> as.character() |> as.factor() |> as.numeric(),
    x = xd
  )
}

generate_dir_dep_stan_data2 <- function(imputed_full_df, time_res = c("daily", "hourly", "10mins"), fct = c("dir", "dep")) {
  # Create a reference data frame
  ref_df <- imputed_full_df |>
    filter(dir == "S" & dep == 2) |>
    select(year, date, time, tree, k) |>
    filter(k > 1e-6) |>
    rename(k_ref = k)

  if (time_res == "daily") {
    imp_df <- imputed_full_df |>
      left_join(ref_df, by = c("year", "date", "time", "tree")) |>
      group_by(date, tree, dir, dep) |>
      summarise(k = mean(k, na.rm = TRUE),
                k_ref = mean(k_ref, na.rm = TRUE),
                .groups = "drop") |>
      filter(dir != "S" | dep != 2) |>
      mutate(dir_fct = factor(dir, levels = c("S", "N", "E", "W"))) |>
      mutate(dep_fct = as.factor(dep)) |>
      mutate(time2 = date |> as.factor() |> as.numeric()) |>
      filter(k > 1e-6) |>
      filter(!is.na(k_ref))
  } else if (time_res == "hourly") {
    imp_df <- imputed_full_df |>
      left_join(ref_df, by = c("year", "date", "time", "tree")) |>
      mutate(hour = as.integer(substr(time, 1, 2))) |>
      group_by(year, date, hour, tree, dir, dep) |>
      summarise(k = mean(k, na.rm = TRUE),
                k_ref = mean(k_ref, na.rm = TRUE),
                .groups = "drop") |>
      filter(dir != "S" | dep != 2) |>
      mutate(dir_fct = factor(dir, levels = c("S", "N", "E", "W"))) |>
      mutate(dep_fct = as.factor(dep)) |>
      mutate(time2 = paste(date, hour) |> as.factor() |> as.numeric()) |>
      filter(k > 1e-6) |>
      filter(!is.na(k_ref))
  } else {
    imp_df <- imputed_full_df |>
      left_join(ref_df, by = c("year", "date", "time", "tree")) |>
      filter(dir != "S" | dep != 2) |>
      mutate(dir_fct = factor(dir, levels = c("S", "N", "E", "W"))) |>
      mutate(dep_fct = as.factor(dep)) |>
      mutate(time2 = paste(date, time) |> as.factor() |> as.numeric()) |>
      filter(k > 1e-6) |>
      filter(!is.na(k_ref))
  }

  if (fct == "dep") {
    imp_df <- imp_df |>
      filter(dir == "S")
    xd <- model.matrix(log(k) ~  dep_fct, data = imp_df)
  } else {
    imp_df <- imp_df |>
      filter(dep == 2)
    xd <- model.matrix(log(k) ~  dir_fct, data = imp_df)
  }

  # no intercept
  xd <- xd[, -1]

  list(
    N = nrow(imp_df),
    K = ncol(xd),
    M = imp_df |> pull(tree) |> unique() |> length(),
    log_k = log(imp_df$k),
    log_k_ref = log(imp_df$k_ref),
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

generate_post_ab_each <- function(draws, species = TRUE) {
  s1 <- draws |>
    filter(species == "Hevea brasiliensis")
  if (species) {
    s1 <- s1$fit_species[[1]]$draws
  } else {
    s1 <- s1$fit_segments[[1]]$draws
  }
   s1 |>
    posterior::as_draws_df() |>
    dplyr::select(log_a, b, sigma)
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
    dep = c(2, 4, 6), tree = sprintf("t%02d", 1:16)) |>
    mutate(year = str_sub(date, 1, 4) |> as.numeric())

  ref_df <- imputed_full_df |>
    filter(dir == "S" & dep == 2) |>
    select(year, date, time, tree, k) |>
    rename(k_ref = k)

  imp_df <- imputed_full_df |>
    full_join(tmp, by = c("year", "date", "time", "tree", "dir", "dep")) |>
    full_join(ref_df, by = c("year", "date", "time", "tree")) |>
    mutate(dir_fct = factor(dir, levels = c("S", "N", "E", "W"))) |>
    mutate(dep_fct = as.factor(dep)) |>
    mutate(time2 = paste(date, time) |> as.factor() |> as.numeric())  |>
    dplyr::select(year, date, time, k, k_ref, tree, dir, dep) |>
    mutate(tree = as.character(tree))

  imp_df
}

pred_sap <- function(x, para) {
  log_mu <- para$alpha + para$beta * log(x)
  exp(log_mu)
}

# sap-area for 0-2, 2-4, and 4-6 cm depth
calc_s <- function(dbh, depth, bark = 0.77) {
  r <- dbh / 2 - bark
  4 * pi * (r + 1 - depth)
}

calc_s_check <- function(dbh, depth, bark = 0.77) {
  r <- dbh / 2 - bark
  b <- depth - 2
  pi * (r - b)^2 - pi * (r - depth)^2
}

# sap-area deeper than 6 cm
# r1: outer
# r2: inner
calc_s_plus <- function(dbh, sap, bark = 0.77) {
  r1 <- dbh / 2 - 6 - bark
  r2 <- dbh / 2 - sap - bark
  pi * (r1^2 - r2^2)
}

# sap-area between 4 - x cm depth
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

calc_fd <- function(log_k, post) {
  mu <- post$log_a + post$b * log_k
  exp(mu)
}

calc_fd_granier_a <- function(log_k, post) {
  mu <- log(119) + post$b * log_k
  exp(mu)
}

calc_fd_granier_b <- function(log_k, post) {
  mu <- post$log_a + 1.23 * log_k
  exp(mu)
}

calc_quantiles <- function(x, na.rm = TRUE) {
  q <- quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = na.rm)
  return(list(ll = q[1], l = q[2], m = q[3], h = q[4], hh = q[5]))
}

generate_sarea_df <- function(dbh_imp_df, post_slen_m) {
  dbh_imp_df |>
    mutate(
      slen = pred_sap(dbh, post_slen_m),
      s_0_2 = calc_s(dbh, 2),
      s_2_4 = calc_s(dbh, 4),
      s_4_6 = ifelse(slen < 6, calc_s_cut(dbh, slen), calc_s(dbh, 6)),
      s_6_c = ifelse(slen >= 6, calc_s_plus(dbh, slen), 0)
    ) |>
    mutate(across(starts_with("s"), ~ .x * 0.25)) # quarter
}

# 600 is the number of seconds in 10 min
# dividing by 10^4 (1e-4) to convert from cm2 to m2
# dividing by 10^3 (1e-3) to convert from g to Kg
fd_scaling <- function(ab_uncertainty_df) {
  ab_uncertainty_df |>
    group_by(tree) |>
    summarize(across(c(fd_m, fd_l, fd_h, fd_mean, fd_sd), mean)) |>
    mutate(
      across(
        .cols = starts_with("fd"),
        .fns = ~ .x * 600 * 1e-4 * 1e-3,
        .names = "scaled_{.col}"
      )
    ) |>
    summarize(across(starts_with("scaled"), sum)) |>
    ungroup() |>
    map_dbl(~ .x / 800) # divide each column by 800
}

# A common function to mutate and join the dir_dep_imp_df with post_dir_dep_mid_df
prepare_dir_dep_imp_df <- function(dir_dep_imp_df, post_dir_dep_df) {
  dir_dep_imp_df |>
    # head(1000) |>
    mutate(dir_dep = str_c(str_to_lower(dir), dep, sep = "_")) |>
    full_join(post_dir_dep_df, by = "dir_dep") |>
    mutate(
      log_k = if_else(is.na(k), log(k_ref) + dir_dep_effects, log(k)),
      log_k = if_else(is.infinite(log_k), NA_real_, log_k)
    )
}

# A common function to summarize data
summarize_df_data <- function(data) {
  data |>
    group_by(year, tree) |>
    summarize(
      fd_m = median(fd),
      fd_l = quantile(fd, 0.025),
      fd_h = quantile(fd, 0.975),
      fd_mean = mean(fd),
      fd_sd = sd(fd),
      .groups = "drop"
    )
}

# A common function to process full_df with sapwood area calculations
process_full_df <- function(dir_dep_imp_df_re, sarea_df2) {
  full_df <- full_join(dir_dep_imp_df_re, sarea_df2, by = c("date", "tree"))

  s_df <- full_df |>
    mutate(s = case_when(
      dep == 2 ~ s_0_2,
      dep == 4 ~ s_2_4,
      dep == 6 ~ s_4_6
    )) |>
    select(-s_0_2, -s_2_4, -s_4_6, -s_6_c)

  s6_df <- full_df |>
    filter(dep == 6, s_6_c > 0) |>
    mutate(dep = 7, s = s_6_c) |>
    select(-s_0_2, -s_2_4, -s_4_6, -s_6_c) |>
    mutate(log_k = log_k - log(2))

  bind_rows(s_df, s6_df)
}

calc_summary <- function(row, s_df) {
  s_df2 <- s_df |>
    mutate(log_fd = row$log_a + row$b * log_k) |>
    mutate(fd = exp(log_fd) * s) # scale by sapwood area

  # Summarizing to get the desired percentiles for log_fd
  s_df3 <- s_df2 |>
    group_by(year, tree) |>
    summarise(fd = sum(fd, na.rm = TRUE), .groups = "drop")

  return(s_df3)
}

generate_ab_uncertainty <- function(post_ab_fit_draws, post_dir_dep_mid, sarea_df, dir_dep_imp_df, n_draws = 1000) {
  post_dir_dep_mid_df <- tibble(dir_dep = names(post_dir_dep_mid),
    dir_dep_effects = as.numeric(post_dir_dep_mid))
  dir_dep_imp_df_re <- prepare_dir_dep_imp_df(dir_dep_imp_df, post_dir_dep_mid_df)
  full_df_processed <- process_full_df(dir_dep_imp_df_re, sarea_df)

  if (is.null(post_ab_fit_draws)) {
    granier_df <- tibble(log_a = log(119), b = 1.23)
    summary_stats <- calc_summary(granier_df, full_df_processed )
  } else {
    # Convert post_ab to a list of lists
    post_ab_list <- lapply(1:n_draws, function(i) as.list(post_ab_fit_draws[i, ]))
    # pmap
    summary_stats <- pmap_dfr(list(post_ab_list, list(full_df_processed)), calc_summary)
  }

  summarize_df_data(summary_stats)

}

dir_dep_map <- function(i, post_dir_dep_draws, dir_dep_imp_df, sarea_df, post_ab_mid) {
  post_dir_dep <- post_dir_dep_draws[i,]
  post_dir_dep_df <- tibble(
    dir_dep = names(post_dir_dep),
    dir_dep_effects = as.numeric(post_dir_dep)
  )

  dir_dep_imp_df_re <- prepare_dir_dep_imp_df(dir_dep_imp_df, post_dir_dep_df)
  full_df_processed <- process_full_df(dir_dep_imp_df_re, sarea_df) |>
    mutate(log_fd = post_ab_mid["log_a"] + post_ab_mid["b"] * log_k) |>
    mutate(fd = exp(log_fd) * s) |> # scale by sapwood area
    group_by(year, tree) |>
    summarise(fd = sum(fd, na.rm = TRUE), .groups = "drop")

  return(full_df_processed)
}

generate_dir_dep_uncertainty <- function(post_ab_mid, post_dir_dep_draws, sarea_df, dir_dep_imp_df, n_draws = 1000) {

  summary_stats <- map_df(1:n_draws, dir_dep_map, post_dir_dep_draws, dir_dep_imp_df, sarea_df, post_ab_mid)

  summarize_df_data(summary_stats)
}

sarea_map <- function(i, dbh_imp_df, post_slen_draws, dir_dep_imp_df_re, post_ab_mid) {
  sarea_df <- generate_sarea_df(dbh_imp_df, post_slen_draws[i, ])
  full_df_processed <- process_full_df(dir_dep_imp_df_re, sarea_df) |>
    mutate(log_fd = post_ab_mid["log_a"] + post_ab_mid["b"] * log_k) |>
    mutate(fd = exp(log_fd) * s) |> # scale by sapwood area
    group_by(year, tree) |>
    summarise(fd = sum(fd, na.rm = TRUE), .groups = "drop")
  return(full_df_processed)
}

generate_sarea_uncertainty <- function(post_ab_mid, post_dir_dep_mid, dbh_imp_df, slen_draws, dir_dep_imp_df, n_draws = 1000) {
  post_dir_dep_mid_df <- tibble(dir_dep = names(post_dir_dep_mid),
    dir_dep_effects = as.numeric(post_dir_dep_mid))

  dir_dep_imp_df_re <- prepare_dir_dep_imp_df(dir_dep_imp_df, post_dir_dep_mid_df)

  summary_stats <- map_df(1:n_draws, sarea_map, dbh_imp_df, slen_draws, dir_dep_imp_df_re, post_ab_mid)

  summarize_df_data(summary_stats)

}
