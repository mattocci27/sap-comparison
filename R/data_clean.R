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

