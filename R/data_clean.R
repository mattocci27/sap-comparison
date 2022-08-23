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
