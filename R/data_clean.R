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
  # d_pres <- d2 |>
  #   filter(pres_type != "tens") |>
  #   rename(pres_calib = "ks") |>
  #   dplyr::select(!pres_type)

  # d_tens <- d2 |>
  #   filter(pres_type == "tens") |>
  #   rename(tens_calib = "ks") |>
  #   dplyr::select(!pres_type)

  # full_join(d_pres, d_tens) |>
  #   mutate(pressure = case_when(
  #     pressure == "p_2" ~ 0.02,
  #     pressure == "p_5"  ~ 0.05,
  #     pressure == "p_8"  ~ 0.08
  #   )) |>
  # write_csv("data/ks_pres_tens_trees.csv")
  # paste("data/ks_pres_tens_trees.csv")
}
