library(tidyverse)
library(targets)

tar_load(xylem_lab)

segments <- tar_read(fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.08)

pool <- tar_read(fit_ab_summary_granier_without_traits_full_pool_sap_all_clean_0.08)

names_data <- xylem_lab |>
  mutate(alpha1 = str_replace_all(sp_num1, "_", ",")) |>
  mutate(alpha1 = str_c("alpha[", alpha1, "]")) |>
  mutate(alpha2 = str_replace_all(alpha1, "1,", "2,")) |>
  pivot_longer(alpha1:alpha2, values_to = "variable")

segments2 <- segments |>
  filter(str_detect(variable, "alpha")) |>
  full_join(names_data) |>
  dplyr::select(species, xylem_long_fct, q50, q2.5, q97.5, name) |>
  arrange(xylem_long_fct, species)

segments_a <- segments2 |>
  filter(name == "alpha1") |>
  mutate_if(is.numeric, \(x) exp(x)) |>
  mutate_if(is.numeric, \(x) round(x, digits = 2)) |>
  mutate_if(is.numeric, \(x) format(x, nsmall = 2, trim = TRUE)) |>
  mutate(multilevel_a = str_c(q50, " [", q2.5, ", ", q97.5, "]")) |>
  dplyr::select(-q50:-name)

segments_b <- segments2 |>
  filter(name == "alpha2") |>
  mutate_if(is.numeric, \(x) round(x, digits = 2)) |>
  mutate_if(is.numeric, \(x) format(x, nsmall = 2)) |>
  mutate(multilevel_b = str_c(q50, " [", q2.5, ", ", q97.5, "]")) |>
  dplyr::select(-q50:-name)

pool2 <- pool |>
  filter(str_detect(variable, "alpha")) |>
  full_join(names_data) |>
  dplyr::select(species, xylem_long_fct, q50, q2.5, q97.5, name) |>
  arrange(xylem_long_fct, species)

pool_a <- pool2 |>
  filter(name == "alpha1") |>
  mutate_if(is.numeric, \(x) exp(x)) |>
  mutate_if(is.numeric, \(x) round(x, digits = 2)) |>
  mutate_if(is.numeric, \(x) format(x, nsmall = 2, trim = TRUE)) |>
  mutate(pool_a = str_c(q50, " [", q2.5, ", ", q97.5, "]")) |>
  dplyr::select(-q50:-name)

pool_b <- pool2 |>
  filter(name == "alpha2") |>
  mutate_if(is.numeric, \(x) round(x, digits = 2)) |>
  mutate_if(is.numeric, \(x) format(x, nsmall = 2)) |>
  mutate(pool_b = str_c(q50, " [", q2.5, ", ", q97.5, "]")) |>
  dplyr::select(-q50:-name)

ab_data <- full_join(segments_a, segments_b) |>
  full_join(pool_a) |>
  full_join(pool_b) |>
  write_csv("ab_008_without_traits.csv")
