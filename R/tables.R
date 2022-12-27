library(tidyverse)
library(targets)

tar_load(xylem_lab)

segments <- tar_read(fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.08)

pool <- tar_read(fit_ab_summary_granier_without_traits_full_pool_sap_all_clean_0.08)

tmp <- tar_read(fit_ab_each_sap_sp_clean_0.08)

segments_sep <- map_dfr(tmp$fit_segments,
  \(x)x$summary |> filter(variable %in% c("log_a", "b", "a")))  |>
  mutate(species = rep(tmp$species, each = 3))

segments_sep |>
  filter(variable == "a") |>
  mutate_if(is.numeric, \(x) round(x, digits = 2)) |>
  mutate_if(is.numeric, \(x) format(x, nsmall = 2, trim = TRUE)) |>
  mutate(multilevel_sep_a = str_c(q50, " [", q2.5, ", ", q97.5, "]")) |>
  dplyr::select(-variable:-ess_tail)

clean_sep <- function(data, var, name) {
  data <- data |>
    filter(variable == var) |>
    mutate_if(is.numeric, \(x) round(x, digits = 2)) |>
    mutate_if(is.numeric, \(x) format(x, nsmall = 2, trim = TRUE)) |>
    mutate(hoge = str_c(q50, " [", q2.5, ", ", q97.5, "]")) |>
    dplyr::select(-variable:-ess_tail)
  colnames(data)[2] <- name
  data
}

clean_all <- function(data, var, name) {
   data <- data |>
    filter(str_detect(variable, "^a")) |>
    full_join(names_data) |>
    dplyr::select(species, xylem_long_fct, variable, q50, q2.5, q97.5) |>
    mutate_if(is.numeric, \(x) round(x, digits = 2)) |>
    mutate_if(is.numeric, \(x) format(x, nsmall = 2, trim = TRUE)) |>
    mutate(multilevel_a = str_c(q50, " [", q2.5, ", ", q97.5, "]")) |>
    arrange(xylem_long_fct, species)

   if (var == "ln_a") {
    data <- data |>
      filter(str_detect(variable, "alpha\\[1"))
   } else if (var == "b") {
    data <- data |>
      filter(str_detect(variable, "alpha\\[2"))
   } else if (var == "a") {
    data <- data |>
      filter(str_detect(variable, "^a\\["))
   }
   colnames(data)[ncol(data)] <- name
   data |>
    dplyr::select(-variable:-q97.5)
}


clean_all(segments, "a", "multi_full_a") |>
  full_join(clean_all(segments, "ln_a", "multi_full_ln_a")) |>
  full_join(clean_all(segments, "b", "multi_full_b")) |>
  full_join(clean_all(pool, "a", "pool_full_a")) |>
  full_join(clean_all(pool, "ln_a", "pool_full_ln_a")) |>
  full_join(clean_all(pool, "b", "pool_full_b")) |>
  full_join(clean_sep(segments_sep, "a", "multi_sep_a")) |>
  full_join(clean_sep(segments_sep, "log_a", "multi_sep_ln_a")) |>
  full_join(clean_sep(segments_sep, "b", "multi_sep_a")) |>
  full_join(clean_sep(pool_sep, "a", "pool_sep_a")) |>
  full_join(clean_sep(pool_sep, "log_a", "pool_sep_ln_a")) |>
  full_join(clean_sep(pool_sep, "b", "pool_sep_b")) |>
  DT::datatable()

ab_data <- full_join(segments_a, segments_b) |>
  full_join(pool_a) |>
  full_join(pool_b) |>
  write_csv("ab_008_without_traits.csv")
