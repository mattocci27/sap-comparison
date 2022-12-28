# library(tidyverse)
# library(targets)

# tar_load(xylem_lab)

# segments <- tar_read(fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.08)

# pool <- tar_read(fit_ab_summary_granier_without_traits_full_pool_sap_all_clean_0.08)

# tmp <- tar_read(fit_ab_each_sap_sp_clean_0.08)

clean_sep <- function(data, names_data, var, name) {
  data <- data |>
    filter(variable == var) |>
    mutate_if(is.numeric, \(x) round(x, digits = 2)) |>
    mutate_if(is.numeric, \(x) format(x, nsmall = 2, trim = TRUE)) |>
    mutate(hoge = str_c(q50, " [", q2.5, ", ", q97.5, "]")) |>
    dplyr::select(-variable:-ess_tail)
  colnames(data)[2] <- name
  data
}

clean_all <- function(data, names_data, var, name) {
   data <- data |>
    filter(str_detect(variable, "^a")) |>
    full_join(names_data, by = "variable") |>
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

write_ab_csv2 <- function(summary_segments, summary_pool, summary_sep, xylem_lab, out) {
  segments_sep <- map_dfr(summary_sep$fit_segments,
    \(x)x$summary |> filter(variable %in% c("log_a", "b", "a")))  |>
    mutate(species = rep(summary_sep$species, each = 3))

  pool_sep <- map_dfr(summary_sep$fit_pool,
    \(x)x$summary |> filter(variable %in% c("log_a", "b", "a")))  |>
    mutate(species = rep(summary_sep$species, each = 3))

  names_data <- xylem_lab |>
    mutate(alpha1 = str_replace_all(sp_num1, "_", ",")) |>
    mutate(alpha1 = str_c("alpha[", alpha1, "]")) |>
    mutate(alpha2 = str_replace_all(alpha1, "1,", "2,")) |>
    mutate(a = str_c("a[", sp_num, "]")) |>
    pivot_longer(alpha1:a, values_to = "variable")

  clean_all(summary_segments, names_data, "a", "multi_full_a") |>
    full_join(clean_all(summary_segments, names_data, "ln_a", "multi_full_ln_a")) |>
    full_join(clean_all(summary_segments, names_data, "b", "multi_full_b")) |>
    full_join(clean_all(summary_pool, names_data, "a", "pool_full_a")) |>
    full_join(clean_all(summary_pool, names_data, "ln_a", "pool_full_ln_a")) |>
    full_join(clean_all(summary_pool, names_data, "b", "pool_full_b")) |>
    full_join(clean_sep(segments_sep, names_data, "a", "multi_sep_a")) |>
    full_join(clean_sep(segments_sep, names_data, "log_a", "multi_sep_ln_a")) |>
    full_join(clean_sep(segments_sep, names_data, "b", "multi_sep_a")) |>
    full_join(clean_sep(pool_sep, names_data, "a", "pool_sep_a")) |>
    full_join(clean_sep(pool_sep, names_data, "log_a", "pool_sep_ln_a")) |>
    full_join(clean_sep(pool_sep, names_data, "b", "pool_sep_b")) |>
    my_write_csv(out)
}

# ab_data <- full_join(segments_a, segments_b) |>
#   full_join(pool_a) |>
#   full_join(pool_b) |>
#   write_csv("ab_008_without_traits.csv")
