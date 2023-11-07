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
  # data <- summary_segments
  data <- data |>
    filter(str_detect(variable, "^a")) |>
    full_join(names_data, by = "variable") |>
    dplyr::select(species, xylem_long_fct, variable, q50, q2.5, q97.5) |>
    mutate(across(c(q50, q2.5, q97.5),
                  ~ ifelse(str_detect(variable, "alpha\\[1"), exp(.x), .x))) |>
    mutate(across(where(is.numeric), round, digits = 2)) |>
    mutate(across(where(is.numeric),
                  ~ format(.x, nsmall = 2, trim = TRUE))) |>
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

# summary_segments <- tar_read(fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.08)
# summary_pool <- tar_read(fit_ab_summary_granier_without_traits_full_pool_sap_all_clean_0.08)
# summary_sep <- tar_read(fit_ab_each_sap_sp_clean_0.08)
# tar_load(xylem_lab)

write_ab_csv <- function(summary_segments, summary_species, summary_fit_each, xylem_lab, out) {
  seg_inc <- map_dfr(summary_fit_each$fit_segments,
    \(x)x$summary |> filter(variable %in% c("log_a", "b", "a")))  |>
    mutate(species = rep(summary_fit_each$species, each = 3)) |>
    mutate(across(c(q50, q2.5, q97.5),
                  ~ ifelse(str_detect(variable, "log_a"), exp(.x), .x)))

  sp_only <- map_dfr(summary_fit_each$fit_species,
    \(x)x$summary |> filter(variable %in% c("log_a", "b", "a")))  |>
    mutate(species = rep(summary_fit_each$species, each = 3)) |>
    mutate(across(c(q50, q2.5, q97.5),
                  ~ ifelse(str_detect(variable, "log_a"), exp(.x), .x)))

  names_data <- xylem_lab |>
    mutate(alpha1 = str_replace_all(sp_num1, "_", ",")) |>
    mutate(alpha1 = str_c("alpha[", alpha1, "]")) |>
    mutate(alpha2 = str_replace_all(alpha1, "1,", "2,")) |>
    mutate(a = str_c("a[", sp_num, "]")) |>
    pivot_longer(alpha1:a, values_to = "variable")

  clean_all(summary_segments, names_data, "a", "multi_full_a") |>
    full_join(clean_all(summary_segments, names_data, "ln_a", "multi_full_ln_a")) |>
    full_join(clean_all(summary_segments, names_data, "b", "multi_full_b")) |>
    full_join(clean_all(summary_species, names_data, "a", "pool_full_a")) |>
    full_join(clean_all(summary_species, names_data, "ln_a", "pool_full_ln_a")) |>
    full_join(clean_all(summary_species, names_data, "b", "pool_full_b")) |>
    full_join(clean_sep(seg_inc, names_data, "a", "multi_sep_a")) |>
    full_join(clean_sep(seg_inc, names_data, "log_a", "multi_sep_ln_a")) |>
    full_join(clean_sep(seg_inc, names_data, "b", "multi_sep_b")) |>
    full_join(clean_sep(sp_only, names_data, "a", "pool_sep_a")) |>
    full_join(clean_sep(sp_only, names_data, "log_a", "pool_sep_ln_a")) |>
    full_join(clean_sep(sp_only, names_data, "b", "pool_sep_b")) |>
    my_write_csv(out)
}

# ab_data <- full_join(segments_a, segments_b) |>
#   full_join(pool_a) |>
#   full_join(pool_b) |>
#   write_csv("ab_008_without_traits.csv")
# para <- tar_read(ab_var_clean_008)

write_varpart_notrait_table <- function(para, out) {
  tmp <- read_csv(para) |>
    filter(variable %in%
      c("var_a_segment", "var_a_sp", "var_a_xylem",
         "var_b_segment", "var_b_sp", "var_b_xylem")) |>
      mutate_if(is.numeric, \(x) round(x, digits = 2)) |>
      mutate_if(is.numeric, \(x) format(x, nsmall = 2, trim = TRUE)) |>
      mutate(var = str_c(q50, "% [", q25, ", ", q75, "]"))  |>
      dplyr::select(variable, var)

  tmp_a <- tmp  |>
    filter(str_detect(variable, "_a")) |>
    rename(a = var)
  tmp_b <- tmp  |>
    filter(str_detect(variable, "_b")) |>
    rename(b = var)

  tmp_a |>
    mutate(b = tmp_b$b) |>
    mutate(variable = case_when(
      str_detect(variable, "segment") ~ "Segments",
      str_detect(variable, "sp") ~ "Species",
      str_detect(variable, "xylem") ~ "Xylem types",
    )) |>
    rename(Levels = variable) |>
    my_write_csv(out)
}

write_dir_dep_table <- function(dir_dep_summary, out) {
  dir_dep_summary |>
    janitor::clean_names() |>
    filter(str_detect(variable, "beta")) |>
    mutate(variable = case_when(
      variable == "beta[1]" ~ "2-4 cm",
      variable == "beta[2]" ~ "4-6 cm",
      variable == "beta[3]" ~ "N",
      variable == "beta[4]" ~ "E",
      variable == "beta[5]" ~ "W",
    )) |>
    dplyr::select(variable, q2_5, q50, q97_5) |>
    mutate(across(where(is.numeric), exp)) |>
    mutate(across(where(is.numeric), round, 2)) |>
    rename(Factor = variable) |>
    mutate(Posterior = paste0(q50, " [", q2_5, ", ", q97_5, "]")) |>
    dplyr::select(Factor, Posterior) |>
    my_write_csv(out)

  # tmp2 <- dep_summary
  # tmp2 <- tmp2 |>
  #   janitor::clean_names() |>
  #   filter(str_detect(variable, "beta")) |>
  #   mutate(variable = case_when(
  #     variable == "beta[1]" ~ "2-4 cm",
  #     variable == "beta[2]" ~ "4-6 cm",
  #   ))

  # tmp3 <- bind_rows(tmp1, tmp2) |>
  #   dplyr::select(variable, q2_5, q50, q97_5) |>
  #   mutate(across(where(is.numeric), exp))

  # tmp3 |>
  #   mutate(across(where(is.numeric), round, 2)) |>
  #   rename(Factor = variable) |>
  #   mutate(Posterior = paste0(q50, " [", q2_5, ", ", q97_5, "]")) |>
  #   dplyr::select(Factor, Posterior) |>
  #   my_write_csv(out)

    # kbl(escape = FALSE) |>
    # kable_classic()
}

write_sap_table <- function(sap_summary, out) {
  sap_summary |>
    janitor::clean_names() |>
    filter(variable %in% c("alpha", "beta")) |>
    mutate(variable = case_when(
      variable == "alpha" ~ "Intercept",
      variable == "beta" ~ "Slope",
    )) |>
    mutate(across(where(is.numeric), round, 2)) |>
    rename(Factor = variable) |>
    mutate(Posterior = paste0(q50, " [", q2_5, ", ", q97_5, "]")) |>
    dplyr::select(Factor, Posterior) |>
    my_write_csv(out)
}

generate_species_ab_table <- function(summary_sep) {
  map_dfr(summary_sep$fit_species,
    \(x)x$summary |> filter(variable %in% c("log_a", "b")))  |>
    mutate(species = rep(summary_sep$species, each = 2)) |>
    mutate(across(c(q50, q2.5, q97.5),
                  ~ ifelse(str_detect(variable, "log_a"), exp(.x), .x))) |>
    mutate(level = "species", target = species) |>
    mutate(variable = ifelse(variable == "log_a", "a", "b")) |>
    mutate(variable_meaning = case_when(
      variable == "a" ~ "coefficient a",
      variable == "b" ~ "coefficient b",
    )) |>
    dplyr::select(
      variable, level, target, variable_meaning, q50, q2.5, q97.5, rhat,
      effective_sample_size = ess_tail
    )
}

generate_segments_ab_table <- function(summary_sep) {
  summary_sep |>
    dplyr::select(species, fit_segments) |>
    unnest(cols = c(fit_segments)) |>
    filter(names(fit_segments) == "summary") |>
    unnest(cols = c(fit_segments)) |>
    filter(str_detect(variable, "log_a$|b$|A")) |>
    mutate(tmp = str_extract(variable, "(?<=,)[0-9]+(?=\\])")) |>
    mutate(across(c(q50, q2.5, q97.5),
                  ~ ifelse(str_detect(variable, "log_a|A\\[1"), exp(.x), .x))) |>
    mutate(level = ifelse(str_detect(variable, "A"),
      "segment",
      "species")) |>
    mutate(target = ifelse(str_detect(variable, "A"),
      paste("segment", tmp, "of", species),
      species)) |>
    mutate(variable = ifelse(variable == "log_a", "a", variable)) |>
    mutate(variable_meaning = case_when(
      variable == "a" ~ "coefficient a",
      variable == "b" ~ "coefficient b",
      str_detect(variable, "A\\[1") ~ "coefficients a",
      str_detect(variable, "A\\[2") ~ "coefficients b",
    )) |>
    dplyr::select(
      variable, level, target, variable_meaning, q50, q2.5, q97.5, rhat,
      effective_sample_size = ess_tail
    )
}

generate_species_segments_ab_csv <- function(df, out) {
  df |>
    dplyr::select(
      variable_name = variable,
      level, target,
      variable_meaning,
      max_pg,
      q50, q2.5, q97.5,
      rhat,
      effective_sample_size) |>
    mutate_if(is.numeric, \(x) round(x, digits = 2)) |>
    mutate(effective_sample_size = round(effective_sample_size, digits = 0)) |>
    mutate_if(is.numeric, \(x) format(x, nsmall = 2)) |>
    mutate(effective_sample_size = format(effective_sample_size, nsmall = 0)) |>
    my_write_csv(out)
}
