library(targets)
library(tarchetypes)
library(tidyverse)
library(stantargets)
# library(cmdstanr)
library(furrr)
library(clustermq)
library(bayesplot)
# library(doParallel)

source("R/data_clean.R")
source("R/stan.R")
source("R/figs.R")
source("R/tables.R")
source("R/scale.R")
source("R/traits.R")

plan(multicore)
options(clustermq.scheduler = "multicore")

set.seed(123)

tar_option_set(packages = c(
  "tidyverse",
  "patchwork",
  "cowplot",
  "bayesplot",
  "cmdstanr",
  "smatr",
  "ggsma",
  "ggpubr",
  "ggridges",
  "RColorBrewer",
  "scales",
  "loo",
  "here",
  "VIM",
  "lubridate",
  "foreach",
  "doParallel",
  "missForest",
  "data.table",
  "gtable",
  "ggrepel",
  "viridis",
  "ggpointdensity"
))

# tar_option_set(
#   garbage_collection = TRUE,
#   memory = "transient"
# )

pg <- c(seq(0.02, 0.08, by = 0.01), 0.025, 0.035)

# Register the parallel backend
n_cores <- parallel::detectCores(logical = FALSE)  # Detect the number of available CPU cores
# cl <- parallel::makeCluster(n_cores - 1)  # Create a cluster with one less core than available
cl <- 8
doParallel::registerDoParallel(cl)  # Register the parallel backend

# Target to list all Stan files
stan_files <- list.files(path = "stan", pattern = "\\.stan$", full.names = TRUE)

# Dynamic branching target to compile Stan models
compile_stan_models <- tar_map(
  values = tibble(stan_file = stan_files),
  # names = "model",
  tar_target(
    stan_model,
    {
      model_path <- stan_file
      # model_name <- sub("\\.stan$", "", basename(model_path))
      cmdstan_model(model_path)
    }
  )
)

# list(compile_stan_models)

 # raw data ----------------------------------
raw_data_list <- list(
  tar_target(
    five_spp_csv,
    "data-raw/pres_tens_five_spp.csv",
    format = "file"
  ),
  tar_target(
    ks_five_spp_csv,
    "data-raw/ks_pres_tens_five_spp.csv",
    format = "file"
  ),
  tar_target(
    ks_five_trees_raw_csv,
    "data-raw/ks_pres_tens_trees_raw.csv",
    format = "file"
  ),
  tar_target(
    cond_count_raw_csv,
    "data-raw/cond_count_raw.csv",
    format = "file"
  ),
  tar_target(
    ks_trees_csv,
    clean_ks_trees(ks_five_trees_raw_csv),
    format = "file"
  ),
  tar_target(
    ks_spp_err_csv,
    clean_ks_trees_err(ks_trees_csv),
    format = "file"
  ),
  tar_target(
    cond_count_csv,
    clean_cond_count(cond_count_raw_csv,
    "data/cond_count.csv"),
    format = "file"
  ),
  tar_target(
    calibration_raw_data_csv,
    "data-raw/calibration_raw_data.csv",
    format = "file"
  ),
  tar_target(
    rubber_raw_data_csv,
    "data-raw/rubber_raw_data.csv",
    format = "file"
  ),
  tar_target(
    girth_increment_csv,
    "data-raw/girth_increment.csv",
    format = "file"
  ),
  tar_target(
    initial_dbh_csv,
    "data-raw/initial_dbh.csv",
    format = "file"
  ),
  tar_target(
    sapwood_depth_csv,
    "data-raw/sapwood_depth.csv",
    format = "file"
  ),
  tar_target(
    pub_ab_csv,
    "data-raw/pub_ab.csv",
    format = "file"
  ),
  NULL
)

 # main analysis ----------------------------------
main_list <- list(
  tar_target(
    sma_scatter_plot,
    sma_scatter(five_spp_csv)
  ),
  # tar_target(
  #   ks_box_plot,
  #   ks_box(ks_trees_csv)
  # ),
  tar_target(
    sma_scatter_log_plot,
    sma_scatter(five_spp_csv, log = TRUE)
  ),
  tar_target(
    sma_scatter_ks_plot,
    sma_scatter(ks_five_spp_csv)
  ),
  tar_target(
    sma_scatter_ks_log_plot,
    sma_scatter(ks_five_spp_csv, log = TRUE)
  ),

  tar_target(
    dummy_data,
    generate_dummy_data(n = 3, mu_hat = -2, n_beta = 3, sigma_alpha = .8, sigma_beta = 2, sigma_gamma = 0.1, sigma = 0.3, seed = 123)
  ),

  tar_map(
    values = list(log = c("log", "nolog")),
    tar_target(
      anova_data,
      generate_anova_data(ks_spp_err_csv, log = ifelse(log == "log", TRUE, FALSE))
    ),
    tar_stan_mcmc(
       fit,
       c("stan/anova.stan", "stan/anova_noint.stan"),
       data = anova_data,
       refresh = 0,
       chains = 4,
       parallel_chains = 1,
       iter_warmup = 1000,
       iter_sampling = 1000,
       adapt_delta = 0.99,
       max_treedepth = 15,
       seed = 123,
       return_draws = TRUE,
       return_diagnostics = TRUE,
       return_summary = TRUE,
       summaries = list(
         mean = ~mean(.x),
         sd = ~sd(.x),
         mad = ~mad(.x),
         ~posterior::quantile2(.x, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)),
         posterior::default_convergence_measures()
      )
    )
  ),
  tar_stan_mcmc(
     fit_dummy,
     c("stan/anova.stan", "stan/anova_noint.stan"),
     data = dummy_data,
     refresh = 0,
     chains = 4,
     parallel_chains = 1,
     iter_warmup = 1000,
     iter_sampling = 1000,
     adapt_delta = 0.9,
     max_treedepth = 15,
     seed = 123,
     return_draws = TRUE,
     return_diagnostics = TRUE,
     return_summary = TRUE,
     summaries = list(
       mean = ~mean(.x),
       sd = ~sd(.x),
       mad = ~mad(.x),
       ~posterior::quantile2(.x, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)),
       posterior::default_convergence_measures()
    )
  ),

  tar_target(
    loo_anova,
    lapply(
      list(
           fit_mcmc_anova_log,
           fit_mcmc_anova_nolog,
           fit_mcmc_anova_noint_log,
           fit_mcmc_anova_noint_nolog
        ),
    \(x)x$loo(cores = parallel::detectCores())
    )
  ),

  tar_target(
    sma_ks_plot, {
      p <- sma_ks(five_spp_csv, ks_trees_csv, log = TRUE)
      my_ggsave(
        "figs/sma_ks",
        p,
        dpi = 600,
        width = 6.81,
        height = 4.4
      )
    },
    format = "file"
  ),

  tar_target(
    ks_box_plot, {
      p <- ks_box(ks_trees_csv)
      my_ggsave(
        "figs/ks_box",
        p,
        dpi = 600,
        width = 82,
        height = 100,
        units = "mm"
      )
    },
    format = "file"
  ),

  tar_target(
    coef_intervals_pres_tens_plot, {
      p1 <- coef_intervals_sd(fit_draws_anova_noint_log)
      p2 <- coef_intervals_mean(fit_draws_anova_noint_log)
      p3 <- coef_intervals_diff(fit_draws_anova_noint_log)
      p <- p1 + p2 + p3 + plot_spacer() +
        plot_layout(ncol = 2) +
        plot_annotation(tag_levels = "A")
      my_ggsave(
        "figs/coef_intervals_pres_tens",
        p,
        dpi = 600,
        width = 173,
        height = 173,
        units = "mm"
      )
    },
    format = "file"
  ),


  tar_target(
    anova_yml,
    write_anova_yml(
      "yml/anova.yml",
      fit_draws_anova_noint_log,
      ll = 0.25, hh = 0.75),
    format = "file"
  ),

  # fig2 ------------------------------------------
  tar_target(
   piecewise_logistic_stan_data,
   generate_piecewise_logistic_stan_data(cond_count_csv)
  ),
  tar_stan_mcmc(
    piecewise_logistic,
    "stan/piecewise_logistic.stan",
    data = piecewise_logistic_stan_data,
    refresh = 0,
    chains = 4,
    parallel_chains = 1,
    iter_warmup = 2000,
    iter_sampling = 2000,
    adapt_delta = 0.9,
    max_treedepth = 15,
    seed = 123,
    return_draws = TRUE,
    return_diagnostics = TRUE,
    return_summary = TRUE,
    summaries = list(
      mean = ~mean(.x),
      sd = ~sd(.x),
      mad = ~mad(.x),
      ~posterior::quantile2(.x, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)),
      posterior::default_convergence_measures()
    )
  ),

  tar_stan_mcmc(
    quad_logistic,
    "stan/hierarchical_logistic.stan",
    data = generate_logistic_stan_data(cond_count_csv, quad = TRUE),
    refresh = 0,
    chains = 4,
    parallel_chains = 1,
    iter_warmup = 2000,
    iter_sampling = 2000,
    adapt_delta = 0.9,
    max_treedepth = 15,
    seed = 123,
    return_draws = TRUE,
    return_diagnostics = TRUE,
    return_summary = TRUE,
    summaries = list(
      mean = ~mean(.x),
      sd = ~sd(.x),
      mad = ~mad(.x),
      ~posterior::quantile2(.x, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)),
      posterior::default_convergence_measures()
    )
  ),
  tar_stan_mcmc(
    simple_logistic,
    "stan/hierarchical_logistic.stan",
    data = generate_logistic_stan_data(cond_count_csv, quad = FALSE),
    refresh = 0,
    chains = 4,
    parallel_chains = 1,
    iter_warmup = 2000,
    iter_sampling = 2000,
    adapt_delta = 0.9,
    max_treedepth = 15,
    seed = 123,
    return_draws = TRUE,
    return_diagnostics = TRUE,
    return_summary = TRUE,
    summaries = list(
      mean = ~mean(.x),
      sd = ~sd(.x),
      mad = ~mad(.x),
      ~posterior::quantile2(.x, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)),
      posterior::default_convergence_measures()
    )
  ),

  tar_target(
    logistic_sp_plot, {
      p <- plot_logistic_sp(quad_logistic_draws_hierarchical_logistic, "data/cond_count.csv")
      my_ggsave(
        "figs/count_pressure_quadratic",
        p,
        dpi = 600,
        width = 4.3,
        height = 14,
        units = "cm"
      )
    },
    format = "file"
  ),

  tar_target(
    coef_intervals_logistic_plot, {
      p <- coef_intervals_logistic(quad_logistic_draws_hierarchical_logistic)
      my_ggsave(
        "figs/coef_intervals_logistic",
        p,
        dpi = 600,
        width = 16,
        height = 8.5,
        units = "cm"
      )
    },
    format = "file"
  ),

# dummy for sap ---------
  tar_target(
    dummy_sap_stan_data,
    generate_dummy_data_ab()
  ),

  # simple -------------------
  tar_target(
    fd_k_traits_csv,
    clean_sap_data(calibration_raw_data_csv, "data/fd_k_traits.csv"),
    format = "file"
  ),

  tar_target(
    sap_stan_data,
    generate_sap_stan_data(fd_k_traits_csv, upper_pressure = 0.03)
  ),


  tar_target(
    species_only_file,
    compile_model("stan/species_only.stan"),
    format = "file"
  ),
  tar_target(
    segments_inclusive_file,
    compile_model("stan/segments_inclusive.stan"),
    format = "file"
  ),
  tar_target(
    model1_like_file,
    compile_model("stan/model1_like.stan"),
    format = "file"
  ),
  tar_target(
    model2_like_file,
    compile_model("stan/model2_like.stan"),
    format = "file"
  ),


#   tar_target(
#     all_seg_table,
#     generate_summary_trait_table(
#       fit_abt_summary_granier_with_traits_sap_trait_clean_all,
#       fd_k_traits_csv
#     )
#   ),

#   tar_target(
#    varpart_notrait_table,
#    write_varpart_notrait_table(
#      ab_var_clean_008,
#      "data/varpart_notrait.csv"),
#    format = "file"
#   ),

#   tar_target(
#     d2015_imputed, {
#     rubber_impute(rubber_raw_data_csv, y2015 = TRUE)
#     }
#   ),
#   tar_target(
#     d2016_imputed, {
#     rubber_impute(rubber_raw_data_csv, y2015 = FALSE)
#     }
#   ),
  NULL
)

# without traits -------------------------------------------------------
granier_without_traits_mapped <- tar_map(
    list(p = c(0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.025, 0.035)),
    tar_target(sap_all_clean,
      generate_sap_stan_data(fd_k_traits_csv,
        remove_abnormal_values = TRUE,
        upper_pressure = p)),
    tar_target(sap_sp_clean,
      generate_sap_stan_data_sp(fd_k_traits_csv,
        remove_abnormal_values = TRUE,
        upper_pressure = p)),
    tar_stan_mcmc(
      fit,
      c("stan/species_xylem.stan",
        "stan/segments_xylem.stan"),
      data = sap_all_clean,
      refresh = 0,
      chains = 4,
      parallel_chains = getOption("mc.cores", 1),
      iter_warmup = 2000,
      iter_sampling = 2000,
      adapt_delta = 0.95,
      max_treedepth = 15,
      seed = 123,
      return_draws = TRUE,
      return_diagnostics = TRUE,
      return_summary = TRUE,
      summaries = list(
        mean = ~mean(.x),
        sd = ~sd(.x),
        mad = ~mad(.x),
        ~posterior::quantile2(.x, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)),
        posterior::default_convergence_measures()
      )
    ),
    tar_target(
      fit_each, {
      sap_sp_clean |>
          mutate(fit_species = map(stan_data, fit_model,
            species_only_file,
            adapt_delta = 0.95,
            iter_warmup = 2000,
            iter_sampling = 2000)) |>
          mutate(fit_segments = map(stan_data, fit_model,
            segments_inclusive_file,
            adapt_delta = 0.95,
            iter_warmup = 2000,
            iter_sampling = 2000))
      }),
    tar_target(
      species_ab_table,
      generate_species_ab_table(fit_each) |>
        mutate(max_pg = p)
    ),
    tar_target(
      segments_ab_table,
      generate_segments_ab_table(fit_each) |>
        mutate(max_pg = p)
    ),
    tar_target(
      species_only_post_ab, {
        generate_post_ab_each(fit_each) |> sample_n(1000)
      }
    ),
    tar_target(
      segments_inclusive_post_ab, {
        generate_post_ab_each(fit_each) |> sample_n(1000)
      }
    )
  )

like_check_list <- list(
  tar_stan_mcmc(
    fit_like,
    c("stan/model3_like.stan",
      "stan/model4_like.stan"),
    data = sap_all_clean_0.08,
    refresh = 0,
    chains = 4,
    parallel_chains = getOption("mc.cores", 4),
    iter_warmup = 2000,
    iter_sampling = 2000,
    adapt_delta = 0.95,
    max_treedepth = 15,
    seed = 123,
    return_draws = TRUE,
    return_diagnostics = FALSE,
    return_summary = FALSE
  ),
  tar_target(sap_sp_clean_like,
    generate_sap_stan_data_sp(fd_k_traits_csv,
      remove_abnormal_values = TRUE,
      upper_pressure = 0.8)),
  tar_target(
    fit_each_like, {
    sap_sp_clean_0.08 |>
        mutate(fit_species = map(stan_data, fit_model,
          model1_like_file,
          adapt_delta = 0.95,
          iter_warmup = 2000,
          iter_sampling = 2000)) |>
        mutate(fit_segments = map(stan_data, fit_model,
          model2_like_file,
          adapt_delta = 0.95,
          iter_warmup = 2000,
          iter_sampling = 2000))
    })
)

segments_xylem_df_mapped <- tar_map(
  list(stan_summary =
    rlang::syms(
    str_c("fit_summary_segments_xylem_",
      c(seq(0.02, 0.08, by = 0.01), 0.025, 0.035)))),
  tar_target(
    fit_summary,
    stan_summary
  ),
  tar_target(
    table,
    generates_segments_xylem_table(stan_summary)
  ),
  tar_target(
    table_re,
    generate_summary_non_trait_table(table)
  )
)

segments_xylem_draws_mapped <- tar_map(
  list(stan_draws =
    rlang::syms(
    str_c("fit_draws_segments_xylem_",
      c(seq(0.02, 0.08, by = 0.01), 0.025, 0.035)))),
  tar_target(
    fit_draws,
    stan_draws
  ),
  tar_target(
    draws, {
      stan_draws  |>
        janitor::clean_names() |>
        dplyr::select(starts_with("alpha"))
    }
  )
)

species_xylem_df_mapped <- tar_map(
  list(stan_summary =
    rlang::syms(
    str_c("fit_summary_species_xylem_",
      c(seq(0.02, 0.08, by = 0.01), 0.025, 0.035)))),
  tar_target(
    fit_summary,
    stan_summary
  ),
  tar_target(
    table,
    generates_segments_xylem_table(stan_summary)
  ),
  tar_target(
    table_re,
    generate_summary_non_trait_table(table)
  )
)

species_xylem_post_ab_mapped <- tar_map(
  list(stan_draws =
    rlang::syms(
    str_c("fit_draws_species_xylem_",
      c(seq(0.02, 0.08, by = 0.01), 0.025, 0.035)))),
    tar_target(
      species_xylem_post_ab,
      generate_post_ab(stan_draws) |> sample_n(1000)
    )
)

segments_xylem_post_ab_mapped <- tar_map(
  list(stan_draws =
    rlang::syms(
    str_c("fit_draws_segments_xylem_",
      c(seq(0.02, 0.08, by = 0.01), 0.025, 0.035)))),
    tar_target(
      segments_xylem_post_ab,
      generate_post_ab(stan_draws) |> sample_n(1000)
    )
)

tar_combined_segments_xylem_df <- tar_combine(
  segments_xylem_df_combined,
  segments_xylem_df_mapped[["table_re"]],
  command = dplyr::bind_rows(!!!.x, .id = "id")
)
tar_combined_segments_xylem_draws <- tar_combine(
  segments_xylem_draws_combined,
  segments_xylem_draws_mapped[["draws"]],
  command = dplyr::bind_rows(!!!.x, .id = "id")
)

tar_combined_segments_xylem_summary <- tar_combine(
  segments_xylem_summary_combined,
  segments_xylem_df_mapped[["fit_summary"]],
  command = dplyr::bind_rows(!!!.x, .id = "id")
)

tar_combined_species_xylem_df <- tar_combine(
  species_xylem_df_combined,
  species_xylem_df_mapped[["table_re"]],
  command = dplyr::bind_rows(!!!.x, .id = "id")
)
tar_combined_species_xylem_summary <- tar_combine(
  species_xylem_summary_combined,
  species_xylem_df_mapped[["fit_summary"]],
  command = dplyr::bind_rows(!!!.x, .id = "id")
)

# with single traits -----------------------------------------
  granier_with_traits_mapped <- tar_map(
    values = list(trait_name = rlang::syms(c(
      "log_dh", "log_vf", "wood_density",
      "log_ks", "log_vaf", "log_swc"))),
    tar_target(
      stan_data_noxylem,
      generate_sap_each_trait_no_xylem_stan_data(
        fd_k_traits_csv,
        trait_name = trait_name,
        remove_abnormal_values = TRUE)
    ),
    tar_target(
      stan_data_xylem,
      generate_sap_each_trait_xylem_stan_data(
        fd_k_traits_csv,
        trait_name = trait_name,
        remove_abnormal_values = TRUE)
    ),
    tar_stan_mcmc(
      fit,
      "stan/segments_noxylem_traits.stan",
      data = stan_data_noxylem,
      refresh = 0,
      chains = 4,
      parallel_chains = 4,
      iter_warmup = 2000,
      iter_sampling = 2000,
      adapt_delta = 0.95,
      max_treedepth = 15,
      seed = 123,
      return_draws = TRUE,
      return_diagnostics = TRUE,
      return_summary = TRUE,
      summaries = list(
        mean = ~mean(.x),
        sd = ~sd(.x),
        mad = ~mad(.x),
        ~posterior::quantile2(.x, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)),
        posterior::default_convergence_measures()
        )
     ),
   tar_stan_mcmc(
      fit2,
      "stan/segments_xylem_traits.stan",
      data = stan_data_xylem,
      refresh = 0,
      chains = 4,
      parallel_chains = 4,
      iter_warmup = 4000,
      iter_sampling = 2000,
      adapt_delta = 0.99,
      max_treedepth = 15,
      seed = 123,
      return_draws = TRUE,
      return_diagnostics = TRUE,
      return_summary = TRUE,
      summaries = list(
        mean = ~mean(.x),
        sd = ~sd(.x),
        mad = ~mad(.x),
        ~posterior::quantile2(.x, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)),
        posterior::default_convergence_measures()
        )
    ),
   tar_stan_mcmc(
      fit3,
      "stan/segments_noxylem_traits_sp.stan",
      data = stan_data_noxylem,
      refresh = 0,
      chains = 4,
      parallel_chains = 1,
      iter_warmup = 2000,
      iter_sampling = 2000,
      adapt_delta = 0.95,
      max_treedepth = 15,
      seed = 123,
      return_draws = TRUE,
      return_diagnostics = TRUE,
      return_summary = TRUE,
      summaries = list(
        mean = ~mean(.x),
        sd = ~sd(.x),
        mad = ~mad(.x),
        ~posterior::quantile2(.x, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)),
        posterior::default_convergence_measures()
        )
    ),
   tar_stan_mcmc(
      fit4,
      "stan/segments_noxylem_traits_sp2.stan",
      data = stan_data_noxylem,
      refresh = 0,
      chains = 4,
      parallel_chains = 1,
      iter_warmup = 2000,
      iter_sampling = 2000,
      adapt_delta = 0.95,
      max_treedepth = 15,
      seed = 123,
      return_draws = TRUE,
      return_diagnostics = TRUE,
      return_summary = TRUE,
      summaries = list(
        mean = ~mean(.x),
        sd = ~sd(.x),
        mad = ~mad(.x),
        ~posterior::quantile2(.x, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)),
        posterior::default_convergence_measures()
        )
    )
  )

tar_combined_species_ab_table <- tar_combine(
  species_ab_table_combined,
  granier_without_traits_mapped[["species_ab_table"]],
  command = dplyr::bind_rows(!!!.x)
)
tar_combined_segments_ab_table <- tar_combine(
  segments_ab_table_combined,
  granier_without_traits_mapped[["segments_ab_table"]],
  command = dplyr::bind_rows(!!!.x)
)

segments_noxylem_traits_post_ab_mapped <- tar_map(
  list(stan_summary =
    rlang::syms(
    str_c("fit_summary_segments_noxylem_traits_",
      c("log_dh", "log_vf", "wood_density", "log_ks", "log_vaf"))),
    key = c("log_dh", "log_vf", "wood_density", "log_ks", "log_vaf")),
    tar_target(
      post,
      generate_summary_trait_table(
        stan_summary, fd_k_traits_csv, xylem = FALSE) |>
          mutate(trait = key) |>
          dplyr::select(trait, everything())
    )
)
segments_xylem_traits_post_ab_mapped <- tar_map(
  list(stan_summary =
    rlang::syms(
    str_c("fit2_summary_segments_xylem_traits_",
      c("log_dh", "log_vf", "wood_density", "log_ks", "log_vaf"))),
    key = c("log_dh", "log_vf", "wood_density", "log_ks", "log_vaf")),
    tar_target(
      post,
      generate_summary_trait_table(
        stan_summary, fd_k_traits_csv, xylem = TRUE) |>
          mutate(trait = key) |>
          dplyr::select(trait, everything())
    )
)

tar_combined_segments_noxylem_traits_table <- tar_combine(
  segments_noxylem_traits_table_combined,
  segments_noxylem_traits_post_ab_mapped[["post"]],
  command = dplyr::bind_rows(!!!.x)
)
tar_combined_segments_xylem_traits_table <- tar_combine(
  segments_xylem_traits_table_combined,
  segments_xylem_traits_post_ab_mapped[["post"]],
  command = dplyr::bind_rows(!!!.x)
)


# granier analysis -----------------------------------
granier_list <- list(
  granier_with_traits_mapped,
  granier_without_traits_mapped,
  tar_combined_species_ab_table,
  tar_combined_segments_ab_table,
  tar_combined_segments_xylem_summary,
  tar_combined_species_xylem_summary,
  species_xylem_post_ab_mapped,
  segments_xylem_post_ab_mapped,
  segments_xylem_traits_post_ab_mapped,
  segments_noxylem_traits_post_ab_mapped,
  tar_combined_segments_noxylem_traits_table,
  tar_combined_segments_xylem_traits_table,
  segments_xylem_draws_mapped,
  tar_combined_segments_xylem_draws,
  tar_target(
    fit_draws_segments_xylem_combined, {
      segments_xylem_draws_combined |>
      mutate(max_pg = str_extract(id, "\\d+\\.\\d+"))
    }
  ),
  tar_target(
    segments_inclusive_ab_csv,
    generate_species_segments_ab_csv(segments_ab_table_combined, "data/segments_inclusive_ab.csv"),
    format = "file"
  ),
  tar_target(
    species_only_ab_csv,
    generate_species_segments_ab_csv(species_ab_table_combined, "data/species_only_ab.csv"),
    format = "file"
  ),
  segments_xylem_df_mapped,
  species_xylem_df_mapped,
  tar_combined_segments_xylem_df,
  tar_combined_species_xylem_df,
  tar_target(
    segments_xylem_csv, {
       segments_xylem_df_combined |>
       mutate(max_pg = str_extract(id, "\\d+\\.\\d+")) |>
       dplyr::select(-id) |>
       my_write_csv("data/segments_xylem_post.csv")
    },
    format = "file"
  ),
  tar_target(
    species_xylem_csv, {
       species_xylem_df_combined |>
       mutate(max_pg = str_extract(id, "\\d+\\.\\d+")) |>
       dplyr::select(-id) |>
       my_write_csv("data/species_xylem_post.csv")
    },
    format = "file"
  ),
  tar_target(
    traits_xylem_table,
    generate_summary_trait_table(
      fit2_summary_segments_xylem_traits_log_ks, fd_k_traits_csv, xylem = TRUE)
  ),
  tar_target(
    traits_noxylem_table,
    generate_summary_trait_table(
      fit_summary_segments_noxylem_traits_log_ks, fd_k_traits_csv, xylem = FALSE)
  ),
  tar_target(
    ab_csv,
    write_ab_csv(
      fit_summary_segments_xylem_0.08,
      fit_summary_species_xylem_0.08,
      fit_each_0.08,
      xylem_lab,
      "data/all_ab.csv"),
    format = "file"
  ),
  tar_target(
    xylem_lab,
    generate_xylem_lab(fd_k_traits_csv)
  ),
  tar_target(
    traits_xylem_table_csv,
    my_write_csv(segments_xylem_traits_table_combined, "data/traits_xylem_post.csv"),
    format = "file"
  ),
  tar_target(
    traits_noxylem_table_csv,
    my_write_csv(segments_noxylem_traits_table_combined, "data/traits_noxylem_post.csv"),
    format = "file"
  ),

  tar_target(
    pool_multi_plot, {
      p <- line_pool_multi(fd_k_traits_csv,
       xylem_lab,
       fit_summary_segments_xylem_0.08,
       fit_summary_species_xylem_0.08)
      my_ggsave(
        "figs/pool_multi",
        p,
        dpi = 600,
        width = 8,
        height = 12
      )
    },
    format = "file"
  ),
  tar_target(
    ab_points_four_models_plot_all, {
      p <- ab_comp_four_models_points(
        fit_each_0.08,
        fit_summary_species_xylem_0.08,
        fit_summary_segments_xylem_0.08,
        xylem_lab, rm_dip = FALSE)
      my_ggsave(
        "figs/ab_points_four_models_all",
        p,
        dpi = 600,
        width = 110,
        height = 165,
        units = "mm"
      )
    },
    format = "file"
  ),
  tar_target(
    ab_points_four_models_plot, {
      p <- ab_comp_four_models_points(
        fit_each_0.08,
        fit_summary_species_xylem_0.08,
        fit_summary_segments_xylem_0.08,
        xylem_lab, rm_dip = TRUE)
      my_ggsave(
        "figs/ab_points_four_models",
        p,
        dpi = 600,
        width = 110,
        height = 165,
        units = "mm"
      )
    },
    format = "file"
  ),
  tar_target(
    k_range,
    generate_k_range(fd_k_traits_csv)
  ),
  tar_target(
    fit_summary_segments_xylem_combined, {
      segments_xylem_summary_combined |>
      mutate(max_pg = str_extract(id, "\\d+\\.\\d+")) |>
      janitor::clean_names()
    }
  ),
  tar_target(
    pg_multi_plot, {
      p <- line_pg_multi(
        fd_k_traits_csv,
        xylem_lab,
        k_range,
        fit_summary_segments_xylem_combined)
      my_ggsave(
        "figs/pg_multi",
        p,
        dpi = 600,
        width = 8,
        height = 12
      )
    },
    format = "file"
  ),
  tar_target(
    pg_ribbon_a, {
      p <- ab_pg_ribbon(xylem_lab, k_range, fit_summary_segments_xylem_combined)
      my_ggsave(
        "figs/pg_ribbon_a",
        p,
        dpi = 600,
        width = 8,
        height = 12
      )
    },
    format = "file"
  ),
  tar_target(
    pg_ribbon_b, {
      p <- ab_pg_ribbon(xylem_lab, k_range, fit_summary_segments_xylem_combined,
        coef_a = FALSE)
      my_ggsave(
        "figs/pg_ribbon_b",
        p,
        dpi = 600,
        width = 8,
        height = 12
      )
    },
    format = "file"
  ),
  tar_target(
    coef_density_plot, {
      p <- coef_density(xylem_lab,
        fit_draws_segments_xylem_0.08)
      my_ggsave(
        "figs/coef_density",
        p,
        dpi = 600,
        width = 6.8,
        height = 5
      )
    },
    format = "file"
  ),
  tar_target(
    # tar_read(trait_fig_data_combined)
    trait_pred_data_noxylem_combined,
    generate_combined_trait_fig_data(
      summary = list(
        fit_summary_segments_noxylem_traits_log_vaf,
        fit_summary_segments_noxylem_traits_log_ks,
        fit_summary_segments_noxylem_traits_wood_density,
        fit_summary_segments_noxylem_traits_log_swc,
        fit_summary_segments_noxylem_traits_log_dh,
        fit_summary_segments_noxylem_traits_log_vf),
      draws = list(
        fit_draws_segments_noxylem_traits_log_vaf,
        fit_draws_segments_noxylem_traits_log_ks,
        fit_draws_segments_noxylem_traits_wood_density,
        fit_draws_segments_noxylem_traits_log_swc,
        fit_draws_segments_noxylem_traits_log_dh,
        fit_draws_segments_noxylem_traits_log_vf),
      fd_k_traits_csv,
      xylem_lab,
      no_xylem = TRUE,
      single_trait = TRUE,
      sp_level = FALSE)
  ),
  tar_target(
    trait_pred_data_xylem_combined,
    generate_combined_trait_fig_data(
      summary = list(
        fit2_summary_segments_xylem_traits_log_vaf,
        fit2_summary_segments_xylem_traits_log_ks,
        fit2_summary_segments_xylem_traits_wood_density,
        fit2_summary_segments_xylem_traits_log_swc,
        fit2_summary_segments_xylem_traits_log_dh,
        fit2_summary_segments_xylem_traits_log_vf),
      draws = list(
        fit2_draws_segments_xylem_traits_log_vaf,
        fit2_draws_segments_xylem_traits_log_ks,
        fit2_draws_segments_xylem_traits_wood_density,
        fit2_draws_segments_xylem_traits_log_swc,
        fit2_draws_segments_xylem_traits_log_dh,
        fit2_draws_segments_xylem_traits_log_vf),
      fd_k_traits_csv,
      xylem_lab,
      no_xylem = FALSE,
      single_trait = TRUE)
  ),
  tar_target(
    trait_pred_data_noxylem_sp_combined,
    generate_combined_trait_fig_data(
      summary = list(
        fit3_summary_segments_noxylem_traits_sp_log_vaf,
        fit3_summary_segments_noxylem_traits_sp_log_ks,
        fit3_summary_segments_noxylem_traits_sp_wood_density,
        fit3_summary_segments_noxylem_traits_sp_log_swc,
        fit3_summary_segments_noxylem_traits_sp_log_dh,
        fit3_summary_segments_noxylem_traits_sp_log_vf),
      draws = list(
        fit3_draws_segments_noxylem_traits_sp_log_vaf,
        fit3_draws_segments_noxylem_traits_sp_log_ks,
        fit3_draws_segments_noxylem_traits_sp_wood_density,
        fit3_draws_segments_noxylem_traits_sp_log_swc,
        fit3_draws_segments_noxylem_traits_sp_log_dh,
        fit3_draws_segments_noxylem_traits_sp_log_vf),
      fd_k_traits_csv,
      xylem_lab,
      no_xylem = TRUE,
      single_trait = TRUE,
      sp_level = TRUE)
  ),
  tar_target(
    traits_points_main_plot, {
      p <- traits_points_main(trait_pred_data_noxylem_combined)
      my_ggsave(
        "figs/traits_points_main",
        p,
        dpi = 600,
        width = 6.8,
        height = 2.7
      )
    },
    format = "file"
  ),
  # tar_target(
  #   traits_points_sp_main_plot, {
  #     p <- traits_points_main(trait_pred_data_noxylem_sp_combined)
  #     my_ggsave(
  #       "figs/traits_points_sp_main",
  #       p,
  #       dpi = 600,
  #       width = 6.8,
  #       height = 2.7
  #     )
  #   },
  #   format = "file"
  # ),
  tar_target(
    traits_points_main_re_plot, {
      p <- traits_points_main_re(trait_pred_data_noxylem_combined)
      my_ggsave(
        "figs/traits_points_main_re",
        p,
        dpi = 600,
        width = 110,
        height = 110,
        units = "mm"
      )
    },
    format = "file"
  ),
  tar_target(
    traits_sp_points_main_plot, {
      p <- traits_sp_points_main(
        trait_pred_data_noxylem_combined,
        trait_pred_data_noxylem_sp_combined,
        vaf_r2,
        ks_r2,
        vaf_sp_r2,
        ks_sp_r2
        )
      my_ggsave(
        "figs/traits_sp_points_main",
        p,
        dpi = 600,
        width = 173,
        height = 110,
        units = "mm"
      )
    },
    format = "file"
  ),
  tar_target(
    traits_seg_points_si_plot, {
      p <- traits_points_si(
        trait_pred_data_noxylem_combined,
        r2_list = list(
          wd_r2,
          swc_r2,
          dh_r2,
          vf_r2
        ),
        title = "Segments"
        )
      my_ggsave(
        "figs/traits_seg_points_si",
        p,
        dpi = 600,
        width = 173,
        height = 110,
        units = "mm"
      )
    },
    format = "file"
  ),
  tar_target(
    traits_sp_points_si_plot, {
      p <- traits_points_si(
        trait_pred_data_noxylem_sp_combined,
        r2_list = list(
          wd_sp_r2,
          swc_sp_r2,
          dh_sp_r2,
          vf_sp_r2
        ),
        title = "Species",
        sp = TRUE
        )
      my_ggsave(
        "figs/traits_sp_points_si",
        p,
        dpi = 600,
        width = 173,
        height = 110,
        units = "mm"
      )
    },
    format = "file"
  ),
  tar_target(
    ks_r2,
    process_draws_and_calculate_trait_r2(fit_draws_segments_noxylem_traits_log_ks, stan_data_noxylem_log_ks$xj)
  ),
  tar_target(
    ks_sp_r2,
    process_draws_and_calculate_trait_r2(fit3_draws_segments_noxylem_traits_sp_log_ks, stan_data_noxylem_log_ks$xk, sp_level = TRUE)
  ),
  tar_target(
    vaf_r2,
    process_draws_and_calculate_trait_r2(fit_draws_segments_noxylem_traits_log_vaf, stan_data_noxylem_log_vaf$xj)
  ),
  tar_target(
    vaf_sp_r2,
    process_draws_and_calculate_trait_r2(fit3_draws_segments_noxylem_traits_sp_log_vaf, stan_data_noxylem_log_vaf$xk, sp_level = TRUE)
  ),
  tar_target(
    dh_r2,
    process_draws_and_calculate_trait_r2(fit_draws_segments_noxylem_traits_log_dh, stan_data_noxylem_log_dh$xj)
  ),
  tar_target(
    dh_sp_r2,
    process_draws_and_calculate_trait_r2(fit3_draws_segments_noxylem_traits_sp_log_dh, stan_data_noxylem_log_dh$xk, sp_level = TRUE)
  ),
  tar_target(
    vf_r2,
    process_draws_and_calculate_trait_r2(fit_draws_segments_noxylem_traits_log_vf, stan_data_noxylem_log_vf$xj)
  ),
  tar_target(
    vf_sp_r2,
    process_draws_and_calculate_trait_r2(fit3_draws_segments_noxylem_traits_sp_log_vf, stan_data_noxylem_log_vf$xk, sp_level = TRUE)
  ),
  tar_target(
    swc_r2,
    process_draws_and_calculate_trait_r2(fit_draws_segments_noxylem_traits_log_swc, stan_data_noxylem_log_swc$xj)
  ),
  tar_target(
    swc_sp_r2,
    process_draws_and_calculate_trait_r2(fit3_draws_segments_noxylem_traits_sp_log_swc, stan_data_noxylem_log_swc$xk, sp_level = TRUE)
  ),
  tar_target(
    wd_r2,
    process_draws_and_calculate_trait_r2(fit_draws_segments_noxylem_traits_wood_density, stan_data_noxylem_wood_density$xj)
  ),
  tar_target(
    wd_sp_r2,
    process_draws_and_calculate_trait_r2(fit3_draws_segments_noxylem_traits_sp_wood_density, stan_data_noxylem_wood_density$xk, sp_level = TRUE)
  ),
  tar_target(
    traits_points_si_plot, {
      p <- traits_points_si(trait_pred_data_xylem_combined)
      my_ggsave(
        "figs/traits_points_si",
        p,
        dpi = 600,
        width = 6.8,
        height = 2.9
      )
    },
    format = "file"
  ),
  tar_target(
    ab_pg_summary_bars_plot, {
      p <- ab_pg_summary_bars(
        s = fit_summary_segments_xylem_combined,
        d = fit_draws_segments_xylem_combined,
        xylem_lab)
      my_ggsave(
        "figs/ab_pg_summary_bars",
        p,
        dpi = 200,
        width = 173,
        height = 120,
        units = "mm"
      )
    },
    format = "file"
  ),
  tar_target(
    ab_points_model4_sma_list,
    ab_points_model4_sma(
        fit_summary_segments_xylem_0.08,
        fd_k_traits_csv,
        xylem_lab,
        pub_ab_csv
     )
  ),
  tar_target(
    ab_points_model4_plot, {
      p <- ab_points_model4(
        fit_summary_segments_xylem_0.08,
        fd_k_traits_csv,
        xylem_lab,
        pub_ab_csv
        )
      my_ggsave(
        "figs/ab_points_model4",
        p,
        dpi = 600,
        width = 173,
        height = 90,
        units = "mm"
      )
    },
    format = "file"
  ),
  NULL
)

# imputation ------------------------------------------------------------
values <- expand_grid(year = c(2015, 2016), month = 1:12) |>
  filter(!(year == 2016 & (month == 12 | month == 10 | month == 9)))

values2 <- expand_grid(year = 2016, month = c(9, 10, 12))
  # mutate(df_name = rlang::syms(paste0("imputed_df_", year - 1, "_", month)))
impute_mapped <- tar_map(
  values = values,
  tar_target(
    cleaned_df_missforest,
      clean_for_missForest(
        csv = rubber_raw_data_csv,
        year = year,
        month = month
      )
  ),
  tar_target(
    imputed_data,
      missForest::missForest(cleaned_df_missforest,
        parallelize = "forests")
  ),
  tar_target(
    imputed_df, {
      imputed_data$ximp |>
        as_tibble()
    }
    ),
  tar_target(
    nramse,
    imputed_data$OOBerror[1]
    ),
  tar_target(
    imputed_df_btrans,
    backtransform_date(rubber_raw_data_csv, year, month, imputed_df)
    )
  )

tar_combined_imputed_data <- tar_combine(
  combined_imputed_mapped,
  impute_mapped[["imputed_df_btrans"]],
  command = dplyr::bind_rows(!!!.x)
)

tar_combined_nrmse <- tar_combine(
  combined_nrmse,
  impute_mapped[["nramse"]],
  command = dplyr::bind_rows(!!!.x, .id = "id")
)

impute_rest_mapped <- tar_map(
  values = values2,
  tar_target(
    imputed_data,
    {
      if (month == 9) {
        df <- imputed_df_2015_9
      } else if (month == 10) {
        df <- imputed_df_2015_10
      } else if (month == 12) {
        df <- imputed_df_2015_12
      }
      tmp <- clean_for_missForest(
        csv = rubber_raw_data_csv,
        year = year,
        month = month)
      tmp2 <- bind_rows(df, tmp)
      as.data.frame(tmp2[,-1]) |>
        missForest(parallelize = "forests")
    }
  ),
  tar_target(
    imputed_df,
    clean_imputed_df(imputed_data)
   ),
  tar_target(
    nramse,
    imputed_data$OOBerror[1]
    ),
  tar_target(
    imputed_df_btrans,
    backtransform_date(rubber_raw_data_csv, year, month, imputed_df)
  )
)

tar_combined_nrmse_rest <- tar_combine(
  combined_nrmse_rest,
  impute_rest_mapped[["nramse"]],
  command = dplyr::bind_rows(!!!.x, .id = "id")
)
tar_combined_imputed_rest_data <- tar_combine(
  combined_imputed_rest_mapped,
  impute_rest_mapped[["imputed_df_btrans"]],
  command = dplyr::bind_rows(!!!.x)
)

tar_impute <- list(
  impute_mapped,
  tar_combined_imputed_data,
  tar_combined_nrmse,
  impute_rest_mapped,
  tar_combined_imputed_rest_data,
  tar_combined_nrmse_rest,
  tar_target(
    nrmse_df, {
      tmp1 <- combined_nrmse |>
        mutate(id = str_extract(id, "(?<=nramse_).*$"))
      tmp2 <- combined_nrmse_rest |>
        mutate(id = str_extract(id, "(?<=imputed_df_).*$"))
      bind_rows(tmp1, tmp2) |>
        arrange(id)
    }
  ),
  tar_target(
    imputed_full_df, {
      bind_rows(combined_imputed_mapped, combined_imputed_rest_mapped) |>
      mutate(date = ymd(paste(year, "01", "01", sep= "-")) + days(yday - 1)) |>
      arrange(date) |>
      arrange(dir) |>
      arrange(dep) |>
      arrange(tree) |>
      mutate(h = time %/% 60) |>
      mutate(m = time %% 60) |>
      mutate(time = sprintf("%02d:%02d:%02d", h, m, 0)) |>
      dplyr::select(year, date, time, vpd, par, k, tree, dir, dep)
    }
  ),
  tar_target(
    nonimputed_full_df,
    make_long_nonimputed_df(rubber_raw_data_csv)
  ),
  NULL
  )

# scaling ---------------------------------------------------------

tar_dir_dep <- list(
  # put NA values for unmesaured direction and depth combinations
  tar_target(
    dir_dep_imp_df,
    generate_dir_dep_imp_data(
      imputed_full_df)
  ),
  tar_map(
    values = list(data_type = c("dir_only", "dep_only", "dir_dep")),
    tar_target(
      post,
      generate_post_dir_dep(fit2_draws_dir_dep, data_type)
    ),
    tar_target(
      post_1000,
      generate_post_dir_dep(fit2_draws_dir_dep, data_type) |> sample_n(1000)
    )
  ),
  tar_target(
    post_dir_dep_mid,
    apply(post_dir_dep, 2, median)
  ),
  tar_target(
    post_slen, {
      fit_draws_sap_dbh |>
        dplyr::select(alpha, beta, sigma)
    }
  ),
  tar_target(
    post_slen_mid,
    apply(post_slen, 2, median)
  ),
  tar_target(
    post_slen_1000, {
      post_slen |> sample_n(1000)
    },
  ),
  tar_target(
    sarea_df,
    generate_sarea_df(dbh_imp_df,
      post_slen_mid |> as.data.frame() |> t() |> as_tibble())
  ),
  NULL
)

tmp <- c("species_xylem_post_ab_fit_draws_species_xylem",
         "segments_xylem_post_ab_fit_draws_segments_xylem",
         "segments_inclusive_post_ab",
         "species_only_post_ab")
# pg <- 0.08
post_ab_names <- expand_grid(tmp, pg) |>
  mutate(tmp3 = paste(tmp, pg, sep = "_")) |> pull(tmp3)

post_ab_names <- post_ab_names[
  !grepl("^(species_only_post_ab_|segments_inclusive_post_ab_)", post_ab_names) |
  grepl("0.08$", post_ab_names)]

uncertainty_ab_mapped <- tar_map(
  values = tibble(
    post_ab_fit_draws = rlang::syms(post_ab_names)),
  tar_target(
    ab_uncertainty_df,
    generate_ab_uncertainty(
      post_ab_fit_draws,
      post_dir_dep_mid,
      sarea_df,
      dir_dep_imp_df,
      n_draws = 1000
  )),
  tar_target(
    ab_summarized_df,
    summarize_each_uncertainty(ab_uncertainty_df)
  ),
  tar_target(
    ab_scaled_df,
    summarize_stats_uncertainty(ab_summarized_df)
  ),
  NULL
 )

tar_combined_ab_uncertainty <- tar_combine(
  ab_uncertainty_combined_df,
  uncertainty_ab_mapped[["ab_scaled_df"]],
  command = dplyr::bind_rows(!!!.x, .id = "id")
)

uncertainty_dir_dep_mapped <- tar_map(
  values = tibble(
    post_dir_dep_fit_draws =
      rlang::syms(paste0("post_1000_", c("dir_only", "dep_only", "dir_dep")))),
  tar_target(
    dir_dep_uncertainty_df,
    generate_dir_dep_uncertainty(
      post_ab_mid,
      post_dir_dep_fit_draws,
      sarea_df,
      dir_dep_imp_df,
      n_draws = 1000
  )),
  tar_target(
    dir_dep_summarized_df,
    summarize_each_uncertainty(dir_dep_uncertainty_df)
  ),
  tar_target(
    dir_dep_scaled_df,
    summarize_stats_uncertainty(dir_dep_summarized_df)
  ),
  NULL
 )

tar_combined_dir_dep_uncertainty <- tar_combine(
  dir_dep_uncertainty_combined_df,
  uncertainty_dir_dep_mapped[["dir_dep_scaled_df"]],
  command = dplyr::bind_rows(!!!.x, .id = "id")
)

uncertainty_list <- list(
  uncertainty_ab_mapped,
  tar_combined_ab_uncertainty,
  tar_target(
    ab_granier_uncertainty_df,
    generate_ab_uncertainty(
      post_ab_fit_draws = NULL,
      post_dir_dep_mid,
      sarea_df,
      dir_dep_imp_df)
  ),
  tar_target(
    ab_granier_summarized_df,
    summarize_each_uncertainty(ab_granier_uncertainty_df)
  ),
  tar_target(
    ab_granier_uncertainty_combined_df,
    summarize_stats_uncertainty(ab_granier_summarized_df)
  ),
  tar_target(
    post_ab_mid,
    map_dbl(segments_xylem_post_ab_fit_draws_segments_xylem_0.08, median)
  ),
  uncertainty_dir_dep_mapped,
  tar_combined_dir_dep_uncertainty,
  tar_target(
    sarea_uncertainty_df,
    generate_sarea_uncertainty(
      post_ab_mid,
      post_dir_dep_mid,
      dbh_imp_df,
      post_slen_1000,
      dir_dep_imp_df,
      n_draws = 1000)
  ),
  tar_target(
    sarea_summarized_df,
    summarize_each_uncertainty(sarea_uncertainty_df)
  ),
  tar_target(
    sarea_uncertainty_combined_df,
    summarize_stats_uncertainty(sarea_summarized_df)
  ),
  tar_target(
    total_uncertainty_df,
    generate_total_uncertainty(
      segments_xylem_post_ab_fit_draws_segments_xylem_0.08,
      post_1000_dir_dep,
      dbh_imp_df,
      post_slen_1000,
      dir_dep_imp_df,
      n_draws = 1000)
  ),
  tar_target(
    total_summarized_df,
    summarize_each_uncertainty(total_uncertainty_df)
  ),
  tar_target(
    total_uncertainty_combined_df,
    summarize_stats_uncertainty(total_summarized_df)
  ),
  NULL
)

uncertainty_figs_list <- list(
  tar_target(
    tr_bar_ab_df,
    generate_tr_bar_ab_df(ab_uncertainty_combined_df, ab_granier_uncertainty_combined_df)
  ),
  tar_target(
    tr_bar_ab_csv,
    tr_bar_ab_df |> dplyr::select(-id) |>
      my_write_csv("data/tr_bar_ab.csv"),
    format = "file"
  ),
  # tar_target(
  #   tr_bar_ab_plot, {
  #     p <- tr_bar_ab(tr_bar_ab_df, pg, tr_m, fill = model, group = model)
  #     my_ggsave(
  #       "figs/tr_bar_ab",
  #       p,
  #       dpi = 600,
  #       width = 6.81,
  #       height = 4.4
  #     )
  #   },
  #   format = "file"
  # ),
  tar_target(
    tr_bar_all_df,
    generate_tr_bar_all_df(total_uncertainty_combined_df, sarea_uncertainty_combined_df, dir_dep_uncertainty_combined_df)
  ),
  tar_target(
    tr_bar_all_csv,
    tr_bar_all_df |> my_write_csv("data/tr_bar_all.csv"),
    format = "file"
  ),
  # tar_target(
  #   tr_bar_all_plot, {
  #     p <- tr_bar(tr_bar_all_df, id, tr_m, fill = id, group = id) +
  #       scale_y_continuous(breaks = c(0, 250, 500, 750, 1000, 1250)) # Override y-axis breaks for this plot
  #     my_ggsave(
  #       "figs/tr_bar_all",
  #       p,
  #       dpi = 600,
  #       width = 6.81,
  #       height = 4.4
  #     )
  #   },
  #   format = "file"
  # ),
  tar_target(
    rel_bar_df,
    generate_rel_cont_df(total_uncertainty_combined_df,
      ab_scaled_df_segments_xylem_post_ab_fit_draws_segments_xylem_0.08,
      sarea_uncertainty_combined_df,
      dir_dep_uncertainty_combined_df)
  ),
  tar_target(
    rel_bar_csv,
    my_write_csv(rel_bar_df, "data/rel_bar.csv"),
    format = "file"
  ),
  tar_target(
    rel_bar_plot, {
      p <- rel_bar(rel_bar_df)
      my_ggsave(
        "figs/rel_bar",
        p,
        dpi = 600,
        width = 4.33,
        height = 4.33
      )
    },
    format = "file"
  ),
  tar_target(
    tr_bar_comb_plot, {
      tr_bar_all_df2 <- tr_bar_all_df |>
        mutate(id = ifelse(id == "sapwood_aera", "sarea", id))
      id_fct <- c(tr_bar_all_df2$id, as.character(rel_bar_df$id)) |>
        unique() #|>
        # as.factor()
      id_fct <- factor(id_fct, levels = c("dir_only", "dir_dep", "dir_dep_cov", "dep_only", "sarea", "ab", "total"))

      tr_bar_all_df2 <- tr_bar_all_df2 |>
        mutate(id = factor(id, levels = c("dir_only", "dir_dep", "dep_only", "sarea", "total")))
      tmp <- tr_bar_all_df2 |> sample_n(2) |> mutate(id = factor(c("dir_dep_cov", "ab"))) |>
        mutate(across(where(is.numeric), ~ NA_real_))
      # tr_bar_all_df2 <- bind_rows(tr_bar_all_df2, tmp) |>
      #   mutate(id = factor(id, levels = id_fct))

      rel_bar_df2 <- rel_bar_df |>
        mutate(id = factor(id, levels = id_fct))
      tmp <- rel_bar_df2 |> sample_n(1) |> mutate(id = factor(c("total"))) |>
        mutate(across(where(is.numeric), ~ NA_real_))
      rel_bar_df2 <- bind_rows(rel_bar_df2, tmp) |>
        mutate(id = factor(id, levels = c("dir_only", "dir_dep_cov", "dep_only", "sarea", "ab")))

      p1 <- tr_bar_ab(tr_bar_ab_df |> filter(pg == "0.08" | model == "Model 4"),  ab_granier_uncertainty_combined_df, pg, tr_m, fill = model, group = model, model4 = FALSE) +
        coord_cartesian(ylim = c(0, 1250))  +
        scale_y_continuous(breaks = c(0, 250, 500, 750, 1000, 1250)) +
        # annotate("text", x = 1, y = 1250, label = "B", hjust = -1, vjust = 0.5, size = 5) +
        labs(fill = "A) Model") +
        theme(
          axis.text.x = element_text(size = 8, margin = margin(t = 0.5, r = 0, b = 0, l = 0)),
          # axis.text.y = element_text(size = 8,margin = margin(t = 0, r = 0.5, b = 0, l = 0)),
          # axis.title.y.left = element_blank(),
          # axis.text.y.left = element_blank(),
          legend.text = element_text(size = 8),
          legend.key.size = unit(0.4, "cm"),
          # legend.position = c(0.8, 0.25))
          legend.position = "right")
          #  legend.position = "top")
      p3 <- tr_bar(tr_bar_all_df2, ab_granier_uncertainty_combined_df, id, tr_m, fill = id, group = id) +
        scale_y_continuous(breaks = c(0, 250, 500, 750, 1000, 1250)) +
        # annotate("text", x = 1, y = 1250, label = "B", hjust = -0.1, vjust = 0.5, size = 5) +
        xlab("Source of uncertainty") +
        theme(
          axis.text.x = element_blank()
        )
      p4 <- rel_bar(rel_bar_df) +
        xlab("Source of uncertainty") +
        # annotate("text", x = 1, y = 60, label = "B", hjust = -0.1, vjust = 0.5, size = 5) +
        # scale_y_continuous(
        #      sec.axis = sec_axis(~ ., name = "Relative contribution (%)")) +
        theme(
        #   axis.title.y.right = element_text(angle = 90)
          # axis.title.y.left = element_blank(),
          legend.text = element_text(size = 8),
          legend.key.size = unit(0.4, "cm"),
          axis.title.x = element_blank()#,
          # legend.position = "right"
        )

# Define the layout with specific widths
    p1 <- p1 + plot_layout(widths = c(3))
    p3p4 <- p3 + p4 + plot_layout(widths = c(3, 1))

# Combine your plots
    p <- p1 / p3p4 +
      plot_layout(guides = 'collect') +
      plot_annotation(tag_levels = 'A')

      # Combi
      # p <- p1 / (p3 + p4) +
      #   plot_layout(guides = "collect")

      # p <- p + plot_layout(widths = c(3, 1)) +
      #   plot_annotation(tag_levels = "A") &
      #   theme(
      #     plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines")
      #     )

      my_ggsave(
        "figs/tr_bar_comb",
        p,
        dpi = 600,
        width = 6.81,
        height = 6.81
      )
    },
    format = "file"
  ),
  tar_target(
    ab_uncertainty_points_plot, {
      tmp1 <- ab_summarized_df_species_xylem_post_ab_fit_draws_species_xylem_0.08
      tmp2 <- segments_xylem_post_ab_fit_draws_segments_xylem_0.08 %>%
        mutate(id = 1:nrow(.)) |>
        mutate(id = as.character(id))
      tmp3 <- full_join(tmp1, tmp2)
      p <- ggplot(tmp3 |> sample_n(1000), aes(x = exp(log_a), y = b, col = tr)) +
        geom_point(alpha = 0.8) +
        scale_color_viridis_c(name = expression("Transpiration (mm"~y^-1*")" )) +
        scale_x_log10() +
        xlab(expression("Coefficient"~italic(a))) +
        ylab(expression("Coefficient"~italic(b))) +
        my_theme() +
        theme(
          legend.position = "right"
        )
      my_ggsave(
        "figs/ab_uncertainty_points",
        p,
        dpi = 600,
        width = 110,
        height = 65,
        units = "mm"
      )
    },
    format = "file"
  ),
  tar_target(
    tr_example_list,
    generate_tr_example_list(
      post_dir_dep_mid, dir_dep_imp_df, sarea_df,
      segments_xylem_post_ab_fit_draws_segments_xylem_0.02,
      segments_xylem_post_ab_fit_draws_segments_xylem_0.08)
  ),
  tar_target(
    tr_example_panel_plot, {
      p <- tr_example_panel(tr_example_list)
      my_ggsave(
        "figs/tr_example_panel_plot",
        p,
        dpi = 200,
        width = 173,
        height = 58,
        units = "mm"
      )
    },
    format = "file"
  ),
  tar_target(
    full_df_ab_processed,
    process_full_ab_df(sarea_df, post_dir_dep_mid, dir_dep_imp_df)
  ),
  tar_target(
    summary_stats_list,
    prepare_summary_stats(
      full_df_ab_processed,
      segments_xylem_post_ab_fit_draws_segments_xylem_0.02,
      segments_xylem_post_ab_fit_draws_segments_xylem_0.08)
  ),
  tar_target(
    ab_example_panel_plot, {
      p <- ab_example(
        full_df_ab_processed,
        segments_xylem_post_ab_fit_draws_segments_xylem_0.02,
        segments_xylem_post_ab_fit_draws_segments_xylem_0.08,
        summary_stats_list,
        ab_summarized_df_species_xylem_post_ab_fit_draws_species_xylem_0.08,
        segments_xylem_post_ab_fit_draws_segments_xylem_0.08
        )
      my_ggsave(
        "figs/ab_example_panel",
        p,
        dpi = 600,
        width = 173,
        height = 123,
        units = "mm"
      )
    },
    format = "file"
  ),
  tar_target(
    tr_ab_lm, {
      ab_summarized_df <- ab_summarized_df_species_xylem_post_ab_fit_draws_species_xylem_0.08
      segments_post <- segments_xylem_post_ab_fit_draws_segments_xylem_0.08

      segments_post2 <- segments_post %>%
        mutate(id = 1:nrow(.)) |>
        mutate(id = as.character(id))
      tmp <- full_join(ab_summarized_df, segments_post2)

      lm(tr ~ log_a + b, tmp)
    }
  ),
  tar_target(
    dbh_points_plot, {
      p <- dbh_points(dbh_imp_df2, girth_increment_csv)
      my_ggsave(
        "figs/dbh_points",
        p,
        dpi = 600,
        # width = 6.81,
        # height = 6.81
        width = 4.33,
        height = 4.33
      )
    },
    format = "file"
  ),
  tar_target(
    sap_dbh_points_plot, {
      p <- sap_dbh_points(sapwood_depth_csv, fit_draws_sap_dbh)
      my_ggsave(
        "figs/sap_dbh_points",
        p,
        dpi = 600,
        width = 4.33,
        height = 4.33
      )
    },
    format = "file"
  ),
  tar_target(
    imp_points_plot, {
      p <- imp_points(imputed_df_btrans_2016_2, rubber_raw_data_csv, year_1 = 2016, month_1 = 2, day_1 = 12,
                      imputed_df_btrans_2016_5, rubber_raw_data_csv, year_2 = 2016, month_2 = 5, day_2 = 7)
      my_ggsave(
        "figs/imp_points",
        p,
        dpi = 600,
        width = 6.81,
        height = 6.81
      )
    },
    format = "file"
  ),
  tar_target(
    dir_dep_table,
    write_dir_dep_table(fit2_summary_dir_dep, "data/dir_dep_post.csv"),
    format = "file"
  ),
  tar_target(
    dir_dep_comb_table,
    write_dir_dep_comb_table(post_dir_dep, "data/dir_dep_comb_post.csv"),
    format = "file"
  ),
  tar_target(
    sap_table,
    write_sap_table(fit_summary_sap_dbh,
       "data/dbh_sapwood_post.csv"),
    format = "file"
  ),
  NULL
)

sapwood_list <- list(
  tar_target(
    dbh_sap_stan_data,
    generate_dbh_sap_stan_data(sapwood_depth_csv)
  ),
  tar_target(
    dir_dep_stan_data,
    generate_dir_dep_stan_data(imputed_full_df)
  ),
  tar_stan_mcmc(
     fit,
     "stan/sap_dbh.stan",
     data = dbh_sap_stan_data,
     refresh = 0,
     chains = 4,
     parallel_chains = 4,
     iter_warmup = 2000,
     iter_sampling = 2000,
     adapt_delta = 0.9,
     max_treedepth = 15,
     seed = 123,
     return_draws = TRUE,
     return_diagnostics = TRUE,
     return_summary = TRUE,
     summaries = list(
       mean = ~mean(.x),
       sd = ~sd(.x),
       mad = ~mad(.x),
       ~posterior::quantile2(.x, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)),
       posterior::default_convergence_measures()
    )
  ),
  tar_stan_mcmc(
    fit2,
    "stan/dir_dep.stan",
    data = dir_dep_stan_data,
    refresh = 0,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 2000,
    iter_sampling = 2000,
    adapt_delta = 0.9,
    max_treedepth = 15,
    seed = 123,
    return_draws = TRUE,
    return_diagnostics = TRUE,
    return_summary = TRUE,
    summaries = list(
      mean = ~mean(.x),
      sd = ~sd(.x),
      mad = ~mad(.x),
    ~posterior::quantile2(.x, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)),
      posterior::default_convergence_measures()
   )
  ),
  tar_target(
    dbh_imp_df,
    generate_dbh_imp_data(
      girth_increment_csv,
      initial_dbh_csv)
  ),
  tar_target(
    dbh_imp_df2,
    generate_dbh_imp_data2(
      girth_increment_csv,
      initial_dbh_csv)
  ),
  NULL
)

append(raw_data_list, main_list) |>
  append(granier_list) |>
  append(tar_impute) |>
  append(sapwood_list) |>
  append(tar_dir_dep) |>
  append(uncertainty_list) |>
  append(uncertainty_figs_list) |>
  append(like_check_list)
