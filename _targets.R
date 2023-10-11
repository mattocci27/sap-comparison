library(targets)
library(tarchetypes)
library(tidyverse)
library(stantargets)
library(cmdstanr)
library(furrr)
library(clustermq)
library(bayesplot)
# library(doParallel)

source("R/data_clean.R")
source("R/stan.R")
source("R/figs.R")
source("R/tables.R")
source("R/scale.R")

plan(multicore)
options(clustermq.scheduler = "multicore")

set.seed(123)

tar_option_set(packages = c(
  "tidyverse",
  "patchwork",
  "cowplot",
  "bayesplot",
  "httpgd",
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
  "ggrepel"
))

# tar_option_set(
#   garbage_collection = TRUE,
#   memory = "transient"
# )

pg <- c(seq(0.02, 0.08, by = 0.01), 0.025, 0.035)
# check if it's inside a container
# if (file.exists("/.dockerenv") | file.exists("/.singularity.d/startscript")) {
#   Sys.setenv(CMDSTAN = "/opt/cmdstan/cmdstan-2.29.2")
#   set_cmdstan_path("/opt/cmdstan/cmdstan-2.29.2")
# }

# cmdstan_version()

# Register the parallel backend
n_cores <- parallel::detectCores(logical = FALSE)  # Detect the number of available CPU cores
# cl <- parallel::makeCluster(n_cores - 1)  # Create a cluster with one less core than available
cl <- 9
doParallel::registerDoParallel(cl)  # Register the parallel backend

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
  # tar_target(
  #   ks_spp_err_csv,
  #   "data/ks_pres_tens_spp_err.csv",
  #   format = "file"
  # ),

  NULL
)

 # main analysis ----------------------------------
main_list <- list(
  tar_target(
    sma_scatter_plot,
    sma_scatter(five_spp_csv)
  ),
  tar_target(
    ks_box_plot,
    ks_box(ks_trees_csv)
  ),
  tar_target(
    ks_box_plot2,
    ks_box2(ks_trees_csv)
  ),
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
    values = list(log = c("log", "non-log")),
    tar_target(
      anova_data,
      generate_anova_data(ks_spp_err_csv, log = ifelse(log == "log", TRUE, FALSE))
    ),
    tar_stan_mcmc(
       fit,
       c("stan/anova.stan", "stan/anova_noint.stan"),
       data = anova_data,
       refresh = 0,
       chains = 1,
       parallel_chains = getOption("mc.cores", 1),
       iter_warmup = 1,
       iter_sampling = ,
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
     parallel_chains = getOption("mc.cores", 4),
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

  # tar_target(
  #   loo_,
  #   lapply(
  #     list(
  #          fit_anova_mcmc_anova_int_err,
  #          fit_anova_mcmc_anova_noint_err,
  #          fit_anova_log_mcmc_anova,
  #          fit_anova_log_mcmc_anova_noint
  #       ),
  #   \(x)x$loo(cores = parallel::detectCores())
  #   )
  # ),

  tar_target(
    sma_ks_plot, {
      p <- sma_ks(sma_scatter_log_plot, ks_box_plot)
      my_ggsave(
        "figs/sma_ks",
        p,
        dpi = 300,
        width = 6.81,
        height = 4.4
      )
    },
    format = "file"
  ),
  tar_target(
    sma_ks_plot2, {
      p <- sma_ks2(sma_scatter_log_plot, ks_box_plot2)
      my_ggsave(
        "figs/sma_ks2",
        p,
        dpi = 300,
        width = 6.81,
        height = 4.4
      )
    },
    format = "file"
  ),

  # tar_target(
  #   coef_intervals_sd_plot, {
  #     p <- coef_intervals_sd(fit_anova_draws_anova_noint_err)
  #     my_ggsave(
  #       "figs/coef_intervals_sd",
  #       p,
  #       dpi = 300,
  #       width = 8.5,
  #       height = 8.5,
  #       units = "cm"
  #     )
  #   },
  #   format = "file"
  # ),

  # tar_target(
  #   coef_intervals_mean_plot, {
  #     p <- coef_intervals_mean(fit_anova_draws_anova_noint_err)
  #     my_ggsave(
  #       "figs/coef_intervals_mean",
  #       p,
  #       dpi = 300,
  #       width = 8.5,
  #       height = 8.5,
  #       units = "cm"
  #     )
  #   },
  #   format = "file"
  # ),

  # tar_target(
  #   coef_intervals_pres_tens_plot, {
  #     p1 <- coef_intervals_sd(fit_anova_log_draws_anova_noint)
  #     p2 <- coef_intervals_mean(fit_anova_log_draws_anova_noint)
  #     p3 <- coef_intervals_diff(fit_anova_log_draws_anova_noint)
  #     p <- p1 + p2 + p3 + plot_spacer() +
  #       plot_layout(ncol = 2) +
  #       plot_annotation(tag_levels = "A")
  #     my_ggsave(
  #       "figs/coef_intervals_pres_tens",
  #       p,
  #       dpi = 300,
  #       width = 173,
  #       height = 173,
  #       units = "mm"
  #     )
  #   },
  #   format = "file"
  # ),

  # tar_target(
  #   coef_intervals_diff_plot, {
  #     p <- coef_intervals_diff(fit_anova_draws_anova_noint_err)
  #     my_ggsave(
  #       "figs/coef_intervals_diff",
  #       p,
  #       dpi = 300,
  #       width = 8.5,
  #       height = 8.5,
  #       units = "cm"
  #     )
  #   },
  #   format = "file"
  # ),

  # tar_target(
  #   coef_density_sd_plot, {
  #     p <- coef_density_sd(fit_anova_draws_anova_noint_err)
  #     my_ggsave(
  #       "figs/coef_density_sd",
  #       p,
  #       dpi = 300,
  #       width = 8.5,
  #       height = 8.5,
  #       units = "cm"
  #     )
  #   },
  #   format = "file"
  # ),

  # tar_target(
  #   anova_yml,
  #   write_anova_yml(
  #     "yml/anova.yml",
  #     fit_anova_log_draws_anova_noint,
  #     ll = 0.25, hh = 0.75),
  #   format = "file"
  # ),

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
    parallel_chains = getOption("mc.cores", 4),
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
    parallel_chains = getOption("mc.cores", 4),
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
    parallel_chains = getOption("mc.cores", 4),
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
        dpi = 300,
        width = 4.3,
        height = 14,
        units = "cm"
      )
    },
    format = "file"
  ),
  # tar_target(
  #   logistic_sp_plot_linear, {
  #     p <- plot_logistic_sp(simple_logistic_draws_hierarchical_logistic, "data/cond_count.csv", quad = FALSE)
  #     my_ggsave(
  #       "figs/count_pressure_simple",
  #       p,
  #       dpi = 300,
  #       width = 4.3,
  #       height = 16.2,
  #       units = "cm"
  #     )
  #   },
  #   format = "file"
  # ),

  tar_target(
    coef_intervals_logistic_plot, {
      p <- coef_intervals_logistic(quad_logistic_draws_hierarchical_logistic)
      my_ggsave(
        "figs/coef_intervals_logistic",
        p,
        dpi = 300,
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


  # tar_target(
  #   stan_data_no_xylem,
  #   generate_sap_traits_no_xylem_stan_data(fd_k_traits_csv,
  #                            remove_abnormal_values = TRUE)
  # ),
  # tar_stan_mcmc(
  #  fit_abt2,
  #  c("stan/granier_with_traits_no_xylem.stan",
  #  "stan/granier_with_traits_no_xylem_sp.stan"),
  #  data = stan_data_no_xylem,
  #  refresh = 0,
  #  chains = 4,
  #  parallel_chains = getOption("mc.cores", 4),
  #  iter_warmup = 1000,
  #  iter_sampling = 1000,
  #  adapt_delta = 0.9,
  #  max_treedepth = 15,
  #  seed = 123,
  #  return_draws = TRUE,
  #  return_diagnostics = TRUE,
  #  return_summary = TRUE,
  #  summaries = list(
  #    mean = ~mean(.x),
  #    sd = ~sd(.x),
  #    mad = ~mad(.x),
  #    ~posterior::quantile2(.x, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)),
  #    posterior::default_convergence_measures()
  #    )
  #  ),

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

  # tar_target(
  #   traits_loo,
  #   lapply(
  #     list(
  #           all = fit_abt_mcmc_granier_with_traits_sap_trait_clean_all,
  #           vaf = fit_abt_mcmc_granier_with_traits_sap_trait_clean_vaf,
  #           vf = fit_abt_mcmc_granier_with_traits_sap_trait_clean_vf),
  #   \(x)x$loo(cores = parallel::detectCores())
  #   )
  # ),

#   tar_map(
#     values = tibble(pg = c(0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.025, 0.035)) |>
#       mutate(post_ab_pool = c(paste0("fit_ab_draws_granier_without_traits_full_pool_sap_all_clean_", pg))) |>
#       mutate(post_ab_segments = str_replace_all(post_ab_pool, "pool", "segments")) |>
#       mutate(post_ab_pool = rlang::syms(post_ab_pool)) |>
#       mutate(post_ab_segments = rlang::syms(post_ab_segments)),
#     tar_target(
#       post_ab_pool_mc,
#       generate_post_ab(post_ab_pool) |> sample_n(1000)
#     ),
#     tar_target(
#       post_ab_segments_mc,
#       generate_post_ab(post_ab_segments) |> sample_n(1000)
#     )
#   ),




#   tar_target(
#     segments_ab_table_full, {
#       bind_rows(
#         segments_ab_table_sap_sp_clean_0.02,
#         segments_ab_table_sap_sp_clean_0.025,
#         segments_ab_table_sap_sp_clean_0.03,
#         segments_ab_table_sap_sp_clean_0.035,
#         segments_ab_table_sap_sp_clean_0.04,
#         segments_ab_table_sap_sp_clean_0.05,
#         segments_ab_table_sap_sp_clean_0.06,
#         segments_ab_table_sap_sp_clean_0.07,
#         segments_ab_table_sap_sp_clean_0.08
#       )
#     }
#   ),
#   tar_target(
#     pool_ab_table_full, {
#       bind_rows(
#         pool_ab_table_sap_sp_clean_0.02,
#         pool_ab_table_sap_sp_clean_0.025,
#         pool_ab_table_sap_sp_clean_0.03,
#         pool_ab_table_sap_sp_clean_0.035,
#         pool_ab_table_sap_sp_clean_0.04,
#         pool_ab_table_sap_sp_clean_0.05,
#         pool_ab_table_sap_sp_clean_0.06,
#         pool_ab_table_sap_sp_clean_0.07,
#         pool_ab_table_sap_sp_clean_0.08
#       )
#     }
#   ),
#   tar_target(
#     segments_inclusive_ab_csv,
#     generate_species_ab_table_csv(segments_ab_table_full, "data/segments_inclusive_ab.csv"),
#     format = "file"
#   ),
#   tar_target(
#     species_only_ab_csv,
#     generate_species_ab_table_csv(pool_ab_table_full, "data/species_only_ab.csv"),
#     format = "file"
#   ),

#   tar_target(
#     pool_multi_plot, {
#       p <- line_pool_multi(fd_k_traits_csv,
#        xylem_lab,
#        fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.08,
#        fit_ab_summary_granier_without_traits_full_pool_sap_all_clean_0.08)
#       my_ggsave(
#         "figs/pool_multi",
#         p,
#         dpi = 300,
#         width = 8,
#         height = 12
#       )
#     },
#     format = "file"
#   ),
#   tar_target(
#     ab_points_plot, {
#       p <- ab_comp_points(
#         pool_csv = without_traits_pool_0.08_csv,
#         seg_csv =
#         without_traits_fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.08_data.without_traits_segments_0.08.csv,
#         xylem_lab)
#       my_ggsave(
#         "figs/ab_points",
#         p,
#         dpi = 300,
#         width = 173,
#         height = 86,
#         units = "mm"
#       )
#     },
#     format = "file"
#   ),
#   tar_target(
#     ab_points_four_models_plot_all, {
#       p <- ab_comp_four_models_points(
#         fit_ab_each_sap_sp_clean_0.08,
#         fit_ab_summary_granier_without_traits_full_pool_sap_all_clean_0.08,
#         fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.08,
#         xylem_lab, rm_dip = FALSE)
#       my_ggsave(
#         "figs/ab_points_four_models_all",
#         p,
#         dpi = 300,
#         width = 110,
#         height = 165,
#         units = "mm"
#       )
#     },
#     format = "file"
#   ),
#   tar_target(
#     ab_points_four_models_plot, {
#       p <- ab_comp_four_models_points(
#         fit_ab_each_sap_sp_clean_0.08,
#         fit_ab_summary_granier_without_traits_full_pool_sap_all_clean_0.08,
#         fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.08,
#         xylem_lab, rm_dip = TRUE)
#       my_ggsave(
#         "figs/ab_points_four_models",
#         p,
#         dpi = 300,
#         width = 110,
#         height = 165,
#         units = "mm"
#       )
#     },
#     format = "file"
#   ),

#   tar_target(
#     pg_multi_plot, {
#       p <- line_pg_multi(
#         fd_k_traits_csv,
#         xylem_lab,
#         k_range,
#         fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.02,
#         fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.025,
#         fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.03,
#         fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.035,
#         fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.04,
#         fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.05,
#         fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.06,
#         fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.07,
#         fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.08)
#       my_ggsave(
#         "figs/pg_multi",
#         p,
#         dpi = 300,
#         width = 8,
#         height = 12
#       )
#     },
#     format = "file"
#   ),

#   tar_target(
#     vaf_pred_data_each,
#     generate_trait_fig_data(
#       fit_abt2_summary_granier_with_traits_no_xylem_log_vaf,
#       fit_abt2_draws_granier_with_traits_no_xylem_log_vaf,
#       fd_k_traits_csv,
#       xylem_lab,
#       "log_vaf",
#       no_xylem = TRUE,
#       single_trait = TRUE
#     )
#   ),
#   tar_target(
#     ks_pred_data_each,
#     generate_trait_fig_data(
#       fit_abt2_summary_granier_with_traits_no_xylem_log_ks,
#       fit_abt2_draws_granier_with_traits_no_xylem_log_ks,
#       fd_k_traits_csv,
#       xylem_lab,
#       "log_ks",
#       no_xylem = TRUE,
#       single_trait = TRUE
#     )
#   ),
#   tar_target(
#     dh_pred_data_each,
#     generate_trait_fig_data(
#       fit_abt2_summary_granier_with_traits_no_xylem_log_dh,
#       fit_abt2_draws_granier_with_traits_no_xylem_log_dh,
#       fd_k_traits_csv,
#       xylem_lab,
#       "log_dh",
#       no_xylem = TRUE,
#       single_trait = TRUE
#     )
#   ),
#   tar_target(
#     vf_pred_data_each,
#     generate_trait_fig_data(
#       fit_abt2_summary_granier_with_traits_no_xylem_log_vf,
#       fit_abt2_draws_granier_with_traits_no_xylem_log_vf,
#       fd_k_traits_csv,
#       xylem_lab,
#       "log_vf",
#       no_xylem = TRUE,
#       single_trait = TRUE
#     )
#   ),
#   tar_target(
#     wd_pred_data_each,
#     generate_trait_fig_data(
#       fit_abt2_summary_granier_with_traits_no_xylem_wood_density,
#       fit_abt2_draws_granier_with_traits_no_xylem_wood_density,
#       fd_k_traits_csv,
#       xylem_lab,
#       "wood_density",
#       no_xylem = TRUE,
#       single_trait = TRUE
#     )
#   ),

#   tar_target(
#     vaf_pred_data,
#     generate_trait_fig_data(
#       fit_abt_summary_granier_with_traits_sap_trait_clean_all,
#       fit_abt_draws_granier_with_traits_sap_trait_clean_all,
#       fd_k_traits_csv,
#       xylem_lab,
#       "log_vaf"
#     )
#   ),
#   tar_target(
#     ks_pred_data,
#     generate_trait_fig_data(
#       fit_abt_summary_granier_with_traits_sap_trait_clean_all,
#       fit_abt_draws_granier_with_traits_sap_trait_clean_all,
#       fd_k_traits_csv,
#       xylem_lab,
#       "log_ks"
#     )
#   ),
#   tar_target(
#     dh_pred_data,
#     generate_trait_fig_data(
#       fit_abt_summary_granier_with_traits_sap_trait_clean_all,
#       fit_abt_draws_granier_with_traits_sap_trait_clean_all,
#       fd_k_traits_csv,
#       xylem_lab,
#       "log_dh"
#     )
#   ),
#   tar_target(
#     wd_pred_data,
#     generate_trait_fig_data(
#       fit_abt_summary_granier_with_traits_sap_trait_clean_all,
#       fit_abt_draws_granier_with_traits_sap_trait_clean_all,
#       fd_k_traits_csv,
#       xylem_lab,
#       "wood_density"
#     )
#   ),
#   tar_target(
#     vf_pred_data,
#     generate_trait_fig_data(
#       fit_abt_summary_granier_with_traits_sap_trait_clean_all,
#       fit_abt_draws_granier_with_traits_sap_trait_clean_all,
#       fd_k_traits_csv,
#       xylem_lab,
#       "log_vf"
#     )
#   ),
#   tar_target(
#     vaf_pred_data_no_xylem,
#     generate_trait_fig_data(
#       fit_abt2_summary_granier_with_traits_no_xylem,
#       fit_abt2_draws_granier_with_traits_no_xylem,
#       fd_k_traits_csv,
#       xylem_lab,
#       "log_vaf",
#       no_xylem = TRUE
#     )
#   ),
#   tar_target(
#     ks_pred_data_no_xylem,
#     generate_trait_fig_data(
#       fit_abt2_summary_granier_with_traits_no_xylem,
#       fit_abt2_draws_granier_with_traits_no_xylem,
#       fd_k_traits_csv,
#       xylem_lab,
#       "log_ks",
#       no_xylem = TRUE
#     )
#   ),
#   tar_target(
#     dh_pred_data_no_xylem,
#     generate_trait_fig_data(
#       fit_abt2_summary_granier_with_traits_no_xylem,
#       fit_abt2_draws_granier_with_traits_no_xylem,
#       fd_k_traits_csv,
#       xylem_lab,
#       "log_dh",
#       no_xylem = TRUE
#     )
#   ),
#   tar_target(
#     wd_pred_data_no_xylem,
#     generate_trait_fig_data(
#       fit_abt2_summary_granier_with_traits_no_xylem,
#       fit_abt2_draws_granier_with_traits_no_xylem,
#       fd_k_traits_csv,
#       xylem_lab,
#       "wood_density",
#       no_xylem = TRUE
#     )
#   ),
#   tar_target(
#     vf_pred_data_no_xylem,
#     generate_trait_fig_data(
#       fit_abt2_summary_granier_with_traits_no_xylem,
#       fit_abt2_draws_granier_with_traits_no_xylem,
#       fd_k_traits_csv,
#       xylem_lab,
#       "log_vf",
#       no_xylem = TRUE
#     )
#   ),

#   # tar_target(
#   #   traits_points_plot, {
#   #     p <- traits_points(vaf_pred_data, ks_pred_data, xylem_lab)
#   #     my_ggsave(
#   #       "figs/traits_points",
#   #       p,
#   #       dpi = 300,
#   #       width = 7,
#   #       height = 7
#   #     )
#   #   },
#   #   format = "file"
#   # ),

#   tar_target(
#     traits_points_main_plot, {
#       p <- traits_points_main(
#               vaf_pred_data, log_vaf,
#               ks_pred_data, log_ks,
#               wd_pred_data, wood_density,
#               dh_pred_data, log_dh,
#               vf_pred_data, log_vf)
#       my_ggsave(
#         "figs/traits_points_main",
#         p,
#         dpi = 300,
#         # width = 6.8,
#         # height = 10
#         width = 4.33,
#         height = 10
#       )
#     },
#     format = "file"
#   ),
#   tar_target(
#     traits_points_all_si_plot, {
#       p <- traits_points_si(
#         vaf_pred_data_no_xylem, log_vaf,
#         ks_pred_data_no_xylem, log_ks,
#         wd_pred_data_no_xylem, wood_density,
#         dh_pred_data_no_xylem, log_dh,
#         vf_pred_data_no_xylem, log_vf,
#         use_color = FALSE)
#       my_ggsave(
#         "figs/traits_points_all_si",
#         p,
#         dpi = 300,
#         width = 4.33,
#         height = 10
#       )
#     },
#     format = "file"
#   ),
#   tar_target(
#     traits_points_si_plot, {
#       p <- traits_points_si(
#         vaf_pred_data_each, log_vaf,
#         ks_pred_data_each, log_ks,
#         wd_pred_data_each, wood_density,
#         dh_pred_data_each, log_dh,
#         vf_pred_data_each, log_vf,
#         use_color = FALSE)
#       my_ggsave(
#         "figs/traits_points_si",
#         p,
#         dpi = 300,
#         width = 6.8,
#         height = 2.7
#         # width = 4.33,
#         # height = 10
#       )
#     },
#     format = "file"
#   ),

#   tar_target(
#     coef_density_plot, {
#       p <- coef_density(xylem_lab,
#         fit_ab_draws_granier_without_traits_full_segments_sap_all_clean_0.08)
#       my_ggsave(
#         "figs/coef_density",
#         p,
#         dpi = 300,
#         width = 6.8,
#         height = 5
#       )
#     },
#     format = "file"
#   ),
#   tar_target(
#     xylem_lab,
#     generate_xylem_lab(fd_k_traits_csv)
#   ),

#   tar_target(
#     pg_ribbon_a, {
#       p <- ab_pg_ribbon(xylem_lab,
#                         k_range,
#         fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.02,
#         fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.025,
#         fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.03,
#         fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.035,
#         fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.04,
#         fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.05,
#         fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.06,
#         fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.07,
#         fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.08)
#       my_ggsave(
#         "figs/pg_ribbon_a",
#         p,
#         dpi = 300,
#         width = 8,
#         height = 12
#       )
#     },
#     format = "file"
#   ),
#   tar_target(
#     pg_ribbon_b, {
#       p <- ab_pg_ribbon(xylem_lab,
#                         k_range,
#         fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.02,
#         fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.025,
#         fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.03,
#         fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.035,
#         fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.04,
#         fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.05,
#         fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.06,
#         fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.07,
#         fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.08,
#         coef_a = FALSE)
#       my_ggsave(
#         "figs/pg_ribbon_b",
#         p,
#         dpi = 300,
#         width = 8,
#         height = 12
#       )
#     },
#     format = "file"
#   ),

#   tar_target(
#     k_range,
#     generate_k_range(fd_k_traits_csv)
#   ),

#   tar_target(
#     ab_var_clean_008,
#     generate_ab_var_data(
#       "data/ab_var_clean_008.csv",
#       fit_ab_draws_granier_without_traits_full_segments_sap_all_clean_0.08),
#     format = "file"
#   ),


#   tar_map(
#     list(p = c(0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08)),
#     tar_target(sap_segment_raw,
#       generate_sap_stan_data_segment(fd_k_traits_csv,
#         upper_pressure = p)),
#     tar_target(sap_segment_clean,
#       generate_sap_stan_data_segment(fd_k_traits_csv,
#         remove_abnormal_values = TRUE,
#         upper_pressure = p))
#   ),

#   tar_map(
#     list(stan_summary =
#       rlang::syms(
#       str_c("fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_",
#         c(seq(0.02, 0.08, by = 0.01), 0.025, 0.035))),
#       output = str_c("data/without_traits_segments_",
#         c(seq(0.02, 0.08, by = 0.01), 0.025, 0.035), ".csv")),
#     tar_target(without_traits,
#       write_without_traits_csv(stan_summary, output),
#      format = "file")
#   ),

#   tar_target(
#     without_traits_pool_0.08_csv,
#     write_without_traits_csv(
#       fit_ab_summary_granier_without_traits_full_pool_sap_all_clean_0.08,
#       "data/without_traits_pool_0.08.csv"),
#     format = "file"
#   ),

#   tar_target(
#     ab_csv,
#     write_ab_csv2(
#       fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_0.08,
#       fit_ab_summary_granier_without_traits_full_pool_sap_all_clean_0.08,
#       fit_ab_each_sap_sp_clean_0.08,
#       xylem_lab,
#       "data/all_ab.csv"),
#     format = "file"
#   ),
#   tar_target(
#     ks_seg_table,
#     generate_summary_trait_table(
#       fit_abt_summary_granier_with_traits_sap_trait_clean_ks,
#       fd_k_traits_csv
#     )
#   ),
#   tar_target(
#     vaf_seg_table,
#     generate_summary_trait_table(
#       fit_abt_summary_granier_with_traits_sap_trait_clean_vaf,
#       fd_k_traits_csv
#     )
#   ),
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
    # tar_target(sap_all_raw,
    #   generate_sap_stan_data(fd_k_traits_csv,
    #     upper_pressure = p)),
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
    )
  )

segments_xylem_df_mapped <- tar_map(
  list(stan_summary =
    rlang::syms(
    str_c("fit_summary_segments_xylem_",
      c(seq(0.02, 0.08, by = 0.01), 0.025, 0.035))),
    pg = pg),
  tar_target(
    table,
    generates_segments_xylem_table(stan_summary)
  ),
  tar_target(
    table_re,
    generate_summary_non_trait_table(table, pg = pg)
  )
)
species_xylem_df_mapped <- tar_map(
  list(stan_summary =
    rlang::syms(
    str_c("fit_summary_species_xylem_",
      c(seq(0.02, 0.08, by = 0.01), 0.025, 0.035))),
    pg = pg),
  tar_target(
    table,
    generates_segments_xylem_table(stan_summary)
  ),
  tar_target(
    table_re,
    generate_summary_non_trait_table(table, pg = pg)
  )
)

tar_combined_segments_xylem_df <- tar_combine(
  segments_xylem_df_combined,
  segments_xylem_df_mapped[["table_re"]],
  command = dplyr::bind_rows(!!!.x)
)
tar_combined_species_xylem_df <- tar_combine(
  species_xylem_df_combined,
  species_xylem_df_mapped[["table_re"]],
  command = dplyr::bind_rows(!!!.x)
)



# with single traits -----------------------------------------
  granier_with_traits_mapped <- tar_map(
    values = list(trait_name = rlang::syms(c(
      "log_dh", "log_vf", "wood_density",
      "log_ks", "log_vaf", "int"))),
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
      parallel_chains = getOption("mc.cores", 2),
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
      parallel_chains = getOption("mc.cores", 2),
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

granier_list <- list(
  granier_with_traits_mapped,
  granier_without_traits_mapped,
  tar_combined_species_ab_table,
  tar_combined_segments_ab_table,
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
    segments_xylem_csv,
    my_write_csv(segments_xylem_df_combined, "data/segments_xylem_post.csv"),
    format = "file"
  ),
  tar_target(
    species_xylem_csv,
    my_write_csv(species_xylem_df_combined, "data/species_xylem_post.csv"),
    format = "file"
  ),
  # tar_target(
  #   full_segments_traits_csv,
  #   my_write_csv(all_seg_table, "data/full_segments_traits_post.csv"),
  #   format = "file"
  # ),
  NULL
)

#------------------------------------------------------------

values <- expand_grid(year = c(2015, 2016), month = 1:12) |>
  filter(!(year == 2016 & (month == 12 | month == 10 | month == 9)))

values2 <- expand_grid(year = 2016, month = c(9, 10, 12)) |>
  mutate(df_name = rlang::syms(paste0("imputed_df_", year - 1, "_", month)))

impute_mapped <- tar_map(
  values = values,
  tar_target(
    imputed_data,
      missForest_long(
        csv = rubber_raw_data_csv,
        year = year,
        month = month
      )
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
    imputed_rest,
    {
      tmp <- missForest_clean(
        csv = rubber_raw_data_csv,
        year = year,
        month = month)
      bind_rows(tmp, df_name) |>
        missForest(parallelize = "forests")
    }
   ),
  tar_target(
    imputed_df2,
    clean_imputed_df(imputed_rest)
   ),
  tar_target(
    nramse2,
    imputed_rest$OOBerror[1]
    ),
   tar_target(
     imputed_df_btrans2,
     backtransform_date(rubber_raw_data_csv, year, month, imputed_df2)
   )
)

tar_combined_nrmse_rest <- tar_combine(
  combined_nrmse_rest,
  impute_rest_mapped[["nramse2"]],
  command = dplyr::bind_rows(!!!.x, .id = "id")
)

tar_combined_imputed_rest_data <- tar_combine(
  combined_imputed_rest_mapped,
  impute_rest_mapped[["imputed_df_btrans2"]],
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
      dplyr::select(year, date, time, vpd, par, ks, tree, dir, dep)
    }
  ),
  tar_target(
    nonimputed_full_df,
    make_long_nonimputed_df(rubber_raw_data_csv)
  ),
  NULL
  )

# tar_dir_dep <- list(
#   tar_target(
#     dir_dep_imp_df,
#     generate_dir_dep_imp_data(
#       imputed_full_df)
#   ),
#   # tar_target(
#   #   post_ab_pool,
#   #   generate_post_ab(fit_ab_draws_granier_without_traits_full_pool_sap_all_clean_0.08)
#   # ),
#   # tar_target(
#   #   post_ab_segments,
#   #   generate_post_ab(fit_ab_draws_granier_without_traits_full_segments_sap_all_clean_0.08)
#   # ),
#   tar_target(
#     post_dir_dep,
#     generate_post_dir_dep(
#       fit_dir_dep_draws_no_temporal_hourly_dir,
#       fit_dir_dep_draws_no_temporal_hourly_dep)
#   ),
#   tar_target(
#     post_slen, {
#       fit_dbh_sapwood_draws_normal |>
#         dplyr::select(alpha, beta, sigma)
#     }
#   ),
#   # tar_target(
#   #   post_ab_pool_mc,
#   #   post_ab_pool |> sample_n(1000)
#   # ),
#   # tar_target(
#   #   post_ab_segments_mc,
#   #   post_ab_segments |> sample_n(1000)
#   # ),
#   tar_target(
#     post_dir_dep_mc,
#     post_dir_dep |> mutate(beta = map(beta, sample, 1000))
#   ),
#   tar_target(
#     post_slen_mc,
#     post_slen |> sample_n(1000)
#   ),
#   tar_target(
#     dir_dep_imp_full_df,
#     add_t16(rubber_raw_data_csv, dir_dep_imp_df, post_dir_dep)
#   ),
#   NULL
# )

uncertainty_mapped <- tar_map(
    values = expand_grid(folds = 1:30,
      # pg = c(0.02, 0.03)) |>
      pg = c(0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.025, 0.035)) |>
      mutate(post_ab_pool_mc =
        paste0(
          "fit_ab_draws_granier_without_traits_full_pool_sap_all_clean_", pg
        )) |>
      mutate(post_ab_segments_mc = str_replace_all(post_ab_pool_mc, "pool", "segments")) |>
      mutate(post_ab_pool_mc = rlang::syms(post_ab_pool_mc)) |>
      mutate(post_ab_segments_mc = rlang::syms(post_ab_segments_mc)),
    tar_target(
      post_ab_pool_mc2, {
        set.seed(123)
        generate_post_ab(post_ab_pool_mc) |> sample_n(1000)
      }
    ),
    tar_target(
      post_ab_segments_mc2, {
        set.seed(123)
        generate_post_ab(post_ab_segments_mc) |> sample_n(1000)
      }
    ),
    tar_target(
      ab_uncertainty_df,
      generate_ab_uncertainty(
        dir_dep_imp_full_df,
        dbh_imp_df,
        post_ab_pool_mc = post_ab_pool_mc2,
        post_ab_segments_mc = post_ab_segments_mc2,
        post_slen, post_dir_dep, k = 30, i = folds) |>
        mutate(pg = pg)
    )
  )
uncertainty_granier_mapped <- tar_map(
    values = expand_grid(folds = 1:60,
      pg = c(0.08)) |>
      mutate(post_ab_pool_mc =
        paste0(
          "fit_ab_draws_granier_without_traits_full_pool_sap_all_clean_", pg
        )) |>
      mutate(post_ab_segments_mc = str_replace_all(post_ab_pool_mc, "pool", "segments")) |>
      mutate(post_ab_pool_mc = rlang::syms(post_ab_pool_mc)) |>
      mutate(post_ab_segments_mc = rlang::syms(post_ab_segments_mc)),
    tar_target(
      post_ab_pool_mc3, {
        set.seed(123)
        generate_post_ab(post_ab_pool_mc) |> sample_n(1000)
      }
    ),
    tar_target(
      post_ab_segments_mc3, {
        set.seed(123)
        generate_post_ab(post_ab_segments_mc) |> sample_n(1000)
      }
    ),
    tar_target(
      ab_uncertainty_granier_df,
      generate_ab_uncertainty_granier(
        dir_dep_imp_full_df,
        dbh_imp_df,
        post_ab_pool_mc = post_ab_pool_mc3,
        post_ab_segments_mc = post_ab_segments_mc3,
        post_slen, post_dir_dep, k = 60, i = folds) |>
        mutate(pg = pg)
    )
  )

uncertainty_each_mapped <- tar_map(
    values = tibble(folds = 1:60),
    tar_target(
      ab_uncertainty_df_each,
      generate_ab_uncertainty(
        dir_dep_imp_full_df,
        dbh_imp_df,
        post_ab_pool_mc = post_ab_pool_mc4,
        post_ab_segments_mc = post_ab_segments_mc4,
        post_slen, post_dir_dep, k = 60, i = folds) |>
        mutate(pg = 0.08)
    )
 )

tar_combined_ab_uncertainty <- tar_combine(
  ab_uncertainty_full_df,
  uncertainty_mapped[["ab_uncertainty_df"]],
  command = dplyr::bind_rows(!!!.x)
)
tar_combined_ab_uncertainty_granier <- tar_combine(
  ab_uncertainty_full_granier_df,
  uncertainty_granier_mapped[["ab_uncertainty_granier_df"]],
  command = dplyr::bind_rows(!!!.x)
)
tar_combined_ab_uncertainty_each <- tar_combine(
  ab_uncertainty_full_each_df,
  uncertainty_each_mapped[["ab_uncertainty_df_each"]],
  command = dplyr::bind_rows(!!!.x)
)

uncertainty_list <- list(
    tar_target(
      post_ab_pool_mc4, {
        set.seed(123)
        generate_post_ab_each(fit_ab_each_sap_sp_clean_0.08) |> sample_n(1000)
      }
    ),
    tar_target(
      post_ab_segments_mc4, {
        set.seed(123)
        generate_post_ab_each(fit_ab_each_sap_sp_clean_0.08) |> sample_n(1000)
      }
    ),
  uncertainty_mapped,
  uncertainty_granier_mapped,
  uncertainty_each_mapped,
  tar_combined_ab_uncertainty,
  tar_combined_ab_uncertainty_granier,
  tar_combined_ab_uncertainty_each,
  tar_target(
    tr_scaled_bars_plot, {
      p <- tr_scaled_bars(ab_uncertainty_full_df)
      my_ggsave(
        "figs/tr_scaled_bars",
        p,
        dpi = 300,
        width = 6.81,
        height = 4.4
      )
    },
    format = "file"
  ),
  tar_target(
    tr_scaled_bars_plot2, {
      p <- tr_scaled_bars2(ab_uncertainty_full_df, ab_uncertainty_full_each_df)
      my_ggsave(
        "figs/tr_scaled_bars2",
        p,
        dpi = 300,
        width = 6.81,
        height = 4.4
      )
    },
    format = "file"
  ),
  tar_target(
    tr_scaled_bars_csv,
      generate_tr_scaled_df(ab_uncertainty_full_df, ab_uncertainty_full_each_df) |>
        my_write_csv("data/ec_scaled.csv"),
    format = "file"
  ),
  tar_target(
    dbh_points_plot, {
      p <- dbh_points(dbh_imp_df2, girth_increment_csv)
      my_ggsave(
        "figs/dbh_points",
        p,
        dpi = 300,
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
      p <- sap_dbh_points(sapwood_depth_csv, fit_dbh_sapwood_draws_normal)
      my_ggsave(
        "figs/sap_dbh_points",
        p,
        dpi = 300,
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
        dpi = 300,
        width = 6.81,
        height = 6.81
      )
    },
    format = "file"
  ),
  tar_target(
   scaled_sapflow_csv, {
      tmp1 <- generate_tr_scaled_bars_data(ab_uncertainty_full_each_df, each = TRUE) |>
        mutate(model = paste("full", model, sep = "_"))
      tmp2 <- generate_tr_scaled_bars_data(ab_uncertainty_full_df) |>
        mutate(model = paste("sep", model, sep = "_"))
      tmp3 <- bind_rows(tmp1, tmp2) |>
        arrange(pg)
      my_write_csv(tmp3, "data/scaled_spflow.csv")
   },
   format = "file"
  ),
  tar_target(
    dir_dep_table,
    write_dir_dep_table(fit_dir_dep_summary_no_temporal_hourly_dir,
      fit_dir_dep_summary_no_temporal_hourly_dep, "data/dir_dep_post.csv"),
    format = "file"
  ),
  tar_target(
    sap_table,
    write_sap_table(fit_dbh_sapwood_summary_normal,
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
  tar_stan_mcmc(
     fit_dbh_sapwood,
     "stan/normal.stan",
     data = dbh_sap_stan_data,
     refresh = 0,
     chains = 4,
     parallel_chains = getOption("mc.cores", 4),
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
  tar_map(
    values = expand_grid(time_res = c("daily", "hourly"), fct = c("dir", "dep")),
    tar_target(
      dir_dep_stan_data,
      generate_dir_dep_stan_data(imputed_full_df, time_res = time_res, fct = fct)
    ),
    tar_stan_mcmc(
      fit_dir_dep,
      "stan/no_temporal.stan",
      data = dir_dep_stan_data,
      refresh = 1,
      chains = 4,
      parallel_chains = getOption("mc.cores", 4),
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

segments_xylem_df_mapped <- tar_map(
  list(stan_summary =
    rlang::syms(
    str_c("fit_ab_summary_granier_without_traits_full_segments_sap_all_clean_",
      c(seq(0.02, 0.08, by = 0.01), 0.025, 0.035))),
    pg = pg),
  tar_target(
    segments_xylem_table,
    generates_segments_xylem_table(stan_summary)
  ),
  tar_target(
    without_traits_table,
    generate_summary_non_trait_table(segments_xylem_table, pg = pg)
  )
)
species_xylem_df_mapped <- tar_map(
  list(stan_summary =
    rlang::syms(
    str_c("fit_ab_summary_granier_without_traits_full_pool_sap_all_clean_",
      c(seq(0.02, 0.08, by = 0.01), 0.025, 0.035))),
    pg = pg),
  tar_target(
    species_xylem_table,
    generates_segments_xylem_table(stan_summary)
  ),
  tar_target(
    without_traits_table,
    generate_summary_non_trait_table(species_xylem_table, pg = pg)
  )
)

 tar_combined_segments_xylem_df <- tar_combine(
   segments_xylem_df_combined,
   segments_xylem_df_mapped[["without_traits_table"]],
   command = dplyr::bind_rows(!!!.x)
  )
 tar_combined_species_xylem_df <- tar_combine(
   species_xylem_df_combined,
   species_xylem_df_mapped[["without_traits_table"]],
   command = dplyr::bind_rows(!!!.x)
  )

 post_csv_list <- list(
   segments_xylem_df_mapped,
   species_xylem_df_mapped,
   tar_combined_segments_xylem_df,
   tar_combined_species_xylem_df,
   tar_target(
     segments_xylem_csv,
     my_write_csv(segments_xylem_df_combined, "data/segments_xylem_post.csv"),
     format = "file"
   ),
   tar_target(
     species_xylem_csv,
     my_write_csv(species_xylem_df_combined, "data/species_xylem_post.csv"),
     format = "file"
   ),
   tar_target(
     full_segments_traits_csv,
     my_write_csv(all_seg_table, "data/full_segments_traits_post.csv"),
     format = "file"
   )
 )

# append(raw_data_list, main_list) |>
#   append(tar_impute) |>
#   append(sapwood_list) |>
#   append(tar_dir_dep) |>
#   append(uncertainty_list) |>
#   append(post_csv_list)
append(raw_data_list, main_list) |>
  append(granier_list) |>
  append(tar_impute) |>
  append(sapwood_list)
