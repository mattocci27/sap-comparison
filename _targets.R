library(targets)
library(tarchetypes)
library(tidyverse)
library(stantargets)
library(cmdstanr)
library(furrr)
library(languageserver)
library(clustermq)
library(quarto)
library(bayesplot)

source("R/data_clean.R")
source("R/stan.R")
source("R/figs.R")

plan(multicore)
options(clustermq.scheduler = "multicore")

tar_option_set(packages = c(
  "tidyverse",
  "patchwork",
  "bayesplot",
  "httpgd",
  "smatr",
  "ggsma",
  "ggpubr",
  "ggridges",
  "RColorBrewer",
  "scales",
  "loo"
))

tar_option_set(
  garbage_collection = TRUE,
  memory = "transient"
)

# check if it's inside a container
# if (file.exists("/.dockerenv") | file.exists("/.singularity.d/startscript")) {
#   Sys.setenv(CMDSTAN = "/opt/cmdstan/cmdstan-2.29.2")
#   set_cmdstan_path("/opt/cmdstan/cmdstan-2.29.2")
# }

# cmdstan_version()

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
    anova_mvn_data,
    generate_anova_mvn_data(ks_trees_csv, model_type = "normal")
  ),
  tar_target(
    anova_lmvn_data,
    generate_anova_mvn_data(ks_trees_csv, model_type = "log-normal")
  ),
  tar_target(
    anova_gmvn_data,
    generate_anova_mvn_data(ks_trees_csv, model_type = "gamma")
  ),
  tar_target(
    anova_data,
    generate_anova_data(ks_spp_err_csv)
  ),
  tar_target(
    anova_data_log,
    generate_anova_data(ks_spp_err_csv, log = TRUE)
  ),
  tar_target(
    anova_data_err,
    generate_anova_data(ks_spp_err_csv, err = TRUE)
  ),
  tar_target(
    dummy_data,
    generate_dummy_data(n = 3, mu_hat = -2, n_beta = 3, sigma_alpha = .8, sigma_beta = 2, sigma_gamma = 0.1, sigma = 0.3, seed = 123)
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
     seed = 123
  ),
  tar_stan_mcmc(
     fit_anova_raw,
     c(
      "stan/anova.stan",
      "stan/anova_noint.stan",
      "stan/anova_gamma.stan"
      ),
     data = anova_data,
     refresh = 0,
     chains = 4,
     parallel_chains = getOption("mc.cores", 1),
     iter_warmup = 1000,
     iter_sampling = 1000,
     adapt_delta = 0.95,
     max_treedepth = 15,
     seed = 123
  ),
  tar_stan_mcmc(
     fit_anova,
     c(
      "stan/anova_int_err.stan",
      "stan/anova_noint_err.stan"
      ),
     data = anova_data_err,
     refresh = 0,
     chains = 4,
     parallel_chains = getOption("mc.cores", 1),
     iter_warmup = 1000,
     iter_sampling = 1000,
     adapt_delta = 0.95,
     max_treedepth = 15,
     seed = 123
  ),
  tar_stan_mcmc(
     fit_anova_log,
     c("stan/anova.stan",
      "stan/anova_noint.stan"
      ),
     data = anova_data_log,
     refresh = 0,
     chains = 4,
     parallel_chains = getOption("mc.cores", 1),
     iter_warmup = 1000,
     iter_sampling = 1000,
     adapt_delta = 0.95,
     max_treedepth = 15,
     seed = 123
  ),


  # tar_stan_mcmc(
  #    fit_mvn,
  #    "stan/anova_mvn.stan",
  #    data = anova_mvn_data,
  #    refresh = 0,
  #    chains = 4,
  #    parallel_chains = getOption("mc.cores", 4),
  #    iter_warmup = 1000,
  #    iter_sampling = 1000,
  #    adapt_delta = 0.9,
  #    max_treedepth = 15,
  #    seed = 123
  # ),
  # tar_stan_mcmc(
  #    fit_lmvn,
  #    "stan/anova_mvn.stan",
  #    data = anova_lmvn_data,
  #    refresh = 0,
  #    chains = 4,
  #    parallel_chains = getOption("mc.cores", 4),
  #    iter_warmup = 1000,
  #    iter_sampling = 1000,
  #    adapt_delta = 0.9,
  #    max_treedepth = 15,
  #    seed = 123
  # ),
  # tar_stan_mcmc(
  #    fit_gmvn,
  #    "stan/anova_mvn.stan",
  #    data = anova_gmvn_data,
  #    refresh = 0,
  #    chains = 4,
  #    parallel_chains = getOption("mc.cores", 4),
  #    iter_warmup = 1000,
  #    iter_sampling = 1000,
  #    adapt_delta = 0.9,
  #    max_treedepth = 15,
  #    seed = 123
  # ),

  tar_target(
    loo_,
    lapply(
      list(
           fit_anova_raw_mcmc_anova,
           fit_anova_raw_mcmc_anova_noint,
           fit_anova_log_mcmc_anova,
           fit_anova_log_mcmc_anova_noint
        ),
    \(x)x$loo(cores = parallel::detectCores())
    )
  ),

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

  tar_target(
    coef_intervals_sd_plot, {
      p <- coef_intervals_sd(fit_anova_draws_anova_noint_err)
      my_ggsave(
        "figs/coef_intervals_sd",
        p,
        dpi = 300,
        width = 8.5,
        height = 8.5,
        units = "cm"
      )
    },
    format = "file"
  ),

  tar_target(
    coef_intervals_mean_plot, {
      p <- coef_intervals_mean(fit_anova_draws_anova_noint_err)
      my_ggsave(
        "figs/coef_intervals_mean",
        p,
        dpi = 300,
        width = 8.5,
        height = 8.5,
        units = "cm"
      )
    },
    format = "file"
  ),

  tar_target(
    coef_intervals_diff_plot, {
      p <- coef_intervals_diff(fit_anova_draws_anova_noint_err)
      my_ggsave(
        "figs/coef_intervals_diff",
        p,
        dpi = 300,
        width = 8.5,
        height = 8.5,
        units = "cm"
      )
    },
    format = "file"
  ),

  tar_target(
    coef_density_sd_plot, {
      p <- coef_density_sd(fit_anova_draws_anova_noint_err)
      my_ggsave(
        "figs/coef_density_sd",
        p,
        dpi = 300,
        width = 8.5,
        height = 8.5,
        units = "cm"
      )
    },
    format = "file"
  ),

  tar_target(
    anova_yml,
    write_anova_yml(
      "yml/anova.yml",
      fit_anova_draws_anova_noint_err,
      ll = 0.25, hh = 0.75),
format = "file"
  ),
  # tar_target(
  #   ks_pred_draws,
  #   create_stan_tab(fit1_draws_anova)
  # ),
  # tar_target(
  #   ks_log_pred_draws,
  #   create_stan_tab(fit2_draws_anova)
  # ),

  # tar_target(
  #   ks_bar_plot,
  #   ks_bars(ks_five_spp_csv, ks_pred_draws)
  # ),
  # tar_target(
  #   ks_log_bar_plot,
  #   ks_bars(ks_five_spp_csv, ks_log_pred_draws)
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
        height = 16.2,
        units = "cm"
      )
    },
    format = "file"
  ),
  tar_target(
    logistic_sp_plot_linear, {
      p <- plot_logistic_sp(simple_logistic_draws_hierarchical_logistic, "data/cond_count.csv", quad = FALSE)
      my_ggsave(
        "figs/count_pressure_simple",
        p,
        dpi = 300,
        width = 4.3,
        height = 16.2,
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
        dpi = 300,
        width = 8.5,
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

  # tar_stan_mcmc(
  #    fit_dummy_sap,
  #    "stan/sap.stan",
  #    data = dummy_sap_stan_data,
  #    refresh = 0,
  #    chains = 4,
  #    parallel_chains = getOption("mc.cores", 4),
  #    iter_warmup = 2000,
  #    iter_sampling = 2000,
  #    adapt_delta = 0.9,
  #    max_treedepth = 15,
  #    seed = 123,
  #    output_dir = "loa"
  # ),


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

  tar_map(
    list(trait_set = c("all", "vaf", "vf", "ks", "noks", "noks2")),
    tar_target(sap_trait_clean,
      generate_sap_traits_stan_data(fd_k_traits_csv,
                             remove_abnormal_values = TRUE,
                             trait_set = trait_set))
  ),

  tar_map(
    values = list(stan_data = rlang::syms(c(
      "sap_trait_clean_all",
      "sap_trait_clean_vaf",
      "sap_trait_clean_vf",
      "sap_trait_clean_ks",
      "sap_trait_clean_noks",
      "sap_trait_clean_noks2"
      ))),
    tar_stan_mcmc(
      fit_abt,
      "stan/granier_with_traits.stan",
      data = stan_data,
      refresh = 0,
      chains = 4,
      parallel_chains = getOption("mc.cores", 2),
      iter_warmup = 2000,
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
    )
  ),

  tar_map(
    values = list(
      mcmc = rlang::syms(c(
        "fit_abt_mcmc_granier_with_traits_sap_trait_clean_all",
        "fit_abt_mcmc_granier_with_traits_sap_trait_clean_vaf",
        "fit_abt_mcmc_granier_with_traits_sap_trait_clean_vf",
        "fit_abt_mcmc_granier_with_traits_sap_trait_clean_ks",
        "fit_abt_mcmc_granier_with_traits_sap_trait_clean_noks",
        "fit_abt_mcmc_granier_with_traits_sap_trait_clean_noks2"))),
    tar_target(
      traits_loo,
      my_loo(mcmc)
    )
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

  tar_map(
    list(p = c(0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.015, 0.025, 0.035)),
    tar_target(sap_all_raw,
      generate_sap_stan_data(fd_k_traits_csv,
        upper_pressure = p)),
    tar_target(sap_all_clean,
      generate_sap_stan_data(fd_k_traits_csv,
        remove_abnormal_values = TRUE,
        upper_pressure = p))
  ),

  tar_map(
    values = list(stan_data = rlang::syms(c(
      "sap_all_raw_0.02",
      "sap_all_raw_0.03",
      "sap_all_raw_0.04",
      "sap_all_raw_0.05",
      "sap_all_raw_0.06",
      "sap_all_raw_0.07",
      "sap_all_raw_0.08",
      "sap_all_clean_0.015",
      "sap_all_clean_0.02",
      "sap_all_clean_0.025",
      "sap_all_clean_0.03",
      "sap_all_clean_0.035",
      "sap_all_clean_0.04",
      "sap_all_clean_0.05",
      "sap_all_clean_0.06",
      "sap_all_clean_0.07",
      "sap_all_clean_0.08"
      ))),
    tar_stan_mcmc(
      fit_ab,
      c("stan/granier_without_traits.stan",
       "stan/granier_without_traits2.stan"),
      data = stan_data,
      refresh = 0,
      chains = 4,
      parallel_chains = getOption("mc.cores", 1),
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

  tar_map(
    list(p = c(0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.015, 0.025, 0.035)),
    tar_target(sap_sp_raw,
      generate_sap_stan_data_sp(fd_k_traits_csv,
        upper_pressure = p)),
    tar_target(sap_sp_clean,
      generate_sap_stan_data_sp(fd_k_traits_csv,
        remove_abnormal_values = TRUE,
        upper_pressure = p))
  ),

  # # tar_target(
  # #   sap_stan_data_each,
  # #   generate_sap_stan_data_each(calibration_raw_data_csv)
  # # ),

  tar_map(
    values = list(stan_data_each = rlang::syms(c(
      "sap_sp_raw_0.02",
      "sap_sp_raw_0.03",
      "sap_sp_raw_0.04",
      "sap_sp_raw_0.05",
      "sap_sp_raw_0.06",
      "sap_sp_raw_0.07",
      "sap_sp_raw_0.08",
      "sap_sp_clean_0.015",
      "sap_sp_clean_0.02",
      "sap_sp_clean_0.025",
      "sap_sp_clean_0.03",
      "sap_sp_clean_0.035",
      "sap_sp_clean_0.04",
      "sap_sp_clean_0.05",
      "sap_sp_clean_0.06",
      "sap_sp_clean_0.07",
      "sap_sp_clean_0.08"
      ))),
    tar_target(
      fit_ab_each, {
      stan_data_each |>
          mutate(fit = map(stan_data, fit_model,
            granier_simple_file,
            iter_warmup = 2000,
            iter_sampling = 2000)) |>
          mutate(fit2 = map(stan_data, fit_model,
            granier_simple2_file,
            iter_warmup = 2000,
            iter_sampling = 2000))
      })
  ),

  tar_target(
    granier_simple_file,
    compile_model("stan/granier_simple.stan"),
    format = "file"
  ),

  tar_target(
    granier_simple2_file,
    compile_model("stan/granier_simple2.stan"),
    format = "file"
  ),

  tar_target(
    pool_multi_plot, {
      p <- line_pool_multi(fd_k_traits_csv,
       fit_ab_summary_granier_without_traits_sap_all_clean_0.08,
       fit_ab_summary_granier_without_traits2_sap_all_clean_0.08)
      my_ggsave(
        "figs/pool_multi",
        p,
        dpi = 300,
        width = 8,
        height = 12
      )
    },
    format = "file"
  ),

  tar_target(
    pg_multi_plot, {
      p <- line_pg_multi(
                         fd_k_traits_csv,
        xylem_lab,
                        k_range,
        fit_ab_summary_granier_without_traits_sap_all_clean_0.02,
        fit_ab_summary_granier_without_traits_sap_all_clean_0.025,
        fit_ab_summary_granier_without_traits_sap_all_clean_0.03,
        fit_ab_summary_granier_without_traits_sap_all_clean_0.035,
        fit_ab_summary_granier_without_traits_sap_all_clean_0.04,
        fit_ab_summary_granier_without_traits_sap_all_clean_0.05,
        fit_ab_summary_granier_without_traits_sap_all_clean_0.06,
        fit_ab_summary_granier_without_traits_sap_all_clean_0.07,
        fit_ab_summary_granier_without_traits_sap_all_clean_0.08)
      my_ggsave(
        "figs/pg_multi",
        p,
        dpi = 300,
        width = 8,
        height = 12
      )
    },
    format = "file"
  ),

  # tar_target(
  #   coef_density_plot, {
  #     p <- coef_density(xylem_lab,
  #       fit_ab_draws_granier_without_traits_sap_all_clean_0.08)
  #     my_ggsave(
  #       "figs/coef_density",
  #       p,
  #       dpi = 300,
  #       width = 12,
  #       height = 7
  #     )
  #   },
  #   format = "file"
  # ),

  tar_target(
    coef_density_plot2, {
      p <- coef_density(xylem_lab,
        fit_ab_draws_granier_without_traits_sap_all_clean_0.08,
        looks = "patchwork")
      my_ggsave(
        "figs/coef_density2",
        p,
        dpi = 300,
        width = 6.8,
        height = 5
      )
    },
    format = "file"
  ),

  tar_target(
    xylem_lab,
    generate_xylem_lab(fd_k_traits_csv)
  ),

  tar_target(
    pg_ribbon_a, {
      p <- ab_pg_ribbon(xylem_lab,
                        k_range,
        fit_ab_summary_granier_without_traits_sap_all_clean_0.02,
        fit_ab_summary_granier_without_traits_sap_all_clean_0.025,
        fit_ab_summary_granier_without_traits_sap_all_clean_0.03,
        fit_ab_summary_granier_without_traits_sap_all_clean_0.035,
        fit_ab_summary_granier_without_traits_sap_all_clean_0.04,
        fit_ab_summary_granier_without_traits_sap_all_clean_0.05,
        fit_ab_summary_granier_without_traits_sap_all_clean_0.06,
        fit_ab_summary_granier_without_traits_sap_all_clean_0.07,
        fit_ab_summary_granier_without_traits_sap_all_clean_0.08)
      my_ggsave(
        "figs/pg_ribbon_a",
        p,
        dpi = 300,
        width = 8,
        height = 12
      )
    },
    format = "file"
  ),
  tar_target(
    pg_ribbon_b, {
      p <- ab_pg_ribbon(xylem_lab,
                        k_range,
        fit_ab_summary_granier_without_traits_sap_all_clean_0.02,
        fit_ab_summary_granier_without_traits_sap_all_clean_0.025,
        fit_ab_summary_granier_without_traits_sap_all_clean_0.03,
        fit_ab_summary_granier_without_traits_sap_all_clean_0.035,
        fit_ab_summary_granier_without_traits_sap_all_clean_0.04,
        fit_ab_summary_granier_without_traits_sap_all_clean_0.05,
        fit_ab_summary_granier_without_traits_sap_all_clean_0.06,
        fit_ab_summary_granier_without_traits_sap_all_clean_0.07,
        fit_ab_summary_granier_without_traits_sap_all_clean_0.08,
        coef_a = FALSE)
      my_ggsave(
        "figs/pg_ribbon_b",
        p,
        dpi = 300,
        width = 8,
        height = 12
      )
    },
    format = "file"
  ),

  tar_target(
    k_range,
    generate_k_range(fd_k_traits_csv)
  ),

  tar_target(
    ab_var_clean_008,
    generate_ab_var_data(
      "data/ab_var_clean_008.csv",
      fit_ab_draws_granier_without_traits_sap_all_clean_0.08),
    format = "file"
  ),


  # tnar_map(
  #   list(p = c(0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08)),
  #   tar_target(sap_segment_raw,
  #     generate_sap_stan_data_segment(fd_k_traits_csv,
  #       upper_pressure = p)),
  #   tar_target(sap_segment_clean,
  #     generate_sap_stan_data_segment(fd_k_traits_csv,
  #       remove_abnormal_values = TRUE,
  #       upper_pressure = p))
  # ),

  # tar_target(
  #   fit_ab_each, {
  #     sap_stan_data_each |>
  #       mutate(fit = map(stan_data, fit_model, granier_simple_file))
  #   }
  # ),

  # tar_quarto(
  #   dummy_html,
  #   "docs/dummy.qmd"
  # ),
  # tar_quarto(
  #   ks_ratio_html,
  #   "docs/ks_ratio.qmd"
  # ),

  tar_quarto(
    data_exp_html,
    "docs/data_exp.qmd"
  ),
  tar_quarto(
    report_html,
    "docs/report.qmd"
  ),

  NULL
)

append(raw_data_list, main_list)
