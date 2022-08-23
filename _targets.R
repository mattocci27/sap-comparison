library(targets)
library(tarchetypes)
library(tidyverse)
library(stantargets)
library(cmdstanr)
library(furrr)
library(languageserver)
library(clustermq)
library(quarto)

source("R/data_clean.R")
source("R/stan.R")
source("R/figs.R")

plan(multicore)
options(clustermq.scheduler = "multicore")

tar_option_set(packages = c(
  "tidyverse",
  "httpgd",
  "ggsma",
  "ggpubr",
  "RColorBrewer"
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
    ks_trees_csv,
    clean_ks_trees(ks_five_trees_raw_csv),
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
    generate_anova_data(ks_five_spp_csv)
  ),
  tar_target(
    anova_data_log,
    generate_anova_data(ks_five_spp_csv, log = TRUE)
  ),
  tar_target(
    dummy_data,
    generate_dummy_data(n = 3, mu_hat = -2, n_beta = 3, sigma_alpha = .8, sigma_beta = 2, sigma_gamma = 0.1, sigma = 0.3, seed = 123)
  ),
  tar_stan_mcmc(
     fit_dummy,
     "stan/anova.stan",
     data = dummy_data,
     refresh = 0,
     chains = 4,
     parallel_chains = getOption("mc.cores", 4),
     iter_warmup = 1000,
     iter_sampling = 1000,
     draws = TRUE,
     diagnostics = TRUE,
     summary = TRUE,
     adapt_delta = 0.9,
     max_treedepth = 15,
     seed = 123
  ),
  tar_stan_mcmc(
     fit_dummy_noint,
     "stan/anova_noint.stan",
     data = dummy_data,
     refresh = 0,
     chains = 4,
     parallel_chains = getOption("mc.cores", 4),
     iter_warmup = 1000,
     iter_sampling = 1000,
     draws = TRUE,
     diagnostics = TRUE,
     summary = TRUE,
     adapt_delta = 0.9,
     max_treedepth = 15,
     seed = 123
  ),

  tar_stan_mcmc(
     fit1,
     "stan/anova.stan",
     data = anova_data,
     refresh = 0,
     chains = 4,
     parallel_chains = getOption("mc.cores", 4),
     iter_warmup = 1000,
     iter_sampling = 1000,
     draws = TRUE,
     diagnostics = TRUE,
     summary = TRUE,
     adapt_delta = 0.95,
     max_treedepth = 15,
     seed = 123
  ),

  tar_stan_mcmc(
     fit_anova_inter,
     "stan/anova.stan",
     data = anova_data,
     refresh = 0,
     chains = 4,
     parallel_chains = getOption("mc.cores", 4),
     iter_warmup = 1000,
     iter_sampling = 1000,
     draws = TRUE,
     diagnostics = TRUE,
     summary = TRUE,
     adapt_delta = 0.95,
     max_treedepth = 15,
     seed = 123
  ),
  tar_stan_mcmc(
     fit_anova_nointer,
     "stan/anova_noint.stan",
     data = anova_data,
     refresh = 0,
     chains = 4,
     parallel_chains = getOption("mc.cores", 4),
     iter_warmup = 1000,
     iter_sampling = 1000,
     draws = TRUE,
     diagnostics = TRUE,
     summary = TRUE,
     adapt_delta = 0.95,
     max_treedepth = 15,
     seed = 123
  ),
  tar_stan_mcmc(
     fit_anova_log_inter,
     "stan/anova.stan",
     data = anova_data_log,
     refresh = 0,
     chains = 4,
     parallel_chains = getOption("mc.cores", 4),
     iter_warmup = 1000,
     iter_sampling = 1000,
     draws = TRUE,
     diagnostics = TRUE,
     summary = TRUE,
     adapt_delta = 0.95,
     max_treedepth = 15,
     seed = 123
  ),
  tar_stan_mcmc(
     fit_anova_log_nointer,
     "stan/anova_noint.stan",
     data = anova_data_log,
     refresh = 0,
     chains = 4,
     parallel_chains = getOption("mc.cores", 4),
     iter_warmup = 1000,
     iter_sampling = 1000,
     draws = TRUE,
     diagnostics = TRUE,
     summary = TRUE,
     adapt_delta = 0.95,
     max_treedepth = 15,
     seed = 123
  ),


  tar_stan_mcmc(
     fit3,
     "stan/anova_gamma.stan",
     data = anova_data,
     refresh = 0,
     chains = 4,
     parallel_chains = getOption("mc.cores", 4),
     iter_warmup = 1000,
     iter_sampling = 1000,
     draws = TRUE,
     diagnostics = TRUE,
     summary = TRUE,
     adapt_delta = 0.95,
     max_treedepth = 15,
     seed = 123
  ),
  tar_stan_mcmc(
     fit_mvn,
     "stan/anova_mvn.stan",
     data = anova_mvn_data,
     refresh = 0,
     chains = 4,
     parallel_chains = getOption("mc.cores", 4),
     iter_warmup = 1000,
     iter_sampling = 1000,
     draws = TRUE,
     diagnostics = TRUE,
     summary = TRUE,
     adapt_delta = 0.9,
     max_treedepth = 15,
     seed = 123
  ),
  tar_stan_mcmc(
     fit_lmvn,
     "stan/anova_mvn.stan",
     data = anova_lmvn_data,
     refresh = 0,
     chains = 4,
     parallel_chains = getOption("mc.cores", 4),
     iter_warmup = 1000,
     iter_sampling = 1000,
     draws = TRUE,
     diagnostics = TRUE,
     summary = TRUE,
     adapt_delta = 0.9,
     max_treedepth = 15,
     seed = 123
  ),
  tar_stan_mcmc(
     fit_gmvn,
     "stan/anova_mvn.stan",
     data = anova_gmvn_data,
     refresh = 0,
     chains = 4,
     parallel_chains = getOption("mc.cores", 4),
     iter_warmup = 1000,
     iter_sampling = 1000,
     draws = TRUE,
     diagnostics = TRUE,
     summary = TRUE,
     adapt_delta = 0.9,
     max_treedepth = 15,
     seed = 123
  ),

  tar_target(
    loo_,
    lapply(
      list(
           fit_anova_inter_mcmc_anova,
           fit_anova_nointer_mcmc_anova_noint,
           fit_anova_log_inter_mcmc_anova,
           fit_anova_log_nointer_mcmc_anova_noint
        ),
    \(x)x$loo(cores = parallel::detectCores())
    )
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

  # tar_quarto(
  #   report_html,
  #   "docs/report.qmd"
  # ),

  tar_quarto(
    dummy_html,
    "docs/dummy.qmd"
  ),

  NULL

)

append(raw_data_list, main_list)
