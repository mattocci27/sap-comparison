library(targets)
library(tarchetypes)
library(tidyverse)
library(stantargets)
library(cmdstanr)
library(furrr)
library(languageserver)

source("R/stan.R")
source("R/figs.R")

plan(multicore)
options(clustermq.scheduler = "multicore")

tar_option_set(packages = c(
  "tidyverse",
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

list(
  # data cleaning ----------------------------------
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
    anova_data,
    generate_anova_data(five_spp_csv)
  ),
  tar_target(
    anova_data_log,
    generate_anova_data(five_spp_csv, log = TRUE)
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
     fit2,
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
  tar_quarto(
    report_html,
    "docs/report.qmd"
  )
)
