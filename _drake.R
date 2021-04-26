setwd(here::here())

source("R/packages.R")
source("R/functions.R")
source("R/settings.R")
source("R/plan.R")

options(clustermq.scheduler = "multicore")

plan <- bind_plans(qc_prep, pre_impute_qc, download_impute, post_impute, analysis_prep, analysis_prep_case_only, init_analysis, process_init)

drake_config(plan, verbose = 2)
