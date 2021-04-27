setwd(here::here())

source("R/packages.R")
source("R/functions.R")
source("R/settings.R")
source("R/plan.R")

options(clustermq.scheduler = "multicore")

plan <- bind_plans(qc_prep, pre_impute_qc, download_impute, post_impute, analysis_prep, analysis_prep_case_only, init_analysis, process_init)

drake_config(plan, verbose = 1, log_make = "plan.log", cache_log_file = "cache_log.csv",
  parallelism = "clustermq", jobs = 120, template = list(cpus = 1, partition = "cluster",
  log_file = "/data/sgg2/jenny/projects/PSYMETAB/process_init_%a_clustermq.out"))
