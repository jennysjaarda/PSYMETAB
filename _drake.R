setwd(here::here())

source("R/packages.R")
source("R/functions.R")
source("R/settings.R")
source("R/plan.R")

options(clustermq.scheduler = "slurm", clustermq.template = "slurm_clustermq.tmpl")

ukbb_plan <- bind_plans(ukbb_analysis, ukbb_control, ukbb_replication, ukbb_icd_vs_drug)

plan <- bind_plans(qc_prep, pre_impute_qc, download_impute, post_impute, analysis_prep, analysis_prep_case_only, init_analysis, process_init, ukbb_plan)

drake_config(plan, verbose = 1, log_make = "plan.log", cache_log_file = "cache_log.csv",
  parallelism = "clustermq", jobs = 24, template = list(cpus = 4, partition = "cluster2", targets = "ukbb_psy_subgroups_filter_bmi_slope",
  log_file = "/data/sgg2/jenny/projects/PSYMETAB/plan_%a_clustermq.out"), target = "write_celine_aripiprazole_pheno_data"
)
