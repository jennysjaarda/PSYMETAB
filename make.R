# make.R

# Master script for running PSYMETAB analyses.
# This project uses drake to manage workflows.
# For more information about drake, see
# https://ropensci.github.io/drake/

# Setup ----

# Set working directory
setwd(here::here())

# Load packages ----
source("code/packages.R")

# Load functions and plans  ----
source("code/functions.R")
source("code/settings.R")
source("code/plan.R")

# Create output folders  ---
dir.create(plink_bed_out,showWarnings = F)
dir.create("analysis/QC", showWarnings = F)
dir.create("data/processed/phenotype_data/GWAS_input", showWarnings = F)
create_analysis_dirs("analysis/GWAS")

# Run analyses ----

## Run plans with different requirements seperately because the `future` backend is so slow

make(qc_prep, console_log_file = "qc_prep.log", cache_log_file = "cache_log.csv")

make(pre_impute_qc, console_log_file = "pre_impute_qc.log", cache_log_file = "cache_log.csv",
  parallelism = "clustermq", jobs = 1, template = list(cpus = 16, partition = "sgg",
  log_file = "/data/sgg2/jenny/projects/PSYMETAB/pre_impute_qc_%a_clustermq.out"))

### download files and impute on Michigan server

make(download_impute, console_log_file = "download_impute.log", cache_log_file = "cache_log.csv",
  parallelism = "clustermq", jobs = 1, template = list(cpus = 11, partition = "sgg",
  log_file = "/data/sgg2/jenny/projects/PSYMETAB/download_impute_%a_clustermq.out"))

make(post_impute, console_log_file = "post_impute.log", cache_log_file = "cache_log.csv",
  parallelism = "clustermq", jobs = 1, template = list(cpus = 16, partition = "sgg",
  log_file = "/data/sgg2/jenny/projects/PSYMETAB/post_impute_%a_clustermq.out"))

make(analysis_prep, console_log_file = "analysis_prep.log", cache_log_file = "cache_log.csv",
  parallelism = "clustermq", jobs = 4, template = list(cpus = 1, partition = "sgg",
  log_file = "/data/sgg2/jenny/projects/PSYMETAB/analysis_prep_%a_clustermq.out"))

make(init_analysis, console_log_file = "init_analysis.log", cache_log_file = "cache_log.csv",
  parallelism = "clustermq", jobs = 3, template = list(cpus = 16, partition = "cluster",
  log_file = "/data/sgg2/jenny/projects/PSYMETAB/init_analysis_%a_clustermq.out"))

make(prs, console_log_file = "prs.log", cache_log_file = "cache_log.csv",
  parallelism = "clustermq", jobs = 4, template = list(cpus = 16, partition = "sgg",
  log_file = "/data/sgg2/jenny/projects/PSYMETAB/prs_%a_clustermq.out"))


# Now definied within dynamic target
# baseline_gwas_info = define_baseline_inputs(GWAS_input)
# interaction_gwas_info = define_interaction_inputs(GWAS_input)
# subgroup_gwas_info = define_subgroup_inputs(GWAS_input)
