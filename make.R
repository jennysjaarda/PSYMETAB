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

# make(qc_prep, log_make = "qc_prep.log", cache_log_file = "cache_log.csv")
#
# make(pre_impute_qc, log_make = "pre_impute_qc.log", cache_log_file = "cache_log.csv",
#   parallelism = "clustermq", jobs = 1, template = list(cpus = 16, partition = "sgg",
#   log_file = "/data/sgg2/jenny/projects/PSYMETAB/pre_impute_qc_%a_clustermq.out"))
#
# ### download files and impute on Michigan server
#
# make(download_impute, log_make = "download_impute.log", cache_log_file = "cache_log.csv",
#   parallelism = "clustermq", jobs = 1, template = list(cpus = 11, partition = "sgg",
#   log_file = "/data/sgg2/jenny/projects/PSYMETAB/download_impute_%a_clustermq.out"))
#
# make(post_impute, log_make = "post_impute.log", cache_log_file = "cache_log.csv",
#   parallelism = "clustermq", jobs = 1, template = list(cpus = 16, partition = "cluster",
#   log_file = "/data/sgg2/jenny/projects/PSYMETAB/post_impute_%a_clustermq.out"))

make(analysis_prep, log_make = "analysis_prep.log", cache_log_file = "cache_log.csv",
  parallelism = "clustermq", jobs = 8, template = list(cpus = 1, partition = "cluster",
  log_file = "/data/sgg2/jenny/projects/PSYMETAB/analysis_prep_%a_clustermq.out"))

make(init_analysis, log_make = "init_analysis.log", cache_log_file = "cache_log.csv",
  parallelism = "clustermq", jobs = 16, template = list(cpus = 8, partition = "cluster",
  log_file = "/data/sgg2/jenny/projects/PSYMETAB/init_analysis_%a_clustermq.out"))

make(process_init, log_make = "process_init.log", cache_log_file = "cache_log.csv",
  parallelism = "clustermq", jobs = 80, template = list(cpus = 2, partition = "cluster",
  log_file = "/data/sgg2/jenny/projects/PSYMETAB/process_init_%a_clustermq.out"))

# make(prs, log_make = "prs.log", cache_log_file = "cache_log.csv",
#   parallelism = "clustermq", jobs = 4, template = list(cpus = 16, partition = "sgg",
#   log_file = "/data/sgg2/jenny/projects/PSYMETAB/prs_%a_clustermq.out"))
