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

make(pre_impute,parallelism = "future",jobs = 1, console_log_file = "pre_impute_qc.out", template = list(partition = "sgg"))

baseline_gwas_info = define_baseline_inputs(GWAS_input)
interaction_gwas_info = define_interaction_inputs(GWAS_input)
subgroup_gwas_info = define_subgroup_inputs(GWAS_input)
make(pre_impute,parallelism = "future",jobs = 4, console_log_file = "post_impute.out", template = list(partition = "sgg"))

### download files and impute on Michigan server
