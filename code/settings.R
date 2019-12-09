#####################################################################################
## Author: Jenny Sjaarda
## Project: PSYMETAB_GWAS
##
## Time-stamp: <[settings.r] by JS 2019-07-17 14:27:03 CEST>
##
## Description:
##
##
## History:
##
#####################################################################################

### Set data and project directories

project_dir <- "/data/sgg2/jenny/projects/PSYMETAB/"
SGG_generic <- "/data/sgg2/jenny/SGG_generic/"

source(paste0(SGG_generic,"/scripts/settings.r"))
source("code/packages.R")

options(clustermq.scheduler = "slurm", clustermq.template = "slurm_clustermq.tmpl")
drake_hpc_template_file("slurm_clustermq.tmpl")

### data sources
plink_ped_raw <- "data/raw/PLINK_091019_0920/PSYMETAB_GWAS"
plink_bed_out <- "data/processed/PLINK_091019_0920/PSYMETAB_GWAS"
rsconv_raw_file <- "data/raw/reference_files/GSAMD-24v2-0_20024620_A1_b151_rsids.txt"
qc_pheno_file <- "data/raw/phenotype_data/QC_sex_eth.xlsx"
pheno_file <- "data/raw/phenotype_data/PHENO_GWAS_241019_noaccent.csv"
pc_dir <- "analysis/QC/15_final_processing/final_pca"

### location of codes
pre_imputation_script <- "code/pre_imputation_qc.sh"
post_imputation_script <- "code/post_imputation_qc.sh"
download_imputation_script <- "code/download_imputation.sh"
check_imputation_script <- "code/check_imputation.sh"
final_processing_script <-  "code/final_processing.sh"

### define variables

eths <- c("CEU", "EA", "MIXED", "NA", "YRI")
study_name <- "PSYMETAB_GWAS"
leeway_time <- 90 ## plus or minus on follow-up to double check that the "Mois" column make sense.
# For example, if one visit is marked month 3 on September 5/15; and the previous visit marked month 0 on June 15/19, the script below will check if:
# the difference between the two dates (September 5/15 - June 15/19) is within the difference between the two months (3 = 90 days) +/- the 'leeway_time'

min_follow_up <- 14

# extractions
#Anais_output <-
#Fred_output <-
#Aurelie_output <-

# AP1 pour faible risque de prise de poids :
# Amisulpride, Aripiprazole, Brexpiprazole, Cariprazine, Carbamazepine, Chlorprothixene, Flupentixol, Fluphenazine, Haloperidol, Lurasidone,
# Pipampérone, Promazine, Sertindole, Sulpiride, Tiapride

# AP2 pour risque intermédiaire de prise de poids :
# Amitriptyline, Asenapine, Clomipramine, Dibenzepine, Doxepine, Imipramine, Lévomépromazine, Lithium, Mirtazapine,
# Nortriptyline, Opipramol, Palipéridone, Quétiapine, Rispéridone, Trimipramine,  Zuclopenthixol

# AP3 pour risque élevé de prise de poids :
# Clozapine, Olanzapine, Valproate

high_inducers <- c("Olanzapine", "Clozapine" ,"Valproate")

med_inducers <- c("Amitriptyline", "Asenapine", "Clomipramine" ,  "Doxepine", "Levomepromazine" ,
 "Lithium", "Mirtazapine", "Paliperidone", "Risperidone", "Quetiapine", "Trimipramine", "Zuclopenthixol")

low_inducers <- c("Amisulpride", "Aripiprazole", "Brexpiprazole", "Cariprazine", "Carbamazepine", "Chlorprothixene",
"Flupentixol", "Fluphenazine","Haloperidol","Lurasidone", "Pipamperone", "Sertindole", "Sulpiride", "Tiapride")


test_drugs <- tibble(class=c("all_high_inducers", "olanzapine_clozapine", "valproate"), drugs=list(high_inducers, c("Olanzapine", "Clozapine"), c( "Valproate")))
