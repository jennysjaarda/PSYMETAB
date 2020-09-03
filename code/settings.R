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
consortia_dir <- "/data/sgg2/jenny/data/consortia"
source(paste0(SGG_generic,"/scripts/settings.r"))

## register clustermq and future plans
options(clustermq.scheduler = "slurm", clustermq.template = "slurm_clustermq.tmpl")
future::plan(batchtools_slurm, template = "slurm_batchtools.tmpl")

### data sources
plink_ped_raw <- "data/raw/PLINK_091019_0920/PSYMETAB_GWAS"
plink_bed_out <- str_replace(plink_ped_raw,"raw","processed") #same as input_chip

rsconv_raw_file <- "data/raw/reference_files/GSAMD-24v2-0_20024620_A1_b151_rsids.txt"
qc_pheno_file <- "data/raw/phenotype_data/QC_sex_eth.xlsx"
pheno_file <- "data/raw/phenotype_data/PHENO_GWAS_160420_corr_noaccent.csv"
pc_dir <- "analysis/QC/15_final_processing/final_pca"
caffeine_file <- "data/raw/phenotype_data/CAF_Sleep_Jenny_09_06_2020.xlsx"
#plink_output_dir <- "/data/sgg3/jenny/projects/PSYMETAB/analysis/GWAS"
plink_output_dir <- "analysis/GWAS"

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

anonymization_error <- 7 ## the anonymization procedure adds a random number of days to each visit date and it makes it very difficult to join two datasets
# instead of a exact merge (default R merge function), use `fuzzy_join` which will join more loosely. See example of it's use in `merge_pheno_caffeine` function in `code\functions.R`

follow_up_limit <- 180 ## restrict GWAS to number of days follow-up

follow_up_6mo <- 180
follow_up_3mo <- 90
follow_up_1mo <- 30


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

ukbb_files <- tibble(
  name = c("coffee_consumed"),
  file = c("IVs/clump/100240.gwas.imputed_v3.both_sexes.IVs/100240.gwas.imputed_v3.both_sexes.IVs_unpruned.txt"))

ukbb_prs_analysis <- tibble(
  name = c("coffee_consumed"),
  linear_pheno_columns = list(c("logCAF", "logPARAX", "logTHEOPH", "logTHEOBR", "logCAFPARAX", "logCAFPARAXTHEOPH")),
  glm_pheno_columns = list(c("Sleep_disorder")),
  covars = list(c("Age_caffeine", "Age_caffeine_sq")),
)

high_inducers <- c("Olanzapine", "Clozapine" ,"Valproate")

med_inducers <- c("Amitriptyline", "Asenapine", "Clomipramine" ,  "Doxepine", "Levomepromazine" ,
 "Lithium", "Mirtazapine", "Paliperidone", "Risperidone", "Quetiapine", "Trimipramine", "Zuclopenthixol")

low_inducers <- c("Amisulpride", "Aripiprazole", "Brexpiprazole", "Cariprazine", "Carbamazepine", "Chlorprothixene",
"Flupentixol", "Fluphenazine","Haloperidol","Lurasidone", "Pipamperone", "Sertindole", "Sulpiride", "Tiapride")


drug_classes <- c("all", "olanz_cloz", "valproate", "olanz", "cloza", "risp", "quet", "ocq", "ami")
test_drugs <- tibble(class=drug_classes, drugs=list(high_inducers, c("Olanzapine", "Clozapine"), c("Valproate"), c("Olanzapine"),
  c("Clozapine"), c("Risperidone"), c("Quetiapine"), c("Olanzapine", "Clozapine", "Quetiapine"), c("Amisulpride")))

baseline_vars <- c("BMI", "Creatinine", "Glucose", "CholesterolHDL", "LDL", "Tryglycerides")
interaction_vars <- c("BMI", "LDL")
caffeine_vars <-  c("logCAF", "logTHEOBR", "logPARAX", "Sleep_disorder")

standard_covars <- c(paste0("PC", 1:20), "sex")
baseline_covars <- c("Age_sq_Drug_1","Age_Drug_1")
caffeine_covars <- c("Age_caffeine", "Age_caffeine_sq")

outcomes <- c("_slope", "_slope_6mo", "_slope_weight", "_slope_weight_6mo",
  "_change", "_change_1mo", "_change_3mo", "_change_6mo")


interaction_outcome <- apply(expand.grid(interaction_vars,
                                        c("_slope", "_slope_6mo", "_slope_weight", "_slope_weight_6mo",
                                          "_change", "_change_1mo", "_change_3mo", "_change_6mo")),
                             1, paste, collapse="")

interaction_outcome_combinations <- expand.grid(interaction_outcome, dplyr::pull(test_drugs, class))
GWAS_models <- tibble(outcome_variable=c(as.character(interaction_outcome_combinations[,1]),baseline_vars),
                      interaction_variable=c(as.character(interaction_outcome_combinations[,2]),rep(NA, length(baseline_vars))),
                      model=c(rep("interaction", dim(test_drugs)[1]*length(interaction_outcome)), rep("linear", length(baseline_vars))))


gw_sig <- 5e-08
gw_sig_nominal <- 5e-06

maf_threshold <- 0.05
info_threshold <- 0.8


## UKB follow-up relevant variables



UKBB_dir <- "/data/sgg3/data/UKBB/"

#eid_var<- c("53", "21001")
date_followup <- "53"
bmi_var<- "21001"
sex_var <- "31"
age_var <- "21022"
