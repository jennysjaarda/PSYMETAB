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

project_dir <- "/data/sgg2/jenny/projects/PSYMETAB_GWAS/"
SGG_generic <- "/data/sgg2/jenny/SGG_generic/"
original_plink_data <- "PLINK_091019_0920"
source(paste0(SGG_generic,"/scripts/settings.r"))
source(paste0(project_dir,"/scripts/functions.r"))

### define variables


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
