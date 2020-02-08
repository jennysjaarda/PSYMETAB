#####################################################################################
## Author: Jenny Sjaarda
## Project: PSYMETAB_GWAS
##
## Time-stamp: <[ethnicity_check.r] by JS 2019-10-06 10:42:03 CEST>
##
## Description: creates PCA plot based on inferred sex from snpWeights and reported sex from pheno files
##
## Inputs: study name, output directory, PC file (from pca in plink), snpWeights file using 'NA' file, reported ethnicity file
##
## History:
##
#####################################################################################

########################################
#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
#rmote::start_rmote()

source("code/packages.R")
source("code/functions.R")
source("code/settings.R")

#
# args <- c("PSYMETAB_GWAS",
#  "/data/sgg2/jenny/projects/PSYMETAB_GWAS/pipeline/QC/12_ethnicity_admixture/",
#  "/data/sgg2/jenny/projects/PSYMETAB_GWAS/pipeline/QC/12_ethnicity_admixture/pca/PSYMETAB_GWAS_projections.txt",
#  "/data/sgg2/jenny/projects/PSYMETAB_GWAS/pipeline/QC/12_ethnicity_admixture/snpweights/PSYMETAB_GWAS.NA.predpc",
#  "/data/sgg2/jenny/projects/PSYMETAB_GWAS/data/processed/phenotype_data/PSYMETAB_GWAS_eth.txt"
#
# )


study_name<-as.character(args[1]) # OUTPUT_NAME
output_dir <- as.character(args[2])
pc_data<-read.table(as.character(args[3]),header=T) # pcs from plink --pca head
snp_weights<-read.table(as.character(args[4]),header=F) # snpweights result using snpwt.NA file
eth_file <- as.character(args[5])

data_clean <- munge_snpweights(study_name, pc_data, snp_weights, eth_file)

for(eth_var in unique(data_clean$genetic_eth))
{
write_data <- data_clean %>%
  filter(genetic_eth == eth_var) %>%
  dplyr::select(IID, FID)
write.table(write_data, paste0(output_dir, "/", study_name, "_", str_to_upper(eth_var), "_samples.txt"), row.names=F, col.names=F, quote=F)
}

#
# table(data_clean$genetic_eth, data_clean$reported_eth)
#
#       AFRICAN EAST_ASIAN EUROPEAN LATIN OTHER SOUTH_ASIAN UNKNOWN
# CEU         0          0     1556     0    72           1     541
# EA          0         17        0     0     2           3       7
# mixed      54          4       77     0   100          55     138
# NA          0          0        0     1     3           0       1
# YRI        47          4        1     1     3           0      27
