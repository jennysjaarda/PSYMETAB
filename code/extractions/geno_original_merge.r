#####################################################################################
## Author: Jenny Sjaarda
## Project: PSYMETAB_GWAS
##
## Time-stamp: <[pc_merge.r] by JS 2020-07-21 10:42:03 CEST>
##
## Description: this script merges the data that was extracted from the imputation file with the data extracted from the original PLINK file
##  Some SNPs were genotyped on the original chip but were not retained because they were I/Ds.
##
## Inputs:
##
## History:
##
#####################################################################################

########################################
#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
library(tidyverse)
library(data.table)
library(tidyr)
# args <- c("/data/sgg2/jenny/projects/PSYMETAB/data/processed/extractions/Nicolas","SNP_list_23062020",
 "/data/sgg2/jenny/projects/PSYMETAB/data/raw/extractions/Nicolas/SNP_list_23062020.txt")


workdir <- args[1]
outname <- args[2]
input_file <- args[3]

original_extraction <- fread(paste0(workdir,"/",outname, "/", outname,"_extraction.txt"), header=T, data.table=F)

snps_requested <- fread(input_file, header=F, data.table=F)
snps_missing <- read.table(paste0(workdir,"/",outname, "/", outname,"_missing_snps.txt"), header=F)

missing_extraction <- fread(file.path(workdir, outname, paste0(outname, "_missing_snps_extract.raw")), data.table=F)

snps_extracted <- colnames(missing_extraction)[-c(1:6)]
snps_extracted <- sub("_[^_]+$", "", snps_extracted)

full_extraction <- missing_extraction  %>%
  rename(GPCR = IID) %>%
  dplyr::select(-PAT, -MAT, -SEX, -PHENOTYPE, -FID) %>%
  right_join(original_extraction, by = "GPCR")


snps_missing2 <- snps_missing[which(!snps_missing[,1] %in% snps_extracted),]

write.table(full_extraction, paste0(workdir,"/",outname, "/", outname,"_extraction_geno.txt"), row.names=F, col.names=T, quote=F)
write.table(snps_missing2, paste0(workdir,"/",outname, "/", outname,"_missing_geno_snps.txt"), row.names=F, col.names=F, quote=F)
