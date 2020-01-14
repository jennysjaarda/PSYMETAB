#####################################################################################
## Author: Jenny Sjaarda
## Project: PSYMETAB_GWAS
##
## Time-stamp: <[pc_merge.r] by JS 2019-10-06 10:42:03 CEST>
##
## Description:
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
# args <- c("/data/sgg2/jenny/projects/PSYMETAB_GWAS/pipeline/extraction/Fred","diabetes", "/data/sgg2/jenny/projects/PSYMETAB_GWAS/pipeline/QC/15_final_processing/final_pca/CEU/pcs.PSYMETAB_GWAS_CEU_unrelated.txt")

workdir <- args[1]
outname <- args[2]
pc_file <- args[3]

extraction_list <- list()
for(chr in 1:22)
{
  if(file.exists(file.path(workdir, outname, paste0(outname,"_chr", chr, "_extract.raw"))))
  {
    extraction_temp <- fread(file.path(workdir, outname, paste0(outname,"_chr", chr, "_extract.raw")), data.table=F)
    extraction_temp <- extraction_temp  %>%
      separate(FID, into = c("ID", "GPCR"), sep = "_") %>%
      dplyr::select(-PAT, -MAT, -SEX, -PHENOTYPE, -IID,-ID)
    chr_car <- as.character(chr)
    extraction_list[[chr_car]] <- extraction_temp # add it to your list
  }
}

pc_data <- fread(pc_file, data.table=F) %>% separate(FID, into = c("ID", "GPCR"), sep = "_") %>% dplyr::select(-IID,-ID)
extraction_list[["pcs"]] <-pc_data
## get list of frequencies for each SNP for each eth with n >100
result <- extraction_list %>%
  reduce(right_join, by = c("GPCR")) # right join ensures that only individuals in PC file are in final file
                                     # these individuals will be unrelated if that is the specified pc file.

write.table(result, paste0(workdir,"/",outname, "/", outname,"_extraction.txt"), row.names=F, col.names=T, quote=F)
