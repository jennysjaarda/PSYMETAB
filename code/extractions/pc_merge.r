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
# args <- c("/data/sgg2/jenny/projects/PSYMETAB/data/processed/extractions/Nicolas","SNP_list_23062020",
# "/data/sgg2/jenny/projects/PSYMETAB/analysis/QC/15_final_processing/final_pca/ETH/pcs.PSYMETAB_GWAS_ETH_unrelated.txt",
#  "/data/sgg2/jenny/projects/PSYMETAB/data/raw/extractions/Nicolas/SNP_list_23062020.txt")

cat("Extracting PC data..\n")

workdir <- args[1]
outname <- args[2]
pc_file <- args[3]
input_file <- args[4]
snps_requested <- fread(input_file, header=F, data.table=F)
snps_extracted <- numeric()
extraction_list <- list()
for(chr in 1:22)
{
  if(file.exists(file.path(workdir, outname, paste0(outname,"_chr", chr, "_extract.raw"))))
  {
    extraction_temp <- fread(file.path(workdir, outname, paste0(outname, "_chr", chr, "_extract.raw")), data.table=F)
    extraction_temp <- extraction_temp  %>%
      separate(FID, into = c("ID", "GPCR"), sep = "_") %>%
      dplyr::select(-PAT, -MAT, -SEX, -PHENOTYPE, -IID,-ID)

    snps_temp <- colnames(extraction_temp)[-1]
    snps_temp <- sub("_[^_]+$", "",snps_temp)
    snps_extracted <- c(snps_extracted, snps_temp)
    chr_car <- as.character(chr)
    extraction_list[[chr_car]] <- extraction_temp # add it to your list
  }
}

eth_result <- list()
for (eth in c("CEU", "EA", "MIXED", "NA", "YRI")){


  eth_pc_file <- str_replace_all(pc_file, "ETH", eth)
  if(file.exists(eth_pc_file)){
    pc_data <- fread(eth_pc_file, data.table=F) %>% separate(FID, into = c("ID", "GPCR"), sep = "_") %>% dplyr::select(-IID,-ID)
    extraction_list[["eth"]] <- as.data.frame(cbind("GPCR"=as.character(pc_data[["GPCR"]]), eth)) %>% mutate(GPCR = as.character(GPCR))
    extraction_list[["pcs"]] <-pc_data

    ## get list of frequencies for each SNP for each eth with n >100
    result <- extraction_list %>%
      reduce(right_join, by = c("GPCR")) # right join ensures that only individuals in PC file are in final file
                                         # these individuals will be unrelated if that is the specified pc file.
    eth_result[[eth]] <- result # add it to your list

  }

}

eth_merged_result <- rbindlist(eth_result, fill = TRUE)

snps_missing <- snps_requested[which(!snps_requested[,1] %in% snps_extracted),]

write.table(eth_merged_result, paste0(workdir,"/",outname, "/", outname,"_extraction.txt"), row.names=F, col.names=T, quote=F)
write.table(snps_missing, paste0(workdir,"/",outname, "/", outname,"_missing_snps.txt"), row.names=F, col.names=F, quote=F)

## check the missing SNPs in the original PLINK files
