
#####################################################################################
## Author: Jenny Sjaarda
## Project: PSYMETAB_GWAS
##
## Time-stamp: <[maf_check.r] by JS 2019-10-06 10:42:03 CEST>
##
## Description: creates frequency file (one column for ethnicity with n > threshold) and
## list of SNPs to remove that don't pass MAF threshold in any ethnic group
##
## Inputs: study name, afreq file (with 'ETHNICITY_NAME'as placeholder for ethnic group),
## psam file (to get sample length count), output directory, MAF threshold, n threshold
##
## History:
##
#####################################################################################

args <- commandArgs(trailingOnly = TRUE)
source("code/settings.r")

# args <- c("PSYMETAB_GWAS",
#  "/data/sgg2/jenny/projects/PSYMETAB_GWAS/pipeline/QC/14_mafcheck/ETHNICITY_NAME/PSYMETAB_GWAS.ETHNICITY_NAME.afreq",
#  "/data/sgg2/jenny/projects/PSYMETAB_GWAS/pipeline/QC/14_mafcheck/ETHNICITY_NAME/PSYMETAB_GWAS.ETHNICITY_NAME.psam",
#  "/data/sgg2/jenny/projects/PSYMETAB_GWAS/pipeline/QC/14_mafcheck",
#  0.01,100
# )

output_dir <- args[4]
maf_threshold <- as.numeric(args[5])
n_threshold <- as.numeric(args[6])

out <- list()
eths <- for (eth in c("CEU","EA","MIXED", "NA", "YRI"))
{
  freq_file <- gsub("ETHNICITY_NAME", eth, args[2])
  sample_file <- gsub("ETHNICITY_NAME", eth, args[3])
  samples <- fread(sample_file, data.table=F, header=T)
  n_samples <- dim(samples)[1]
  if (n_samples >=100)
  {
  freqs <- fread(freq_file, data.table=F, header=T)
  new_col <- paste0("FREQ_", eth)
  snp_data <- freqs %>% mutate_at("ALT_FREQS", as.numeric) %>%
    dplyr::select(ID, ALT_FREQS) %>%
    rename(!!new_col := ALT_FREQS)
   out[[eth]] <- snp_data
 }

}

## get list of frequencies for each SNP for each eth with n >100
mafs <- out %>%
  reduce(left_join, by = "ID")
## filter to only SNPs that have at least one MAF > threshold
filter_maf <- mafs %>%
  filter_at(vars(starts_with("FREQ_")), all_vars(.< maf_threshold | .> (1-maf_threshold))) %>%
  dplyr::select(ID)

write.table(mafs, paste0(output_dir, "/", args[1], "_maf.txt"), row.names=F, col.names=T, quote=F)
write.table(filter_maf, paste0(output_dir, "/", args[1], "_low_maf_snps.txt"), row.names=F, col.names=F, quote=F)
