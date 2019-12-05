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

source("code/settings.r")
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
colnames(snp_weights) <- c("FID", "IID", "num_SNPs",  paste0("snpweights_PC",1:3), "YRI_perc", "CEU_perc", "EA_perc", "NA_perc")

pheno_eth <- read.table(as.character(args[5]),header=F) # reported ethnicity
colnames(pheno_eth) <- c(  "FID", "IID", "ETH")
pheno_eth <- pheno_eth %>%
  unite("ID2", FID:IID, remove = FALSE, sep="_")
  #mutate_at("ID2", as.factor)

data <- merge(snp_weights, pc_data, by=c("FID", "IID"))

factor_eths <- function(x) {ifelse(x< 0.8, FALSE, TRUE)}
factor_yri <- function(x) {ifelse(x< 0.7, FALSE, TRUE)}

data_clean <- data %>%
    mutate_at(c("CEU_perc","EA_perc","NA_perc"), list(factor=factor_eths)) %>%
    mutate_at(c("YRI_perc"), list(YRI_perc_factor=factor_yri)) %>%

    rename_at(vars(ends_with("_factor")),
      list( ~ str_replace(., "_perc", ""))) %>%
    dplyr::select(ends_with("_factor"), FID)  %>%
    gather(eth, value, -FID) %>%
    filter(value==TRUE) %>%
    dplyr::select(-value) %>%
    left_join(data %>% dplyr::select(-ends_with("_factor")), .) %>%
    arrange(IID)  %>%
    mutate_at(c("YRI_perc","CEU_perc","EA_perc","NA_perc"), list(factor=factor_eths)) %>%
    rename_at(vars(ends_with("_factor")),
      list( ~ str_replace(., "_perc", "")))  %>%
    replace_na(list(eth = "MIXED"))  %>%
    mutate(genetic_eth = recode(eth, CEU_factor = "CEU",
         EA_factor = "EA",
         NA_factor = "NA",
         YRI_factor = "YRI"
     )) %>%
     mutate_at("IID", as.character) %>%
     left_join(., pheno_eth %>% dplyr::select(ID2,ETH), by = c("IID" = "ID2") ) %>%
     dplyr::rename(reported_eth = ETH)

pca_plot <- function(col, colname, title)
{
  ggplot(data_clean) +
  # add scatter points
  geom_point(aes_string(x = "PC1", y = "PC2", col = col),
             alpha = 0.5) +

  # change labels
  labs(title = title,
       x = "PC1",
       y = "PC2",
       #size = "Percent European",
       col = colname) +
  # change the size scale
  scale_size(range = c(0.1, 10)) +
  # add a nicer theme
  theme_classic() +
  # place legend at top and grey axis lines
  theme(legend.position = "top")
}


### plot according to genetic eth
for(col_type in c("genetic", "reported"))
{


  pdf(paste0(output_dir, "/", study_name,"_PCA_",col_type, "_eth.pdf"))
  out_plot <- pca_plot(paste0(col_type,"_eth"), paste0(str_to_title(col_type), " Ancestry"),
  paste0(study_name, " PCA according to ",col_type, " ancestry"))
  print(out_plot)
  dev.off()

}

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
