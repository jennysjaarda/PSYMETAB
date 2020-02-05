#####################################################################################
## Author: Jenny Sjaarda
## Project: PSYMETAB_GWAS
##
## Time-stamp: <[relatedness_filter.r] by JS 2019-09-25 10:42:03 CEST>
##
## Description: filteres relatedness tables to get a maximum set of related individuals
##
##
## History:
##
#####################################################################################

args <- commandArgs(trailingOnly = TRUE) #args <- c("chr1-22.kin0",0.0884 )
source("code/settings.R")
kin_file <- as.character(args[1]) # "chr1-22.kin0"
cutoff <- as.numeric(args[2])
out_dir <- as.character(args[3])
out_name <- as.character(args[4])

library(data.table)
library(dplyr)
kin_table <- fread(kin_file, data.table=F)



relatedness_filter <- function (kin_table, cutoff = 0.0884)
{

    data <- kin_table %>% dplyr::filter(KINSHIP > cutoff)
    data <- data[,-which(colnames(data) %in% c("#FID1","FID2", "NSNP"))]
    remove_samples <- vector(mode = "integer")
    connections <- data %>% tidyr::gather(key = "Label", value = "ID",
        -HETHET, -IBS0, -KINSHIP) %>% dplyr::count(ID, sort = TRUE)
    while (!all(connections[["n"]] == 1)) {
        df_id_long <- data %>% tidyr::gather(key = "Label", value = "ID",
            -HETHET, -IBS0, -KINSHIP)
        remove_samples <- df_id_long %>% dplyr::left_join(count(df_id_long,
            ID), by = "ID") %>% dplyr::arrange(desc(n)) %>% utils::head(n = 1) %>%
            .[["ID"]] %>% append(remove_samples, values = .)
        data <- data %>% dplyr::filter(!(ID1 %in% remove_samples |
            ID2 %in% remove_samples))
        connections <- data %>% tidyr::gather(key = "Label",
            value = "ID", -HETHET, -IBS0, -KINSHIP) %>% dplyr::count(ID,
            sort = TRUE)
    }
    remove_samples <- data %>% .[["ID2"]] %>% append(remove_samples,
        values = .)
    return(remove_samples)
}


remove_ids <- relatedness_filter(kin_table,cutoff)

fid_iid <- numeric()
for(i in 1:length(remove_ids))
{
  ID <-remove_ids[i]
  if(length(which(kin_table$ID1==ID))!=0)
  {
    FID <- kin_table[which(kin_table$ID1==ID)[1],"#FID1"]
  }
  if(length(which(kin_table$ID1==ID))==0)
  {
    FID <- kin_table[which(kin_table$ID2==ID)[1],"FID2"]
  }
  fid_iid <-rbind(fid_iid, cbind(FID, ID))
}

write.table(fid_iid, paste0(out_dir, "/", out_name), row.names=F, col.names=F, quote=F)
