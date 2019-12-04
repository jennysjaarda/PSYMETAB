#####################################################################################
## Author: Jenny Sjaarda
## Project: PSYMETAB_GWAS
##
## Time-stamp: <[update_pvar.R] by JS 2019-09-25 10:42:03 CEST>
##
## Description: updates pvar file to have rsIDs instead of SNPs in the format of chr:bp:ref:alt
##
##
## History:
##
#####################################################################################
args <- commandArgs(trailingOnly = TRUE)

library(data.table)
library(dtplyr)
library(tidyfast)
library(dplyr)
library(tidyr)

ref_file <- as.character(args[1])
#$data/dbSNP/dbSNP_SNP_list_chr${chr}.txt  #ref_file <- "/data/sgg2/jenny/data/dbSNP/dbSNP_SNP_list_chr1.txt"
pvar_file <- as.character(args[2])
#$processed/${input_basename}.bim #pvar_file <- "/data/sgg2/jenny/projects/PSYMETAB/analysis/QC/08_plink_convert/chr1/chr1.pvar"
out_name <- as.character(args[3])
out_dir <- as.character(args[4])

ref <- fread(ref_file,data.table=F,nThread=16)
pvar <- fread(pvar_file,data.table=F,nThread=16)

colnames(ref) <- c("chr", "bp", "snp", "ref", "alt")
colnames(pvar) <- c("chr", "bp", "ID", "a1" ,"a2", "FILTER","INFO")

pvar$index <- 1:dim(pvar)[1]

#s <- strsplit(as.character(ref$alt), ',')
#ref_update <- data.frame(chr=rep(ref$chr, sapply(s, FUN=length)),alt=unlist(s) )


ref_update <- ref %>% lazy_dt() %>%
  mutate(alt = strsplit(as.character(alt), ",")) %>%
  dt_hoist(col=alt) %>% as_tibble()

# older slower code
# ref_update <- temp %>%
#   mutate(alt = strsplit(as.character(alt), ",")) %>%
#   unnest(alt)

# ref_copy <- ref
ref <- ref_update

simple_joint <- full_join(pvar, ref) # matches by chr, bp

simple_joint$new_id <- NA
simple_joint$ref_match <- NA
simple_joint$duplicate <- NA ## duplicate based on rsid

simple_joint <- simple_joint %>% filter(!is.na(index))

select_rows1 <- which(simple_joint$a1==simple_joint$ref &
  simple_joint$a2==simple_joint$alt)
simple_joint[select_rows1, "new_id"] <- simple_joint$snp[select_rows1]
simple_joint[select_rows1, "ref_match"] <- 1

select_rows2 <- which(simple_joint$a1==simple_joint$alt &
  simple_joint$a2==simple_joint$ref)
simple_joint[select_rows2, "new_id"] <- simple_joint$snp[select_rows2]
simple_joint[select_rows2, "ref_match"] <- 1

select_rows3 <- which(!(simple_joint$a1==simple_joint$alt & simple_joint$a2==simple_joint$ref) &
  !(simple_joint$a1==simple_joint$ref & simple_joint$a2==simple_joint$alt))
simple_joint[select_rows3, "new_id"] <- paste(simple_joint$chr[select_rows3],simple_joint$bp[select_rows3],simple_joint$a1[select_rows3],simple_joint$a2[select_rows3],sep=":")
simple_joint[select_rows3, "ref_match"] <- 0

select_rows4 <- which(is.na(simple_joint$snp))
simple_joint[select_rows4, "new_id"] <- paste(simple_joint$chr[select_rows4],simple_joint$bp[select_rows4],simple_joint$a1[select_rows4],simple_joint$a2[select_rows4],sep=":")
simple_joint[select_rows4, "ref_match"] <- 0

simple_joint <- simple_joint[order(simple_joint$index, -simple_joint$ref_match),]

duplicate_ids <- which(duplicated(simple_joint$index))
simple_joint[duplicate_ids, "duplicate"] <- 1
simple_joint[-duplicate_ids, "duplicate"] <- 0

output <- subset(simple_joint, simple_joint$duplicate==0)
## make sure all SNP names are unique
still_dups <- which(duplicated(output$new_id))
output[still_dups,"new_id"] <- paste(output$new_id,output$a1, output$a2, sep=":")[still_dups]


output <- output[,c("ID", "new_id")]


if(dim(output)[1]==dim(pvar)[1])
{
  cat("The map file length matches the original bim file!\n")

}

if(dim(output)[1]!=dim(pvar)[1])
{
  cat("The map file length does not match the original pvar file, go back and check this pvar file: \n")
  cat(pvar_file)

}

data.table::fwrite(output, paste0(out_dir,"/", out_name), row.names=F, col.names=F, sep="\t", quote=F,nThread=16)
