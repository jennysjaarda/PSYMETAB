output_name=PSYMETAB_GWAS
#rs7903146

prune_threshold=0.001
clump_kb=10000 #Clumping window
clump_p1=1 #Clumping sig level for index SNPs
clump_p2=1 #Clumping sig level for secondary SNPs

  plink2 --pfile FULL/${output_name}.FULL_chr10 \
    --chr 10 \
    --make-pfile \
    --from-kb 114750 --to-kb 114760 \
    --out PRSice_test \
    --threads 16

  plink2 --pfile PRSice_test \
    --export bgen-1.2 id-paste=iid \
    --out PRSice_test \
    --threads 16

  plink2 --pfile PRSice_test \
    --make-bfile \
    --out PRSice_test \
    --threads 16

  plink2 --pfile PRSice_test \
    --recode A \
    --out PRSice_test \
    --threads 16

  plink2 --pfile PRSice_test \
    --export bgen-1.2 bits=8 id-paste=iid \
    --out PRSice_test \
    --threads 16


bgenix -g file.bgen -index

extract_bgen <- function(snps, file){

  num_snps <- length(snps)
  bgen_data <- bgen.load(file=file, rsids = snps, max_entries=4)
  snp_map <- bgen_data$variants
  snp_map$AF <- NA
  geno <- numeric()

  for(i in 1:dim(snp_map)[1])
  {
      snp <- as.character(snp_map[["rsid"]][i])
      snp_temp <- bgen_data$data[snp,,]
      snp_temp <- as.data.frame(snp_temp)
      snp_temp[,4] <- snp_temp[,2]+2*snp_temp[,3]
      snp_map$AF[i] <- mean(snp_temp[,4], rm.na=T)/2
      snp_out <- snp_temp[,4]
      geno <- cbind(geno, snp_out)
  }
  colnames(geno) <- bgen_data$variants$rsid
  return(list(geno_data=geno,snp_map=snp_map))
}
t <- extract_bgen("rs117796541", "PRSice_test.bgen")


--target /data/sgg2/jenny/projects/PSYMETAB/analysis/QC/15_final_processing/PRSice_test,/data/sgg2/jenny/projects/PSYMETAB/analysis/QC/15_final_processing/FULL/PSYMETAB_GWAS.FULL_nosex.sample \

Rscript ${PRSice}/PRSice.R --dir . \
    --prsice ${PRSice}//PRSice_linux \


  Rscript PRSice.R --dir . \
    --prsice ./PRSice_linux \

2076956
gawk '$2==10 && $3 > 114750 && $3 > 114760' /data/sgg2/jenny/data/consortia/formatted/Locke_UKBB_meta_BMI.txt > BMI_sub.txt

gawk '$2==10' /data/sgg2/jenny/data/consortia/formatted/Locke_UKBB_meta_BMI.txt > BMI_sub.txt
awk 'NR==1' /data/sgg2/jenny/data/consortia/formatted/Locke_UKBB_meta_BMI.txt > header.txt
#grep -w rs7903146 /data/sgg2/jenny/data/consortia/formatted/Locke_UKBB_meta_BMI.txt > BMI_sub.txt

cat header.txt BMI_sub.txt > BMI_sub_header.txt
Rscript ${PRSice}/PRSice.R --dir . \
    --prsice ${PRSice}//PRSice_linux \
    --target /data/sgg2/jenny/projects/PSYMETAB/analysis/QC/15_final_processing/FULL/PSYMETAB_GWAS.FULL,/data/sgg2/jenny/projects/PSYMETAB/analysis/QC/15_final_processing/FULL/PSYMETAB_GWAS.FULL_nosex.sample \
    --stat BETA \
    --type bgen \
    --binary-target T \
    --all-score \
    --no-regress \
    --base /data/sgg2/jenny/data/consortia/formatted/Mahajan_T2Dbmiadj.txt \
    --snp SNP \
    --stat BETA \
    --A1 EFFECT_ALLELE \
    --A2 OTHER_ALLELE \
    --chr CHR \
    --pvalue PVAL \
    --thread 16 \
    --out T2D_test
