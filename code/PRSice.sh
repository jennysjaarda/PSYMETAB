#!/bin/bash
# redirect stdout/stderr to a file
exec &> analysis/GRS/GRS.log

#####################################################################################
## Author: Jenny Sjaarda
## Project: PSYMETAB_GWAS
##
## Time-stamp: <[PRSice.sh] by JS 2019-09-04 12:21:15 CEST>
##
## Description: generate GRS using PRSice
##
## run within project folder: /data/sgg2/jenny/projects/PSYMETAB
##
##
## History:
##
#####################################################################################

project_dir="/data/sgg2/jenny/projects/PSYMETAB"
plink_data="PLINK_091019_0920"
raw_data=$project_dir/data/raw/$plink_data
QC_dir=$project_dir/analysis/QC

input_chip=$project_dir/data/processed/${plink_data}/PSYMETAB
output_name=PSYMETAB_GWAS
sex_file=$project_dir/data/processed/phenotype_data/PSYMETAB_GWAS_sex.txt
eth_file=$project_dir/data/processed/phenotype_data/PSYMETAB_GWAS_eth.txt
duplicate_file=$project_dir/data/processed/phenotype_data/PSYMETAB_GWAS_dupIDs.txt
duplicate_file_set=$project_dir/data/processed/phenotype_data/PSYMETAB_GWAS_dupIDs_set.txt
strand_file=$project_dir/data/processed/GSAMD_UPPC-b37.strand
conversion_file=$project_dir/data/processed/reference_files/rsid_conversion.txt
KG_PRUNED="$x1000G/1000G_chrall.ex_maf_05.pruned"
KG_SAMPLES="$x1000G/integrated_call_samples_v3.20130502.ALL.panel"
ref=$data/HRC.r1-1.GRCh37.wgs.mac5.sites.tab

GRS_dir=$project_dir/analysis/GRS
sample_file=${QC_dir}/15_final_processing/FULL/PSYMETAB_GWAS.FULL_nosex.sample
bgen_file="${QC_dir}/15_final_processing/FULL/PSYMETAB_GWAS.FULL"

base="${data}/consortia/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt"

  Rscript ${PRSice}/PRSice.R --dir . \
      --prsice ${PRSice}/PRSice_linux \
      --base ${base} \
      --target ${bgen_file},${sample_file} \
      --thread 16 \
      --snp SNP \
      --chr CHR \
      --A1 Tested_Allele \
      --A2 Other_Allele \
      --stat BETA \
      --pvalue P \
      --type bgen \
      --out ${GRS_dir}/BMI \
      --bar-levels 0.00000005,0.0000005,0.000005,0.00005,0.0005,0.005,0.05,0.1 \
      --no-regress


      #  --fastscore

      #
      # --keep-ambig \
      # --bar-levels 0.00000005,0.0000005,0.000005,0.00005,0.0005,0.005,0.05,0.1 \
      # --quantile 10 \
      # --fastscore   #### Remove if you want PRSice to test its permutation algorithm
                    #### Use '--all-score' if you want an output at multiple thresholds
  #done



echo PRSice complete
