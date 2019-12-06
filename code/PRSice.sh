#!/bin/bash

#####################################################################################
## Author: Jenny Sjaarda
## Project: PSYMETAB_GWAS
##
## Time-stamp: <[PRSice.sh] by JS 2019-09-04 12:21:15 CEST>
##
## Description: generate GRS using PRSice
##
## run within project folder: /data/sgg2/jenny/projects/PSYMETAB_GWAS
##
##
## History:
##
#####################################################################################



project_dir="/data/sgg2/jenny/projects/PSYMETAB_GWAS"
QC_dir=$project_dir/pipeline/QC

GRS_dir=$project_dir/pipeline/GRS
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
