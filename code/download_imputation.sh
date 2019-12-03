#!/bin/bash
# redirect stdout/stderr to a file
exec &> download_imputation.log

### DEFINE INPUTS

project_dir="/data/sgg2/jenny/projects/PSYMETAB"
plink_data="PLINK_091019_0920"
raw_data=$project_dir/data/raw/$plink_data

output_dir=$project_dir/analysis/QC
input_chip=$project_dir/data/processed/${plink_data}/PSYMETAB
output_name=PSYMETAB_GWAS
sex_file=$project_dir/data/processed/phenotype_data/PSYMETAB_GWAS_sex.txt
eth_file=$project_dir/data/processed/phenotype_data/PSYMETAB_GWAS_eth.txt
duplicate_file=$project_dir/data/processed/phenotype_data/PSYMETAB_GWAS_dupIDs.txt
strand_file=$project_dir/data/processed/GSAMD_UPPC-b37.strand
conversion_file=$project_dir/data/processed/reference_files/rsid_conversion.txt
KG_PRUNED="$x1000G/1000G_chrall.ex_maf_05.pruned"
KG_SAMPLES="$x1000G/integrated_call_samples_v3.20130502.ALL.panel"
ref=$data/HRC.r1-1.GRCh37.wgs.mac5.sites.tab
flashpca=$bin/flashpca

#############################################################################################################
#---------------------------------- DOWNLOAD IMPUTATION -------------------------------------------------#
#############################################################################################################


if [ ! -d "${output_dir}/06_imputation_get" ] ; then
  mkdir ${output_dir}/06_imputation_get
fi
cd ${output_dir}/06_imputation_get


## QC report
curl -sL https://imputationserver.sph.umich.edu/get/1640984/764cffab51cb0833521050cef2d98bc | bash

## QC stats
curl -sL https://imputationserver.sph.umich.edu/get/1640988/1fcd9a0e4f8b6c4c52510d0146a720fa | bash

## Logs
curl -sL https://imputationserver.sph.umich.edu/get/1640991/d20e8108c6966e7cc2c76f92c368b2a2 | bash

## Imputation results
curl -sL https://imputationserver.sph.umich.edu/get/1640990/7aaf768775a0b8c7f58afbf109627ac | bash


## password provided in email from Michigan imputation server
for chr in {1..22}
do
  unzip -P 'v35orVrSjBKI6W' chr_$chr.zip
  gunzip chr$chr.info.gz
done
