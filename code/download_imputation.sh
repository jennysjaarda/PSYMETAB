#!/bin/bash
# redirect stdout/stderr to a file
exec &> analysis/QC/download_imputation.log

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
curl -sL https://imputationserver.sph.umich.edu/get/1726437/8094e644a19bb19eb29746a7a6bd49e1 | bash

## QC stats
curl -sL https://imputationserver.sph.umich.edu/get/1726441/b6a64d1b662fdeec38a12018e35ae622 | bash

## Logs
curl -sL https://imputationserver.sph.umich.edu/get/1726444/76a2cd8d988008cb4014ada07012a6ce | bash

## Imputation results
curl -sL https://imputationserver.sph.umich.edu/get/1726443/a8df83675382c503e2ebc746f10009b3 | bash


## password provided in email from Michigan imputation server
for chr in {1..22}
do
  unzip -P 'nlCogWNE68vpI' chr_$chr.zip
  gunzip chr$chr.info.gz
done
