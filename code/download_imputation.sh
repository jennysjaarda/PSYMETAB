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
#curl -sL https://imputationserver.sph.umich.edu/get/1726443/a8df83675382c503e2ebc746f10009b3 | bash

wget -nv https://imputationserver.sph.umich.edu/share/results/c43efa9ddd59f323a0ff8a96965101ac/chr_1.zip &
wget -nv https://imputationserver.sph.umich.edu/share/results/701a2663a79134d6bad1ecc054d0c1ae/chr_2.zip &
wget -nv https://imputationserver.sph.umich.edu/share/results/cc49f514a282c3320dd323c7cc3d884f/chr_3.zip &
wget -nv https://imputationserver.sph.umich.edu/share/results/bb603f3028ad7a543ab81c046dc9101d/chr_4.zip &
wget -nv https://imputationserver.sph.umich.edu/share/results/31381560d0e97e2afa6790179bad364/chr_5.zip &
wget -nv https://imputationserver.sph.umich.edu/share/results/3f68159152474a191458728832bf593/chr_6.zip &

wait

wget -nv https://imputationserver.sph.umich.edu/share/results/2e5ca1e4371006d654791a51ebf6d8a6/chr_7.zip &
wget -nv https://imputationserver.sph.umich.edu/share/results/481423d92193fde5517851d4227bec61/chr_8.zip &
wget -nv https://imputationserver.sph.umich.edu/share/results/8c05e2ec128f6e0eca6eb38193b7f612/chr_9.zip &
wget -nv https://imputationserver.sph.umich.edu/share/results/27d241e922cf5ac63f1fa163f4a0b341/chr_10.zip &
wget -nv https://imputationserver.sph.umich.edu/share/results/b99c30a2919f9aa7fe2133233b741b83/chr_11.zip &
wget -nv https://imputationserver.sph.umich.edu/share/results/6e445a9bd9f4a77617034a4c4b517c7b/chr_12.zip &
wget -nv https://imputationserver.sph.umich.edu/share/results/f90c45d66266fc5ac07b5dd28304f8/chr_13.zip &

wait

wget -nv https://imputationserver.sph.umich.edu/share/results/dcd3568a56cf942ae147063494dd5fe/chr_14.zip &
wget -nv https://imputationserver.sph.umich.edu/share/results/ea4352733eedb74a15eada79229b9997/chr_15.zip &
wget -nv https://imputationserver.sph.umich.edu/share/results/c0468e116797376b66a9c56d9c851baf/chr_16.zip &
wget -nv https://imputationserver.sph.umich.edu/share/results/bbe660a5f667e3d53c3efaec9852a065/chr_17.zip &
wget -nv https://imputationserver.sph.umich.edu/share/results/2a66757a414962e039923e46050d291a/chr_18.zip &
wget -nv https://imputationserver.sph.umich.edu/share/results/a5e2c9a79cd0667dbe109040aa972286/chr_19.zip &
wget -nv https://imputationserver.sph.umich.edu/share/results/9685ea3375c0d0b8581eec46df7d3475/chr_20.zip &
wget -nv https://imputationserver.sph.umich.edu/share/results/d2e5f48fe0d4b2157298cb5f8f7daf01/chr_21.zip &
wget -nv https://imputationserver.sph.umich.edu/share/results/9eb608f0dd90c962b2b70b998eae5dd9/chr_22.zip &

wait



## password provided in email from Michigan imputation server
for chr in {1..22}
do
  unzip -P 'nlCogWNE68vpI' chr_$chr.zip
  #gunzip chr$chr.info.gz
done
