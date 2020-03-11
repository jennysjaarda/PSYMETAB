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
curl -sL https://imputationserver.sph.umich.edu/get/1784307/e8087d06a3c3b1bb3cd452fe13ac47c9 | bash

## QC stats
curl -sL https://imputationserver.sph.umich.edu/get/1784311/f3c6bb9ce94d303567dd58f4fb83d7e5 | bash

## Logs
curl -sL https://imputationserver.sph.umich.edu/get/1784314/b962d12b4ae30d700745d1d7781c6566 | bash

## Imputation results
#curl -sL https://imputationserver.sph.umich.edu/get/1784313/e3510c5c86e0345896c55dcac40896da | bash

wget -nv https://imputationserver.sph.umich.edu/share/results/41d1f8ca320f239463baf20418703ef7/chr_1.zip &
wget -nv https://imputationserver.sph.umich.edu/share/results/126870e3dcb6051fef6883894198cea5/chr_2.zip &
wget -nv https://imputationserver.sph.umich.edu/share/results/bb418fed24c123d5ac5309f65ddf1ce8/chr_3.zip &
wget -nv https://imputationserver.sph.umich.edu/share/results/20991e70e221faa4c3356040edbded24/chr_4.zip &

wait

wget -nv https://imputationserver.sph.umich.edu/share/results/832c967c2d500e17601c9944adc38e36/chr_5.zip &
wget -nv https://imputationserver.sph.umich.edu/share/results/a4b962e0d0bae2a5ee6ff9a0f118180/chr_6.zip &
wget -nv https://imputationserver.sph.umich.edu/share/results/2f5b73b5e20f515644c6e4841480c6e4/chr_7.zip &
wget -nv https://imputationserver.sph.umich.edu/share/results/3c50ec69d8bf3cb9bfe60d8f4502f023/chr_8.zip &

wait

wget -nv https://imputationserver.sph.umich.edu/share/results/18b5f7311c693d1ba6de889bc1cef29d/chr_9.zip &
wget -nv https://imputationserver.sph.umich.edu/share/results/5d68073627516ffb992db397d5744c0a/chr_10.zip &
wget -nv https://imputationserver.sph.umich.edu/share/results/a05903356c56dd591835a9d7f3c138c0/chr_11.zip &
wget -nv https://imputationserver.sph.umich.edu/share/results/97b4ad262ba4858722f8227b6bf6264e/chr_12.zip &
wget -nv https://imputationserver.sph.umich.edu/share/results/5717cb7efd61bce24bd091182aea5672/chr_13.zip &
wget -nv https://imputationserver.sph.umich.edu/share/results/35dfa0c123cfe29f9dfd64ecb5b5c5cf/chr_14.zip &

wait


wget -nv https://imputationserver.sph.umich.edu/share/results/48dfefce72cbe7e37d0aaecff20cce42/chr_15.zip &
wget -nv https://imputationserver.sph.umich.edu/share/results/6360d179134edbffb2959fc4077e6ec3/chr_16.zip &
wget -nv https://imputationserver.sph.umich.edu/share/results/4773933541fec29b3f21337a2a2e068c/chr_17.zip &
wget -nv https://imputationserver.sph.umich.edu/share/results/39af636be98c8ffa6b2b7f7ece87947a/chr_18.zip &
wget -nv https://imputationserver.sph.umich.edu/share/results/d4db95a8d17fb7371636680ce46aceee/chr_19.zip &
wget -nv https://imputationserver.sph.umich.edu/share/results/5276b885fc0e418a103b7d20afe1f8fa/chr_20.zip &
wget -nv https://imputationserver.sph.umich.edu/share/results/75627c3cf3be059e146780f0c3f647e6/chr_21.zip &
wget -nv https://imputationserver.sph.umich.edu/share/results/a7170e95153b5dae15c1712cf33c034a/chr_22.zip &

wait

## password provided in email from Michigan imputation server
for chr in {1..22}
do
  unzip -P 'ngA9kfe0RrKSY' chr_$chr.zip
  #gunzip chr$chr.info.gz
done
