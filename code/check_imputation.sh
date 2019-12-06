#!/bin/bash
# redirect stdout/stderr to a file
exec &> analysis/QC/check_imputation.log
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
#---------------------------------- CHECK IMPUTATION -------------------------------------------------#
#############################################################################################################

## Use Will Rayner's tool
# see : https://www.well.ox.ac.uk/~wrayner/tools/Post-Imputation.html

input=${output_dir}/06_imputation_get
if [ ! -d "${output_dir}/07_imputation_check" ] ; then
  mkdir ${output_dir}/07_imputation_check
fi
output=${output_dir}/07_imputation_check

perl $bin/vcfparse.pl -d $input -o $output
cd $bin/IC
perl ic.pl -d $output -r $ref -o $output


#java -Xmx4g -classpath Junk.jar uk.ac.ox.well.t2d.reports.imputation.Main $output $output/summaryOutput
