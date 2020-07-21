#!/bin/bash
snp_file=${1}
output=${2}
qc_data=${3}
pc_data=${4} #$output/${out_name}/
QC_dir=${5}
input_chip=${6}

file_start=$(echo $snp_file | cut -f1 -d.)
out_name=$(basename $file_start)

clean_file=$output/${out_name}/${out_name}_rsids.txt
snp_info=$output/${out_name}/${out_name}_SNP_info.txt

  #clean input file
tr -d '\r' < $snp_file | sed '/^$/d' > $clean_file

num=$(awk 'NR>0 {ORS=" ";print NR}' $clean_file)
awk 'FNR==NR {a[$1]; next}; $3 in a' $clean_file $data/HRC.r1-1.GRCh37.wgs.mac5.sites.tab > $snp_info


for chr in {1..22}
do


  awk -v var="$chr " '{ if ($1 == var) { print $3 } }' $snp_info > $output/$out_name/${out_name}_chr${chr}.txt
  num_snps=$(cat $output/$out_name/${out_name}_chr${chr}.txt | wc -l)
  if [ ${num_snps} != 0 ]; then
    plink2 --pfile $qc_data --extract $output/${out_name}/${out_name}_chr${chr}.txt --make-pfile --out $output/${out_name}/${out_name}_chr${chr}_extract
    plink2 --pfile $output/${out_name}/${out_name}_chr${chr}_extract --recode A --out $output/${out_name}/${out_name}_chr${chr}_extract
  fi
  if [ ${num_snps} == 0 ]; then
    rm $output/$out_name/${out_name}_chr${chr}.txt
  fi
done

for eth in CEU EA MIXED NA YRI ; do

  pc_data=${QC_dir}/15_final_processing/final_pca/ETH/pcs.PSYMETAB_GWAS_ETH_unrelated.txt
  if [ -f "$pc_data" ] ; then
    Rscript code/extractions/pc_merge.r $output $out_name $pc_data $snp_file
  fi
done

plink2 --bfile $input_chip --extract $output/${out_name}/${out_name}_missing_snps.txt --make-bed --out $output/${out_name}/${out_name}_missing_snps_extract
plink2 --bfile $output/${out_name}/${out_name}_missing_snps_extract --recode A --out $output/${out_name}/${out_name}_missing_snps_extract


if [ -s "$output/${out_name}/${out_name}_missing_snps_extract.raw" ]; then
  Rscript code/extractions/geno_original_merge.r $output $out_name $pc_data $snp_file
fi
