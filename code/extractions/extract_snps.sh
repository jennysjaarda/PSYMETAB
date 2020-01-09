#!/bin/bash
snp_file=${1}
output=${2}
qc_data=${3}
pc_data=${4} #$output/${out_name}/

file_start=$(echo $snp_file | cut -f1 -d.)
out_name=$(basename $file_start)
mkdir $output/$out_name
clean_file=$output/${out_name}/${out_name}_rsids.txt
snp_info=$output/${out_name}/${out_name}_SNP_info.txt

  #clean input file
tr -d '\r' < $file | sed '/^$/d' > $clean_file

num=$(awk 'NR>0 {ORS=" ";print NR}' $clean_file)
awk 'FNR==NR {a[$1]; next}; $3 in a' $clean_file $data/HRC.r1-1.GRCh37.wgs.mac5.sites.tab > $snp_info


for chr in {1..22}
do

  awk -v var="$chr " '{ if ($1 == var) { print $3 } }' $snp_info > $output/$out_name/${out_name}_chr${chr}.txt
  num_snps=$(cat $output/$out_name/${out_name}_chr${chr}.txt | wc -l)
  if [ ${num_snps} != 0 ]; then
    plink2 --pfile $qc_data --extract $output/${out_name}/${out_name}_chr${chr}.txt --recode A --out $output/${out_name}/${out_name}_chr${chr}_extract
  fi
  if [ ${num_snps} == 0 ]; then
    rm $output/$out_name/${out_name}_chr${chr}.txt
  fi
done

Rscript $PSYMETAB/clean_extraction.r $output $outname $pc_data
