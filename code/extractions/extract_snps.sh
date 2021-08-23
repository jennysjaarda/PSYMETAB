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
snp_hrc_info=$output/${out_name}/${out_name}_SNP_HRC_info.txt
snp_impute_info=$output/${out_name}/${out_name}_SNP_impute_info.txt

  #clean input file
tr -d '\r' < $snp_file | sed '/^$/d' > $clean_file

# get CHR/BP information from HRC
num=$(awk 'NR>0 {ORS=" ";print NR}' $clean_file)
awk 'FNR==NR {a[$1]; next}; $3 in a' $clean_file $data/HRC.r1-1.GRCh37.wgs.mac5.sites.tab > $snp_hrc_info
sed  -i '1i #CHROM  POS ID  REF ALT AC  AN  AF  AC_EXCLUDING_1000G  AN_EXCLUDING_1000G  AF_EXCLUDING_1000G  AA' $snp_hrc_info

awk 'FNR==NR {a[$1]; next}; $3 in a' $clean_file $qc_data.pvar > $snp_impute_info
sed  -i '1i #CHROM  POS  ID REF ALT FILTER  INFO' $snp_impute_info

for chr in {1..22}
do


  awk -v var="$chr " '{ if ($1 == var) { print $3 } }' $snp_hrc_info > $output/$out_name/${out_name}_chr${chr}.txt
  num_snps=$(cat $output/$out_name/${out_name}_chr${chr}.txt | wc -l)
  if [ ${num_snps} != 0 ]; then
    plink2 --pfile $qc_data --extract $output/${out_name}/${out_name}_chr${chr}.txt --make-pfile --out $output/${out_name}/${out_name}_chr${chr}_extract
    plink2 --pfile $output/${out_name}/${out_name}_chr${chr}_extract --recode A --out $output/${out_name}/${out_name}_chr${chr}_extract
  fi
  if [ ${num_snps} == 0 ]; then
    rm $output/$out_name/${out_name}_chr${chr}.txt
  fi
done


Rscript code/extractions/pc_merge.r $output $out_name $pc_data $snp_file

if [ -s "$output/${out_name}/${out_name}_missing_snps.txt" ]; then
  plink2 --bfile $input_chip --extract $output/${out_name}/${out_name}_missing_snps.txt --make-bed --out $output/${out_name}/${out_name}_missing_snps_extract
  if [ -s "$output/${out_name}/${out_name}_missing_snps_extract.fam" ]; then
    plink2 --bfile $output/${out_name}/${out_name}_missing_snps_extract --recode A --out $output/${out_name}/${out_name}_missing_snps_extract
    Rscript code/extractions/geno_original_merge.r $output $out_name $snp_file
  fi
fi

if [ ! -f "$output/${out_name}/${out_name}_extraction_geno.txt" ]; then
  cp $output/${out_name}/${out_name}_extraction.txt $output/${out_name}/${out_name}_extraction_geno.txt
fi
