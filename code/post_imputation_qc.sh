#!/bin/bash
# redirect stdout/stderr to a file
exec &> analysis/QC/post_imputation_qc.log

##SBATCH --job-name=QC_part2                                        # Job name (this name will appear in the queue)
##SBATCH --workdir=/data/sgg2/jenny/projects/PSYMETAB_GWAS          # The Working Directory of the job
##SBATCH --mail-type=ALL                                            # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=jennysjaarda@gmail.com                         # Where to send mail
##SBATCH --ntasks=1                                                 # Run on a single core
##SBATCH --time=24:00:00                                            # Time limit d-hrs:min:sec
##SBATCH --output=pipeline/QC_part2.out                             # Standard output and error log (%j: job allocation number)
##SBATCH --account=sgg                                              # runs on the sggg nodes.
##SBATCH --partition=sgg                                            # runs on the sggg nodes


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
#---------------------------------- PLINK CONVERSION  ------------------------------------------------------#
#############################################################################################################

## restrict to info/R2 >0.3
if [ ! -d "${output_dir}/08_plink_convert" ] ; then
  mkdir ${output_dir}/08_plink_convert
fi
cd ${output_dir}/08_plink_convert

for chr in {1..22}
do
  if [ ! -d "chr${chr}" ] ; then
    mkdir chr${chr}
  fi

  input_filename=../06_imputation_get/chr${chr}.dose.vcf.gz
  input_basename=chr${chr}

  plink2 --vcf ../06_imputation_get/chr${chr}.dose.vcf.gz  dosage=HDS \
    --exclude-if-info "R2<=0.3"	\
    --threads 1 \
    --double-id \
    --make-pgen \
    --out chr$chr/${input_basename}_test

  plink2 --pfile chr$chr/$input_basename  \
     --extract-if "R2 >= 0.1" \
     --threads 1 \
     --double-id \
     --make-pgen \
     --out chr$chr/${input_basename}_test

  ### update bim file to include rsIDs instead of chr:bp
  Rscript $project_dir/code/qc/update_pvar.R $data/dbSNP/dbSNP_SNP_list_chr${chr}.txt chr$chr/${input_basename}.pvar chr${chr}_update_name.txt chr$chr

  plink2 --pfile chr$chr/$input_basename  \
    --update-name chr$chr/chr${chr}_update_name.txt \
    --make-pgen \
    --threads 16 \
    --out chr$chr/${input_basename}_update_ID

done


#############################################################################################################
#---------------------------------- EXTRACT TYPED SNPS -----------------------------------------------------#
#############################################################################################################


## hitchhikers guide:
#https://rstudio-pubs-static.s3.amazonaws.com/452627_d519d1c86bd249e6a2d9638ef1ea836c.html

if [ ! -d "${output_dir}/09_extract_typed" ] ; then
  mkdir ${output_dir}/09_extract_typed
fi
cd ${output_dir}/09_extract_typed

## extract typed SNPs
for chr in {1..22}
do
  if [ ! -d "chr${chr}" ] ; then
    mkdir chr${chr}
  fi
  #### create a set of only genotyped data to run on snpweights
  plink2 --pfile ../08_plink_convert/chr$chr/chr${chr}_update_ID \
  --require-info "TYPED" \
  --threads 16 \
  --double-id \
  --make-pgen \
  --out chr$chr/temp1

done

## hard call back to bfile format
for chr in {1..22}
do
  plink2 --pfile chr$chr/temp1 \
  --hard-call-threshold 0.1  \
  --threads 16 \
  --double-id \
  --make-bed \
  --out chr$chr/${output_name}.chr${chr}.typed
done

# Write a list of the 22 autosomal filesets to be merged...
if [ -e "fileset.txt" ] ; then
  rm fileset.txt
fi

for chr in $(seq 1 22)
do
  echo chr$chr/${output_name}.chr${chr}.typed >> fileset.txt
done

# ...then merge them
plink --merge-list fileset.txt \
 --make-bed \
 --threads 16 \
 --out ${output_name}.typed

#############################################################################################################
#---------------------------- MERGE INPUTED DATA ACROSS CHROMOSOMES ----------------------------------------#
#############################################################################################################

## hitchhikers guide:
#https://rstudio-pubs-static.s3.amazonaws.com/452627_d519d1c86bd249e6a2d9638ef1ea836c.html


if [ ! -d "${output_dir}/10_merge_imputed/" ] ; then
 mkdir ${output_dir}/10_merge_imputed/
fi
cd ${output_dir}/10_merge_imputed/

# convert to vcf because there is no merge function in plink
for chr in $(seq 1 22)
do
 if [ ! -d "chr${chr}" ] ; then
   mkdir chr${chr}
 fi
 plink2 --pfile ../08_plink_convert/chr$chr/chr${chr}_update_ID \
  --recode vcf id-paste=iid vcf-dosage=HDS \
  --out chr${chr}/${output_name}.chr$chr \
  --threads 16
done

# Write a list of the 22 autosomal filesets to be merged...
if [ -e "fileset.txt" ] ; then
 rm fileset.txt
fi

for chr in $(seq 1 22)
do
 # Write a list of the 22 autosomal filesets to be merged...
 echo chr${chr}/${output_name}.chr${chr}.vcf >> fileset.txt
done

# ...then merge them
bcftools concat --file-list fileset.txt \
  --threads 16 \
  --output ${output_name}

# convert back to pgen
plink2 --vcf ${output_name} dosage=HDS \
  --threads 16 \
  --double-id \
  --make-pgen \
  --out ${output_name}

mv ${output_name} ${output_name}.vcf

#############################################################################################################
#---------------------------------- RELATEDNESS ------------------------------------------------------------#
#############################################################################################################

if [ ! -d "${output_dir}/11_relatedness" ] ; then
  mkdir ${output_dir}/11_relatedness
fi
cd ${output_dir}/11_relatedness


# Compute missingness statistics on this genome-wide data set...
plink2 --pfile ../10_merge_imputed/${output_name} \
  --missing \
  --threads 16 \
  --out ${output_name}_missing
# ...and check if there are any subjects missing more than 10% of total imputed genotypes
awk '$5 >= 0.1' ${output_name}_missing.smiss | awk '{print $1,$2}' > ${output_name}.filter.smiss
awk '$5 >= 0.1' ${output_name}_missing.vmiss | awk '{print $2}' > ${output_name}.filter.vmiss

plink2 --pfile ../10_merge_imputed/${output_name} \
  --exclude ${output_name}.filter.vmiss \
  --remove ${output_name}.filter.smiss \
  --make-king-table \
  --threads 16 \
  --out ${output_name}

wc -l ${output_name}.filter.vmiss
# 1184213 (minus 1 because of header) variants removed due to missingness
wc -l ${output_name}.filter.smiss
# 1 (minus 1 because of header) individuals removed due to missignness (no individuals)

awk '$8 >= 0.0884' ${output_name}.kin0 > ${output_name}.related.filter.smiss
awk '$8 >= 0.35' ${output_name}.kin0 > ${output_name}.duplicates

## A set of participants 17 was found to be duplicated that were not expected.
## After researching with Celine, it was decided to delete the following:

## These pairs of individuals were found to be in both studies under different IDs (PSYMETAB and Severines):
## We will keep the one on the right, but note that the code on the left is the same person and should be used for Severine

# 139CSM	HLHILECG
# 184CSM	WXPQPUWI
# 198CSM	CHMLAUKH
# 202CSM	OQLLRBMM
# 274UAS	YGGZRBYJ
# JO001	UUHDCMNL
# S008	YLXVRYDT
# S092	0021GE

## These pairs of individuals are duplicates but with different codes and we are keeping the one on the right:

# ABDJGAXW	ZXDVAJUD
# JO280	JO276
# S093	JO426

## These pair of individuals are not the same and the genetic info matches the one on the right (remove left):

# JNQSZGGL	XTVNZTRY
# WIELRZDD	FJNJEXCM
# OINDBCQM	SONUEGRB

## These pair of individuals could not be properly sorted and are therefore both being deleted:

# S060	0062GE
# 0030GE	SJJCZOXT

#####################################################################
### create file of duplicates to be removed  ########################
#####################################################################

if [ -e "duplicate_ids_remove.txt" ] ; then
  rm duplicate_ids_remove.txt
fi

for id in 139CSM 184CSM 198CSM 202CSM 274UAS JO001 S008 S092 ABDJGAXW JO280 S093 JNQSZGGL WIELRZDD OINDBCQM S060 0062GE 0030GE SJJCZOXT
do
echo ${id} >> duplicate_ids_remove.txt
done

if [ -e "matching_ids.txt" ] ; then
  rm matching_ids.txt
fi

echo "139CSM  HLHILECG" >> matching_ids.txt
echo "184CSM  WXPQPUWI" >> matching_ids.txt
echo "198CSM  CHMLAUKH" >> matching_ids.txt
echo "202CSM  OQLLRBMM" >> matching_ids.txt
echo "274UAS  YGGZRBYJ" >> matching_ids.txt
echo "JO001 UUHDCMNL" >> matching_ids.txt
echo "S008  YLXVRYDT" >> matching_ids.txt
echo "S092  0021GE" >> matching_ids.txt


if [ -e "${output_name}.duplicates2.txt" ] ; then
  rm ${output_name}.duplicates2.txt
fi

num_dups=$(wc -l <  duplicate_ids_remove.txt)
for i in $(seq 1 $num_dups)
do
  id=$(awk -v var=$i 'NR == var {print $1}' duplicate_ids_remove.txt)
  grep  $id.*$id ../10_merge_imputed/${output_name}.psam >> ${output_name}.duplicates2.txt
done

awk -v OFS=' ' '{print $1, $2}' ${output_name}.duplicates2.txt > ${output_name}.duplicates2.FID_IID.txt

cat ${output_name}.filter.smiss ${output_name}.duplicates2.FID_IID.txt > ${output_name}.filter ## create filtered file for missingness and including duplicate IDs

plink2 --pfile ../10_merge_imputed/${output_name} \
  --remove ${output_name}.duplicates2.FID_IID.txt \
  --make-pgen \
  --threads 16 \
  --out ${output_name}_remove_dups

plink2 --pfile ${output_name}_remove_dups \
  --exclude ${output_name}.filter.vmiss \
  --remove ${output_name}.filter.smiss \
  --make-king-table \
  --threads 16 \
  --out ${output_name}2

awk '$8 >= 0.0884' ${output_name}2.kin0 > ${output_name}2.related.filter.smiss ## 52 pairs above this threshold
awk '$8 >= 0.35' ${output_name}2.kin0 > ${output_name}2.duplicates ## there should be none now

### R script determines maximal set of unrelated individuals based on kinship results
cd $project_dir
Rscript $project_dir/code/qc/relatedness_filter.R "${output_name}2.kin0" 0.0884 "${output_dir}/11_relatedness" "${output_name}_related_ids.txt"

## remove related individuals
cd ${output_dir}/11_relatedness
plink2 --pfile ${output_name}_remove_dups \
  --remove ${output_name}_related_ids.txt \
  --threads 16 \
  --make-pgen \
  --out ${output_name}.relatedness.step11

# if using some mixed model design you can include all the related IDs, but remove the problematic duplicates

#############################################################################################################
#---------------------------------- ETHNICITY CHECK + ADMIXTURE ESTIMATION ---------------------------------#
#############################################################################################################

# plink projection,see: https://groups.google.com/forum/#!searchin/plink2-users/king$20%7Csort:date/plink2-users/W6DL5-hs_Q4/UKEL6wiCBAAJ
if [ ! -d "${output_dir}/12_ethnicity_admixture" ] ; then
  mkdir ${output_dir}/12_ethnicity_admixture
fi
cd ${output_dir}/12_ethnicity_admixture

plink2 --bfile ../09_extract_typed/${output_name}.typed \
  --remove ../11_relatedness/${output_name}.duplicates2.FID_IID.txt \
  --threads 16 \
  --make-bed \
  --out ${output_name}.typed_remove_dups

plink2 --bfile ${output_name}.typed_remove_dups \
  --geno 0.1 \
  --maf 0.05 \
  --hwe 5e-4 \
  --threads 16 \
  --make-bed \
  --out ${output_name}.typed.QC

## snpweights calculation
mkdir snpweights

# 249654 variants and 2715 people pass filters and QC.
sh $project_dir/code/qc/snpweight_conversion/fam2ind.sh ${output_name}.typed.QC.fam
sh $project_dir/code/qc/snpweight_conversion/bim2snp.sh ${output_name}.typed.QC.bim
sh $project_dir/code/qc/snpweight_conversion/bed2geno.sh ${output_name}.typed.QC.bed

if [ -f "inferancestry.par" ] ; then
  rm snpweights/inferancestry.par
fi
echo geno:  ${output_name}.typed.QC.geno >> snpweights/inferancestry.par
echo snp:   ${output_name}.typed.QC.snp >> snpweights/inferancestry.par
echo ind:   ${output_name}.typed.QC.ind >> snpweights/inferancestry.par
echo snpwt: $bin/SNPweights2.1/snpwt.NA  >> snpweights/inferancestry.par
echo predpcoutput: snpweights/${output_name}.NA.predpc >> snpweights/inferancestry.par

# The 10 columns in the ${output_name}.NA.predpc output file are: sample ID, population label,
# number of SNPs used for inference, predicted PC1, predicted PC2, predicted PC3, % YRI ancestry,
# % CEU ancestry, % East Asian ancestry, % Native American ancestry.

# In the second stage of this analysis, we will focus on subjects with >80% predicted CEU ancestry,
# and determine the extents to which these subjects fall into the three following European subpopulations:
# North-Western Europeans (NWE), South-Eastern Europeans (SEE) and Ashkenazi Jewish (AJW).
# This step will also produce the principal components that will be used as covariates in the GWAS to correct for
# population stratification.

# So first let us write a list of CEU subjects to be retained,
# subset our genotype data, and repeat the conversion to EIGENSTRAT format:

$bin/SNPweights2.1/bin/inferanc -p snpweights/inferancestry.par

awk '$8 >= 0.8' snpweights/${output_name}.NA.predpc | awk '{printf $1"\t"$2"\n"}' > snpweights/${output_name}.CEU80
plink --bfile ${output_name}.typed.QC \
  --keep snpweights/${output_name}.CEU80 \
  --make-bed \
  --out snpweights/${output_name}.CEU80
# 249634 variants and 2170 people pass filters and QC.
# 2697 - 2154 = 543 people removed as non Europeans

cd snpweights
sh $project_dir/code/qc/snpweight_conversion/fam2ind.sh ${output_name}.CEU80.fam
sh $project_dir/code/qc/snpweight_conversion/bim2snp.sh ${output_name}.CEU80.bim
sh $project_dir/code/qc/snpweight_conversion/bed2geno.sh ${output_name}.CEU80.bed

cd ..

if [ -f "snpweights/inferancestry_sub.par" ] ; then
  rm snpweights/inferancestry_sub.par
fi
echo geno:  snpweights/${output_name}.CEU80.geno >> snpweights/inferancestry_sub.par
echo snp:   snpweights/${output_name}.CEU80.snp >> snpweights/inferancestry_sub.par
echo ind:   snpweights/${output_name}.CEU80.ind >> snpweights/inferancestry_sub.par
echo snpwt: $bin/SNPweights2.1/snpwt.EA  >> snpweights/inferancestry_sub.par
echo predpcoutput: snpweights/${output_name}.CEU80.EA.predpc >> snpweights/inferancestry_sub.par

$bin/SNPweights2.1/bin/inferanc -p snpweights/inferancestry_sub.par
rm tmp.*
rm snpweights/tmp.*

# The 8 columns in the ${output_name}.CEU80.EA.predpc output file are: sample ID,
# population label, number of SNPs used for inference, predicted PC1, predicted PC2,
# % Northwest European ancestry, % Southeast European ancestry, % Ashkenazi Jewish ancestry.

## PCA calculation

## first create pruned set of SNP
mkdir pca
plink2 --bfile ${output_name}.typed.QC \
 --indep-pairwise 50 5 0.2 \
 --remove ../11_relatedness/${output_name}_related_ids.txt \
 --threads 16 \
 --out pca/${output_name}.indep

# could exclude regions suggested by flaspca by adding the following line:
#--exclude range $flashpca/exclusion_regions_hg19.txt \

# wc -l pca/PSYMETAB_GWAS.QC.indep.prune.in
# leaves 96746 SNPs

## create an unrelated set for calculating PCs with same set of pruned SNPs
plink2 --bfile ${output_name}.typed.QC \
  --extract pca/${output_name}.indep.prune.in \
  --remove ../11_relatedness/${output_name}_related_ids.txt \
  --make-bed \
  --threads 16 \
  --out pca/${output_name}.pruned.unrelated

## create an full set for calculating projecting PCs with same set of pruned SNPs
plink2 --bfile ${output_name}.typed.QC \
  --extract pca/${output_name}.indep.prune.in \
  --make-bed \
  --threads 16 \
  --out pca/${output_name}.pruned.full
#--keep ../11_relatedness/${output_name}_related_ids.txt \ ## could create an unrelated set using keep

cd pca
$flashpca/flashpca --bfile ${output_name}.pruned.unrelated \
  --suffix .${output_name}_unrelated.txt \
  --ndim 100 \
  --outload loadings.txt \
  --outmeansd meansd.txt \
  --numthreads 16

$flashpca/flashpca --bfile ${output_name}.pruned.full \
  --project \
  --inmeansd meansd.txt \
  --outproj ${output_name}_projections.txt \
  --inload loadings.txt -v \
  --numthreads 16

## create ethnicity plots
cd $project_dir
Rscript $project_dir/code/qc/ethnicity_check.R ${output_name} ${output_dir}/12_ethnicity_admixture/pca \
  ${output_dir}/12_ethnicity_admixture/pca/${output_name}_projections.txt \
  ${output_dir}/12_ethnicity_admixture/snpweights/${output_name}.NA.predpc ${eth_file}


#############################################################################################################
#------------------------------------------------ HWE CHECK ------------------------------------------------#
#############################################################################################################


if [ ! -d "${output_dir}/13_hwecheck" ] ; then
  mkdir ${output_dir}/13_hwecheck
fi
cd ${output_dir}/13_hwecheck

if [ -e "${output_name}.hardy.sig.all" ] ; then
  rm ${output_name}.hardy.sig.all
fi

plink2 --pfile ../10_merge_imputed/${output_name} \
  --make-pfile \
  --remove ../11_relatedness/${output_name}.duplicates2.FID_IID.txt \
  --out ${output_name}.remove_dups \
  --threads 16

for eth in CEU EA MIXED NA YRI ; do

  count=$(wc -l < "../12_ethnicity_admixture/pca/${output_name}_${eth}_samples.txt")

  #1. Separate files by ethnicity
  if [ ! -d "${eth}" ] ; then
    mkdir ${eth}
  fi

  plink2 --pfile ${output_name}.remove_dups \
    --make-pfile \
    --keep ../12_ethnicity_admixture/pca/${output_name}_${eth}_samples.txt \
    --remove ../11_relatedness/${output_name}_related_ids.txt \
    --out $eth/${output_name}.${eth} \
    --threads 16

  #2. Evaluate HWE for each variant
  plink2 --pfile $eth/${output_name}.${eth}  \
    --hardy \
    --out $eth/${output_name}.${eth}  \
    --threads 16

  awk '$10 < 0.
  {print $2}' $eth/${output_name}.${eth}.hardy | sort -u > $eth/${output_name}.${eth}.hardy.sig

  # If size of ethnic group is < 100, don't perform HWE check in that ethnic group
  if [ "$count" -ge 100 ] ; then
    cat $eth/${output_name}.${eth}.hardy.sig >> ${output_name}.hardy.sig.all
  fi

done

for eth in CEU EA MIXED NA YRI ; do

  count=$(wc -l < "../12_ethnicity_admixture/pca/${output_name}_${eth}_samples.txt")
  echo $eth":" $count

done

# CEU: 2154
# EA: 29
# MIXED: 426
# NA: 5
# YRI: 83

#3. Identify SNPs failing HWE in at least one ethnic group
sort -u ${output_name}.hardy.sig.all > ${output_name}.hardy.sig.unique
#12750 variants fail HWE

#4. Remove HWE violations
plink2 --pfile ${output_name}.remove_dups \
  --make-pfile \
  --remove ${output_name}.hardy.sig.unique \
  --out ${output_name}.hwecheck.step13 \
  --threads 16

#############################################################################################################
#------------------------------------------------ MAF CHECK ------------------------------------------------#
#############################################################################################################


if [ ! -d "${output_dir}/14_mafcheck" ] ; then
  mkdir ${output_dir}/14_mafcheck
fi
cd ${output_dir}/14_mafcheck

### identify SNPs which fail MAF in all ethnic groups

for eth in CEU EA MIXED NA YRI ; do

  count=$(wc -l < "../12_ethnicity_admixture/pca/${output_name}_${eth}_samples.txt")

  #1. Separate files by ethnicity
  if [ ! -d "${eth}" ] ; then
    mkdir ${eth}
  fi

  #2. identify freq of variants in each eth
  plink2 --pfile ../13_hwecheck/${output_name}.hwecheck.step13 \
    --keep ../12_ethnicity_admixture/pca/${output_name}_${eth}_samples.txt \
    --freq \
    --make-pfile \
    --out $eth/${output_name}.${eth} \
    --threads 16

  #awk '$5 < 0.01 || $5 > 0.99{print $0}' $eth/${output_name}.${eth}.afreq > $eth/${output_name}.${eth}.afreq.rare
  #awk '{print $2}' $eth/${output_name}.${eth}.afreq.rare | sort -u > $eth/${output_name}.${eth}.rare.snps

  # If size of ethnic group is < 100, don't perform HWE check in that ethnic group
  #if [ "$count" -ge 100 ] ; then
  #  cat $eth/${output_name}.${eth}.rare.snps >> ${output_name}.rare.snps.all
  #fi

done

cd $project_dir
Rscript $project_dir/code/qc/maf_check.R ${output_name} ${output_dir}/14_mafcheck/ETHNICITY_NAME/${output_name}.ETHNICITY_NAME.afreq \
  ${output_dir}/14_mafcheck/ETHNICITY_NAME/${output_name}.ETHNICITY_NAME.psam \
  ${output_dir}/14_mafcheck 0.01 100

cd ${output_dir}/14_mafcheck
plink2 --pfile ../13_hwecheck/${output_name}.hwecheck.step13 \
  --make-pfile \
  --remove ${output_name}_low_maf_snps.txt \
  --out ${output_name}.mafcheck.step14 \
  --threads 16
