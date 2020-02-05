
# !/bin/bash
# redirect stdout/stderr to a file
exec &> analysis/QC/pre_imputation_qc.log

#
## SBATCH --job-name=QC_part1                                        # Job name (this name will appear in the queue)
## SBATCH --workdir=/data/sgg2/jenny/projects/PSYMETAB_GWAS          # The Working Directory of the job
## SBATCH --mail-type=ALL                                            # Mail events (NONE, BEGIN, END, FAIL, ALL)
## SBATCH --mail-user=jennysjaarda@gmail.com                         # Where to send mail
## SBATCH --ntasks=1                                                 # Run on a single core
## SBATCH --time=24:00:00                                            # Time limit d-hrs:min:sec
## SBATCH --output=pipeline/QC_part1.out                             # Standard output and error log (%j: job allocation number)
## SBATCH --account=sgg                                              # runs on the sggg nodes.
## SBATCH --partition=sgg                                            # runs on the sggg nodes

####################################################################################################################
####################################################################################################################
####################################################################################################################
####### PSYMETAB GWAS  #############################################################################################
####### QUALITY CONTROL ############################################################################################
####### Genotype Chromosome information: \PCN\UBPC\ANALYSES_RECHERCHE\Jenny\PSYMETAB\PLINK_091019_0920   ######
####### Phenotype data:                  \PCN\UBPC\ANALYSES_RECHERCHE\Jenny\PSYMETAB\GWAS_PHENO          ######
####### AUTHOR: Jenny Sjaarda ######################################################################################
####################################################################################################################
####################################################################################################################

### DEFINE INPUTS

project_dir="/data/sgg2/jenny/projects/PSYMETAB"
plink_data="PLINK_091019_0920"
raw_data=$project_dir/data/raw/$plink_data

output_dir=$project_dir/analysis/QC
input_chip=$project_dir/data/processed/${plink_data}/PSYMETAB_GWAS
output_name=PSYMETAB_GWAS
sex_file=$project_dir/data/processed/phenotype_data/PSYMETAB_GWAS_sex.txt
eth_file=$project_dir/data/processed/phenotype_data/PSYMETAB_GWAS_eth.txt
duplicate_file=$project_dir/data/processed/phenotype_data/PSYMETAB_GWAS_dupIDs.txt
duplicate_file_set=$project_dir/data/processed/phenotype_data/PSYMETAB_GWAS_dupIDs_set.txt
strand_file=$project_dir/data/processed/GSAMD_UPPC-b37.strand
conversion_file=$project_dir/data/processed/reference_files/rsid_conversion.txt
KG_PRUNED="$x1000G/1000G_chrall.ex_maf_05.pruned"
KG_SAMPLES="$x1000G/integrated_call_samples_v3.20130502.ALL.panel"
ref=$data/HRC.r1-1.GRCh37.wgs.mac5.sites.tab

#############################################################################################################
#---------------------------------------------- PREPROCESSING ----------------------------------------------#
#############################################################################################################

#################################################################
#################################################################

if [ ! -d "${output_dir}" ] ; then
  mkdir ${output_dir}
fi
if [ ! -d "${output_dir}/00_preprocessing" ] ; then
  mkdir ${output_dir}/00_preprocessing
fi
cd ${output_dir}/00_preprocessing

# i) Create binary PLINK files (if necessary) 761641 variants to start
# wc -l $input_chip.bim
if [ ! -s "$input_chip.bed" ] || [ ! -s "$input_chip.bim" ] || [ ! -s "$input_chip.fam" ] ; then
plink --file $input_chip \
  --make-bed \
  --threads 16 \
  --out $input_chip
fi

# ii) exclude 11641 Y and MT variants -->  750000 variants remaining
# wc -l non_autosome_or_x_variants
awk '($1==0 || $1>23){print $2}' $input_chip.bim > non_autosome_or_x_variants
plink --bfile $input_chip \
  --exclude non_autosome_or_x_variants \
  --make-bed \
  --threads 16 \
  --out temp1

# iii) remove 7945 duplicate variants --> 742055 variants remaining
# duplicates by chr pos ref alt
# wc -l duplicate_variants
awk '{dup[$1,$4,$5,$6]++}(dup[$1,$4,$5,$6]>1 || dup[$1,$4,$6,$5]>1){print $2}' temp1.bim > duplicate_variants
plink --bfile temp1 \
  --exclude duplicate_variants \
  --make-bed \
  --threads 16 \
  --out temp2

# iv) Update sex
plink --bfile temp2 \
  --update-sex $sex_file \
  --make-bed \
  --threads 16 \
  --out temp3

# v) remove duplicates --> 2752 individuals remaining
# wc -l ${output_name}.preprocessing.step0.fam
plink --bfile temp3 \
  --remove $duplicate_file \
  --make-bed \
  --threads 16 \
  --out ${output_name}.preprocessing.step0


## sanity check: ensure that duplicates have the same genetic info
plink2 --bfile temp3 \
  --keep $duplicate_file_set \
  --make-bed \
  --out ${output_name}.duplicate.set \
  --threads 16

plink2 --bfile ${output_name}.duplicate.set \
  --make-king-table \
  --threads 16 \
  --out temp4

if [ -e "PSYMETAB_GWAS.duplicates.kin0" ] ; then
  rm ${output_name}.duplicates.kin0
fi

num_dups=$(wc -l < $duplicate_file)
for i in $(seq 1 $num_dups)
do
  id_dup=$(awk -v var=$i 'NR == var {print $2}' $duplicate_file)
  id=${id_dup::-3}
  grep  $id.*$id temp4.kin0 >> ${output_name}.duplicates.kin0
done

#############################################################################################################
#--------------------------------------------- STRAND ALIGNMENT --------------------------------------------#
#############################################################################################################
if [ ! -d "${output_dir}/01_strandflip" ] ; then
  mkdir ${output_dir}/01_strandflip
fi
cd ${output_dir}/01_strandflip
#Cut the strand file into a series of plink slices
strand_file_name=`basename $strand_file`
chr_file=${output_dir}/${strand_file_name}.chr
pos_file=${output_dir}/${strand_file_name}.pos
flip_file=${output_dir}/${strand_file_name}.flip
cat $strand_file | cut -f 1,2 > $chr_file
cat $strand_file | cut -f 1,3 > $pos_file
cat $strand_file | awk '{if ($5=="-") print $0}' | cut -f 1 > $flip_file

#1. Apply the chr
plink --bfile ../00_preprocessing/${output_name}.preprocessing.step0 \
  --update-map $chr_file \
  --update-chr \
  --make-bed \
  --out temp1

#2. Apply the pos
plink --bfile temp1 \
  --update-map $pos_file \
  --make-bed \
  --out temp2

#3. Apply the flip
plink --bfile temp2 \
  --flip $flip_file \
  --make-bed \
  --out temp3

#4. Extract the SNPs in the pos file, we don't want SNPs that aren't in the strand file
plink --bfile temp3 \
  --extract $pos_file \
  --make-bed \
  --out temp4

#5. Remove any non autosome / X chromosomal variants (might have changed due to the pos and chr update)
awk '($1==0 || $1>23){print $2}' temp4.bim > non_autosome_or_x_variants
plink --bfile temp4 \
  --exclude non_autosome_or_x_variants \
  --make-bed \
  --out temp5

#6. Update rsids
plink --bfile temp5 \
  --update-map $conversion_file \
  --update-name \
  --make-bed \
  --out temp6

#7. Remove duplicates
awk 'NR==FNR{count[$2]++}count[$2]>1{print $2}' temp6.bim > duplicate_variants
awk 'NR==FNR{count[$1,$4,$5,$6]++;count[$1,$4,$6,$5]++}(count[$1,$4,$5,$6]>1) || (count[$1,$4,$6,$5]>1){print $2}' temp6.bim >> duplicate_variants
sort -u duplicate_variants > duplicate_variants.sort
# wc -l duplicate_variants: 6767 variants
plink --bfile temp6 \
  --exclude duplicate_variants.sort \
  --make-bed \
  --out ${output_name}.strandflip.step1

# 742055-734328 variants lost to strand flip - 734328 variants remain

#############################################################################################################
#-------------------------------------------------- MAF ZERO -----------------------------------------------#
#############################################################################################################
# Removal of SNPs that have MAF-ZERO after missingness checks
if [ ! -d "${output_dir}/02_maf_zero" ] ; then
mkdir ${output_dir}/02_maf_zero
fi
cd ${output_dir}/02_maf_zero
#1. Calculate frequency
plink --bfile ../01_strandflip/${output_name}.strandflip.step1 \
  --freq \
  --out ${output_name}

#2. Write MAF=0 SNPs to a text file
awk '$5==0{print}' ${output_name}.frq > ${output_name}.maf_zero_snp.txt
# wc -l ${output_name}.maf_zero_snp.txt: 93578 variants MAF=0
#3. Remove
plink --bfile ../01_strandflip/${output_name}.strandflip.step1 \
  --exclude ${output_name}.maf_zero_snp.txt \
  --make-bed \
  --out ${output_name}.maf_zero.step2

# 93578 variants lost to maf zero filter ... 640750 variants remaining

#############################################################################################################
#------------------------------------------------ MISSINGNESS ----------------------------------------------#
#############################################################################################################
if [ ! -d "${output_dir}/03_missingness" ] ; then
mkdir ${output_dir}/03_missingness
fi
cd ${output_dir}/03_missingness

#1. Exclude Variants with >10% missingness
# 6480 SNPs failed missingness test ( GENO > 0.1 )
# 634270 SNPs remain
plink --bfile ../02_maf_zero/${output_name}.maf_zero.step2 \
  --geno 0.1 \
  --make-bed \
  --out ${output_name}.missingness_geno_10

#2. Exclude Individuals with >10% missingness
# 3 people fail missingness test ( MIND > 0.1 )
# 2749 people remain
plink --bfile ${output_name}.missingness_geno_10 \
  --mind 0.10 \
  --make-bed \
  --out ${output_name}.missingness_mind_10

#3. Exclude Variants with >5% missingness
# 8256 SNPs failed missingness test ( GENO > 0.05 )
# 626014 SNPs remain
plink --bfile ${output_name}.missingness_mind_10 \
  --geno 0.05 \
  --make-bed \
  --out ${output_name}.missingness_geno_05

#4. Exclude Individuals with >5% missingness
# 3 people fail missingness test ( MIND > 0.05 )
# 2746 people remain
plink --bfile ${output_name}.missingness_geno_05 \
  --mind 0.05 \
  --make-bed \
  --out ${output_name}.missingness_mind_05

#5. Exclude Variants with >1% missingness
# 35957 SNPs failed missingness test ( GENO > 0.01 )
# 590057 SNPs remain
plink --bfile ${output_name}.missingness_mind_05 \
  --geno 0.01 \
  --make-bed \
  --out ${output_name}.missingness_geno_01

#6. Exclude Individuals with >1% missingness
# 5 people fail missingness test ( MIND > 0.05 )
# 2741 people remain
plink --bfile ${output_name}.missingness_geno_01 \
  --mind 0.01 \
  --make-bed \
  --out ${output_name}.missingness.step3

###### TOTAL REMOVED: 3+3+5=11 (11/2752 = 0.40%) SAMPLES &
#################### 6480+8256+35957=50693 (50283/640750 = 7.8%) VARIANTS
###### TOTAL RETAINED: 2752 SAMPLES & 590057 variants

#############################################################################################################
#------------------------------------------------- SEX CHECK -----------------------------------------------#
#############################################################################################################
if [ ! -d "${output_dir}/04_sexcheck" ] ; then
  mkdir ${output_dir}/04_sexcheck
fi
cd ${output_dir}/04_sexcheck

#1. Perform sex check
plink --bfile ../03_missingness/${output_name}.missingness.step3 \
  --check-sex \
  --out ${output_name}.sexcheck.step4

#2. Remove unambiguous sex violations
grep -w "PROBLEM" ${output_name}.sexcheck.step4.sexcheck | awk '($4!="0"){print $1,$2}' > ${output_name}.sexcheck.problem.ids
#wc -l ${output_name}.sexcheck.problem.ids: 26 sex problems
plink --bfile ../03_missingness/${output_name}.missingness.step3 \
  --remove ${output_name}.sexcheck.problem.ids \
  --make-bed \
  --out ${output_name}.sexcheck.step4

# 26 (0.94%) sex errors
# 2715 individuals retained

#############################################################################################################
#---------------------------------- IMPUTATION PRERPARATION -------------------------------------------------#
#############################################################################################################

if [ ! -d "${output_dir}/05_imputation_prep" ] ; then
  mkdir ${output_dir}/05_imputation_prep
fi
cd ${output_dir}/05_imputation_prep

cp ../04_sexcheck/${output_name}.sexcheck.step4.bim ${output_dir}/05_imputation_prep/${output_name}.geno_QC.output.bim
cp ../04_sexcheck/${output_name}.sexcheck.step4.bed ${output_dir}/05_imputation_prep/${output_name}.geno_QC.output.bed
cp ../04_sexcheck/${output_name}.sexcheck.step4.fam ${output_dir}/05_imputation_prep/${output_name}.geno_QC.output.fam

plink --freq --bfile ${output_dir}/05_imputation_prep/${output_name}.geno_QC.output --out ${output_name}.geno_QC.output

perl $bin/HRC-1000G-check-bim-NoReadKey.pl -b ${output_name}.geno_QC.output.bim -f ${output_name}.geno_QC.output.frq \
-r $ref -h

## Produces:
# A set of plink commands to update or remove SNPs based on the checks
# as well as a file (FreqPlot) of cohort allele frequency vs reference
# panel allele frequency.

### replace underscore since vcf conversion uses underscores between FID and IID
tr '_' '-' < ${output_name}.geno_QC.output.fam > ${output_name}.geno_QC.output.fam.temp
mv ${output_name}.geno_QC.output.fam ${output_name}.geno_QC.output.fam.original
mv ${output_name}.geno_QC.output.fam.temp ${output_name}.geno_QC.output.fam

sh Run-plink.sh

for chr in {1..22}
do
  vcf-sort ${output_name}.geno_QC.output-updated-chr${chr}.vcf | bgzip -c > ${output_name}.imputation_input.chr${chr}.vcf.gz
done


## Download zipped files to personal computation and copy to Michigan Imputation Server:

# /data/sgg2/jenny/projects/PSYMETAB_GWAS/pipeline/QC/05_imputation_prep/PSYMETAB_GWAS.imputation_input.chr1.vcf.gz
# /data/sgg2/jenny/projects/PSYMETAB_GWAS/pipeline/QC/05_imputation_prep/PSYMETAB_GWAS.imputation_input.chr2.vcf.gz
# /data/sgg2/jenny/projects/PSYMETAB_GWAS/pipeline/QC/05_imputation_prep/PSYMETAB_GWAS.imputation_input.chr3.vcf.gz
# /data/sgg2/jenny/projects/PSYMETAB_GWAS/pipeline/QC/05_imputation_prep/PSYMETAB_GWAS.imputation_input.chr4.vcf.gz
# /data/sgg2/jenny/projects/PSYMETAB_GWAS/pipeline/QC/05_imputation_prep/PSYMETAB_GWAS.imputation_input.chr5.vcf.gz
# /data/sgg2/jenny/projects/PSYMETAB_GWAS/pipeline/QC/05_imputation_prep/PSYMETAB_GWAS.imputation_input.chr6.vcf.gz
# /data/sgg2/jenny/projects/PSYMETAB_GWAS/pipeline/QC/05_imputation_prep/PSYMETAB_GWAS.imputation_input.chr7.vcf.gz
# /data/sgg2/jenny/projects/PSYMETAB_GWAS/pipeline/QC/05_imputation_prep/PSYMETAB_GWAS.imputation_input.chr8.vcf.gz
# /data/sgg2/jenny/projects/PSYMETAB_GWAS/pipeline/QC/05_imputation_prep/PSYMETAB_GWAS.imputation_input.chr9.vcf.gz
# /data/sgg2/jenny/projects/PSYMETAB_GWAS/pipeline/QC/05_imputation_prep/PSYMETAB_GWAS.imputation_input.chr10.vcf.gz
# /data/sgg2/jenny/projects/PSYMETAB_GWAS/pipeline/QC/05_imputation_prep/PSYMETAB_GWAS.imputation_input.chr11.vcf.gz
# /data/sgg2/jenny/projects/PSYMETAB_GWAS/pipeline/QC/05_imputation_prep/PSYMETAB_GWAS.imputation_input.chr12.vcf.gz
# /data/sgg2/jenny/projects/PSYMETAB_GWAS/pipeline/QC/05_imputation_prep/PSYMETAB_GWAS.imputation_input.chr13.vcf.gz
# /data/sgg2/jenny/projects/PSYMETAB_GWAS/pipeline/QC/05_imputation_prep/PSYMETAB_GWAS.imputation_input.chr14.vcf.gz
# /data/sgg2/jenny/projects/PSYMETAB_GWAS/pipeline/QC/05_imputation_prep/PSYMETAB_GWAS.imputation_input.chr15.vcf.gz
# /data/sgg2/jenny/projects/PSYMETAB_GWAS/pipeline/QC/05_imputation_prep/PSYMETAB_GWAS.imputation_input.chr16.vcf.gz
# /data/sgg2/jenny/projects/PSYMETAB_GWAS/pipeline/QC/05_imputation_prep/PSYMETAB_GWAS.imputation_input.chr17.vcf.gz
# /data/sgg2/jenny/projects/PSYMETAB_GWAS/pipeline/QC/05_imputation_prep/PSYMETAB_GWAS.imputation_input.chr18.vcf.gz
# /data/sgg2/jenny/projects/PSYMETAB_GWAS/pipeline/QC/05_imputation_prep/PSYMETAB_GWAS.imputation_input.chr19.vcf.gz
# /data/sgg2/jenny/projects/PSYMETAB_GWAS/pipeline/QC/05_imputation_prep/PSYMETAB_GWAS.imputation_input.chr20.vcf.gz
# /data/sgg2/jenny/projects/PSYMETAB_GWAS/pipeline/QC/05_imputation_prep/PSYMETAB_GWAS.imputation_input.chr21.vcf.gz
# /data/sgg2/jenny/projects/PSYMETAB_GWAS/pipeline/QC/05_imputation_prep/PSYMETAB_GWAS.imputation_input.chr22.vcf.gz

### IMPUTATION was run on michigian impuation server: https://imputationserver.sph.umich.edu/index.html#!
# using the HRC reference panel
