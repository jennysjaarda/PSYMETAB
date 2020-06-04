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

  #2. identify freq of variants in each eth excluding related individuals
  plink2 --pfile ../13_hwecheck/${output_name}.hwecheck.step13 \
    --keep ../12_ethnicity_admixture/pca/${output_name}_${eth}_samples.txt \
    --remove ../11_relatedness/${output_name}_related_ids.txt \
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
  --exclude ${output_name}_low_maf_snps.txt \
  --out ${output_name}.mafcheck.step14 \
  --threads 16
