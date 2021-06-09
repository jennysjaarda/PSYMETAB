#!/bin/bash

#SBATCH --job-name=extract_PSYMETAB_AurelieReymond                # Job name (this name will appear in the queue)
#SBATCH --chdir=/data/sgg2/jenny/projects/PSYMETAB                 # The Working Directory of the job
#SBATCH --ntasks=1                                                 # Run on a single core
#SBATCH --time=24:00:00                                            # Time limit d-hrs:min:sec
#SBATCH --output=data/processed/extractions/%x.out                 # Standard output and error log (%j: job allocation number)
#SBATCH --account=sgg                                              # runs on the sggg nodes.
#SBATCH --partition=sgg                                            # runs on the sggg nodes

project_dir="/data/sgg2/jenny/projects/PSYMETAB"
input=$projects/PSYMETAB/data/raw/extractions/AurelieReymond
output=$projects/PSYMETAB/data/processed/extractions/AurelieReymond
QC_dir=$project_dir/analysis/QC
#qc_data=${QC_dir}/15_final_processing/CEU/PSYMETAB_GWAS.CEU
#pc_data=${QC_dir}/15_final_processing/final_pca/CEU/pcs.PSYMETAB_GWAS_CEU_unrelated.txt
qc_data=${QC_dir}/15_final_processing/FULL/PSYMETAB_GWAS.FULL
reference_data=$data/HRC.r1-1.GRCh37.wgs.mac5.sites.tab
plink_data="PLINK_091019_0920"
input_chip=$project_dir/data/processed/${plink_data}/PSYMETAB_GWAS
pc_data=${QC_dir}/15_final_processing/final_pca/ETH/pcs.PSYMETAB_GWAS_ETH_unrelated.txt


if [ ! -d "${output}" ] ; then
mkdir ${output}
fi

FILES=$input/*.txt


for file in $FILES
do
  file_start=$(echo $file | cut -f1 -d.)
  out_name=$(basename $file_start)
  #if [ ! -d "${output}/${out_name}" ] ; then
    echo "Processing the following file: "$file
    echo $'\n'
    mkdir ${output}/${out_name}
    sh code/extractions/extract_snps.sh $file $output $qc_data $pc_data $QC_dir $input_chip
    echo "Finished processing file: "$file
    echo $'\n'
  #   echo $'\n'
  # else
  #   echo "The following SNP list has already been processed: "$file
  #   echo "If the code has been modified, delete the folder located here:"
  #   echo ${output}/${out_name}
  #   echo "and the code will rerun."
  #   echo $'\n'
  #   echo $'\n'
  # fi

done

#sbatch $projects/PSYMETAB/code/extractions/AurelieReymond_extract.sbatch
