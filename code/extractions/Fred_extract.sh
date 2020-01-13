#!/bin/bash

#SBATCH --job-name=extract_PSYMET_Fred                             # Job name (this name will appear in the queue)
#SBATCH --workdir=/data/sgg2/jenny/projects/PSYMETAB          # The Working Directory of the job
#SBATCH --ntasks=1                                                 # Run on a single core
#SBATCH --time=24:00:00                                            # Time limit d-hrs:min:sec
#SBATCH --output=data/processed/extractions/%x.out                 # Standard output and error log (%j: job allocation number)
#SBATCH --account=sgg                                              # runs on the sggg nodes.
#SBATCH --partition=sgg                                            # runs on the sggg nodes

project_dir="/data/sgg2/jenny/projects/PSYMETAB"
input=$projects/PSYMETAB/data/raw/extractions/Fred
output=$projects/PSYMETAB/data/processed/extractions/Fred
QC_dir=$project_dir/analysis/QC
qc_data=${QC_dir}/15_final_processing/CEU/PSYMETAB_GWAS.CEU
pc_data=${QC_dir}/15_final_processing/final_pca/CEU/pcs.PSYMETAB_GWAS_CEU_unrelated.txt
reference_data=$data/HRC.r1-1.GRCh37.wgs.mac5.sites.tab

if [ ! -d "${output}" ] ; then
mkdir ${output}
fi

FILES=$input/*.txt

for file in $FILES
do
  sh extract_snps.sh $file $output $qc_data $pc_data

done

#sbatch $projects/PSYMETAB/code/extractions/Fred_extract.sh
