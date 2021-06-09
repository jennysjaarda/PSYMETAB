#!/bin/bash

#SBATCH --job-name=extract_PSYMETAB_master                         # Job name (this name will appear in the queue)
#SBATCH --chdir=/data/sgg2/jenny/projects/PSYMETAB                 # The Working Directory of the job
#SBATCH --ntasks=1                                                 # Run on a single core
#SBATCH --time=24:00:00                                            # Time limit d-hrs:min:sec
#SBATCH --output=data/processed/extractions/%x.out                 # Standard output and error log (%j: job allocation number)
#SBATCH --account=sgg                                              # runs on the sggg nodes.
#SBATCH --partition=cluster2                                            # runs on the sggg nodes

code_directory=$projects/PSYMETAB/code/extractions
template_file=$code_directory/template_extract.sh

requestor_folders=`ls -d data/raw/extractions/*`

for requestor_folder in $requestor_folders

do

  requestor=$(basename $requestor_folder)
  echo "Processing the SNP extraction lists for: "$requestor

  sed "s/FirstnameLastname/$requestor/g" "$template_file" > $code_directory/${requestor}_extract.sh

  sbatch $projects/PSYMETAB/code/extractions/${requestor}_extract.sh

  wait

done


#sbatch $projects/PSYMETAB/code/extractions/master_extract.sh
