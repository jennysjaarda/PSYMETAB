#!/bin/bash
#SBATCH --partition sgg
#SBATCH --workdir /data/sgg2/jenny/projects/PSYMETAB
#SBATCH --job-name master
#SBATCH --output master.out
#SBATCH --mail-type=ALL                                            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jennysjaarda@gmail.com                         # Where to send mail
#SBATCH --ntasks=1                                                 # Run on a single core


Rscript make.R
