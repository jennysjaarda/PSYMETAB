requestor='FirstnameLastname'
code_directory=$projects/PSYMETAB/code/extractions
master_file=$code_directory/master_extract.sbatch

sed -i "s/FirstnameLastname/$requestor/g" "$master_file" > $code_directory/${requestor}_extract_test.sbatch

sbatch $projects/PSYMETAB/code/extractions/${requestor}_extract.sbatch

## Copy to terminal:
#chmod +x  $code_directory/${requestor}_extract_test.sh
#$code_directory/${requestor}_extract_test.sh
