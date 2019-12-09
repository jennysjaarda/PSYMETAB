

project_dir="/data/sgg2/jenny/projects/PSYMETAB_GWAS"
plink_data="PLINK_091019_0920"
raw_data=$project_dir/data/raw/$plink_data
QC_dir=$project_dir/analysis/QC
GWAS_dir=$project_dir/analysis/GWAS
input_chip=$project_dir/data/processed/${plink_data}/PSYMETAB
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

pheno_dir=$project_dir/data/processed/phenotype_data/GWAS_input
variable=$1

for pc in $(seq 1 22);do printf PC${pc}','; done

sex, Age_Drug_1




cd $GWAS_dir

for eth in CEU EA MIXED NA YRI ; do
  count=$(wc -l < "$QC_dir/12_ethnicity_admixture/pca/${output_name}_${eth}_samples.txt")

  if [ "$count" -ge 100 ] ; then
    if [ ! -d "${eth}" ] ; then
      mkdir ${eth}
    fi

### interaction model

	plink2 --pfile $QC_dir/15_final_processing/${eth}/${output_name}.${eth} \
		--pheno $interaction_pheno_file \
		--linear interaction \
		--pheno-name bmi_change \
		--keep $QC_dir/12_ethnicity_admixture/pca/${output_name}_${eth}_samples.txt \
		--remove $QC_dir/11_relatedness/${output_name}_related_ids.txt \
		--covar $interaction_covar_file \
		--threads 16 \
		--out $GWAS_dir/${eth}/PSYMETAB_GWAS_interaction \
		--parameters 1-27, 49 \
    --covar-variance-standardize

	fi

done
