#!/bin/bash
# redirect stdout/stderr to a file
exec &> analysis/QC/final_processing.log

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
#------------------------------------------------ FINAL CHECK ----------------------------------------------#
#############################################################################################################

if [ ! -d "${output_dir}/15_final_processing" ] ; then
  mkdir ${output_dir}/15_final_processing
fi
cd ${output_dir}/15_final_processing

## only perform for samples > 100
## create final psam file
for eth in CEU EA MIXED NA YRI ; do
  count=$(wc -l < "../12_ethnicity_admixture/pca/${output_name}_${eth}_samples.txt")

  if [ "$count" -ge 100 ] ; then
    if [ ! -d "${eth}" ] ; then
      mkdir ${eth}
    fi

  plink2 --pfile ../14_mafcheck/${output_name}.mafcheck.step14 \
    --make-pfile \
    --keep ../12_ethnicity_admixture/pca/${output_name}_${eth}_samples.txt \
    --out $eth/${output_name}.${eth} \
    --freq \
    --threads 16

  plink2 --pfile $eth/${output_name}.${eth} \
    --export bgen-1.2 id-paste=iid \
    --out $eth/${output_name}.${eth} \
    --threads 16

  plink2 --pfile $eth/${output_name}.${eth} \
    --recode vcf id-paste=iid \
    --out $eth/${output_name}.${eth} \
    --threads 16

  fi

done

if [ ! -d "FULL" ] ; then
  mkdir FULL
fi

plink2 --pfile ../14_mafcheck/${output_name}.mafcheck.step14 \
  --make-pfile \
  --out FULL/${output_name}.FULL \
  --threads 16 \
  --freq

plink2 --pfile FULL/${output_name}.FULL \
  --export bgen-1.2 id-paste=iid \
  --out FULL/${output_name}.FULL \
  --threads 16

plink2 --pfile FULL/${output_name}.FULL \
  --recode vcf id-paste=iid \
  --out FULL/${output_name}.FULL \
  --threads 16

bcftools query -f '%CHROM %POS %ID %REF %ALT %INFO/R2{0}\n' FULL/${output_name}.FULL.vcf > ${output_name}.info
sed  -i '1i CHROM POS ID REF ALT R2' ${output_name}.info

if [ ! -d "final_pca" ] ; then
  mkdir final_pca
fi

cd final_pca
### subset into seperate ethnicity files for final PCA calculation
for eth in CEU EA MIXED NA YRI
do
  if [ ! -d "${eth}" ] ; then
    mkdir ${eth}
  fi

  ## extract ethnicity samples from typed file with duplicates removed
  plink2 --bfile ../../12_ethnicity_admixture/${output_name}.typed_remove_dups \
    --keep ../../12_ethnicity_admixture/pca/${output_name}_${eth}_samples.txt \
    --make-bed \
    --threads 16 \
    --out ${eth}/${output_name}.${eth}

  ## re-clean data within each ethnicity before final PCA
  plink2 --bfile ${eth}/${output_name}.${eth} \
    --geno 0.1 \
    --maf 0.05 \
    --hwe 5e-4 \
    --out ${eth}/${output_name}.${eth}.QC \
    --make-bed

  ## prune, excluding related individuals
  plink2 --bfile ${eth}/${output_name}.${eth}.QC \
    --remove ../../11_relatedness/${output_name}_related_ids.txt \
    --indep-pairwise 50 5 0.2 \
    --out ${eth}/${output_name}.${eth}.indep \
    --threads 16

  ## create an unrelated set for calculating PCs with same set of pruned SNPs
  plink2 --bfile ${eth}/${output_name}.${eth}.QC \
    --extract ${eth}/${output_name}.${eth}.indep.prune.in \
    --remove ../../11_relatedness/${output_name}_related_ids.txt \
    --make-bed \
    --threads 16 \
    --out ${eth}/${output_name}.${eth}.pruned.unrelated

  ## create a full set based on pruning in unrelated samples
  plink2 --bfile ${eth}/${output_name}.${eth} \
    --extract ${eth}/${output_name}.${eth}.indep.prune.in \
    --make-bed \
    --threads 16 \
    --out ${eth}/${output_name}.${eth}.pruned.full

  cd ${eth}
  ## compute PCA components
  $flashpca/flashpca --bfile ${output_name}.${eth}.pruned.unrelated \
    --suffix .${output_name}_${eth}_unrelated.txt \
    --ndim 100 \
    --outload ${eth}_loadings.txt \
    --outmeansd ${eth}_meansd.txt \
    --numthreads 16

  ## project onto full set
  $flashpca/flashpca --bfile ${output_name}.${eth}.pruned.full \
    --project \
    --inmeansd ${eth}_meansd.txt \
    --outproj ${output_name}_${eth}_projections.txt \
    --inload ${eth}_loadings.txt -v \
    --numthreads 16

  cd ${output_dir}/15_final_processing/final_pca
done
