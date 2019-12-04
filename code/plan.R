
# Workflow Plans
#
# ## Data sources
#
# The raw data are in `.xlsx` files. Some of the data are pre-processed
# (mean and sd of count and length of gemmae, 30 min averages of PPFD).
#
# The files include some additional analsyes and other data that won't
# be used here.


qc_prep <- drake_plan (

  # Load and clean data ----

  qc_pheno_raw = read_excel(file_in(qc_pheno_file),sheet = 1),
  id_code = read.csv("data/raw/ID_key.csv", header=T),
  rs_conversion = fread(rsconv_raw_file, data.table=F),
  fam_raw = read_fam(create_bed_out), # read the .fam file
  create_dir = dir.create(dirname(plink_bed_out), showWarnings=F),
  create_bed_out = if(!is.null(create_dir)){run(command="plink",c( "--file", plink_ped_raw, "--make-bed", "--out", plink_bed_out), error_on_status=F)},
  fam_munge = munge_fam(fam_raw, id_code),
  qc_pheno_munge = munge_qc_pheno(qc_pheno_raw,fam_munge),
  dups = find_dups(qc_pheno_munge, fam_munge),
  sex_info = qc_pheno_munge %>% dplyr::select(c(FID,IID,Sexe)),
  eth_info = qc_pheno_munge %>% dplyr::select(c(FID,IID,Ethnie)),
  sex_format = format_sex_file(sex_info),
  eth_format = format_eth_file(eth_info),
  rs_conversion_munge = rs_conversion %>% mutate(RsID = case_when(RsID=="." ~ Name, TRUE ~ RsID)),
  out_sex = write.table(sex_info,  file_out("data/processed/phenotype_data/PSYMETAB_GWAS_sex.txt"),row.names=F, quote=F, col.names=F),
  out_eth = write.table(sex_info,  file_out("data/processed/phenotype_data/PSYMETAB_GWAS_eth.txt"),row.names=F, quote=F, col.names=F),
  out_dups = write.table(dups$dups,file_out("data/processed/phenotype_data/PSYMETAB_GWAS_dupIDs.txt"),row.names=F, quote=F, col.names=F),
  out_dups_set = write.table(dups$dups_set, file_out("data/processed/phenotype_data/PSYMETAB_GWAS_dupIDs_set.txt"),row.names=F, quote=F, col.names=F),
  out_rs_conversion = write.table(rs_conversion_munge, file_out("data/processed/reference_files/rsid_conversion.txt"), row.names=F, quote=F, col.names=F)

)

make(qc_prep)

pre_impute_qc <- drake_plan(
  run_pre_imputation = processx::run(command="sh",c( pre_imputation_script), error_on_status=F)
)

make(pre_impute_qc,parallelism = "clustermq",jobs = 1, console_log_file = "pre_impute_qc.out", template=list(cpus=16, partition="sgg"))

### download files and impute on Michigan server

post_impute <- drake_plan(
  #download_imputation = processx::run(command="sh",c( download_imputation_script), error_on_status=F),
  #run_check_imputation = if(!is.null(download_imputation)){processx::run(command="sh",c( check_imputation_script), error_on_status=F)}
  run_post_imputation = if(!is.null(download_imputation)){processx::run(command="sh",c( post_imputation_script), error_on_status=F)}
)

make(post_impute,parallelism = "clustermq",jobs = 1, console_log_file = "post_impute_qc.out", template=list(cpus=16, partition="sgg"))



visualize <- drake_plan(
  dim(dups) #30 2
  table(sex_info3$Sexe, exclude=NULL)
  #    F    M
  # 1298 1469
  table(eth_info3$Ethnie, exclude=NULL)
  # africain  africain + caucasien       amerique du sud
  #      101                     1                     1
  # Antilles                 arabe     arabe + caucasien
  #        1                    60                     3
  # asiatique asiatique + caucasien                 autre
  #       25                     1                   178
  # caucasien               inconnu                  <NA>
  #     1667                   458                   271

  ## double check no duplicates
  sex_info4 <- unique(sex_info)
  dim(sex_info)==dim(sex_info4)
  sex_info4[which(duplicated(sex_info4[,2])),]
  duplicate_IDs <- sex_info4[which(duplicated(sex_info4[,2])),2]
  sex_info4[which(sex_info4[,2]==duplicate_IDs),]
  #TRUE

  eth_info4 <- unique(eth_info)
  dim(eth_info3)==dim(eth_info4)
  #TRUE
  eth_info4[which(duplicated(eth_info4[,2])),]
  #NONE

    ## anyone missign?
    no_sex_eth <- fam[which(!fam[,2] %in% sex_info4[,2]),c(1,2)]
    no_sex_eth
    # empty
)
