
# Workflow Plans
#

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
  #### run pre-imupation quality control
  run_pre_imputation = processx::run(command="sh",c( pre_imputation_script), error_on_status=F)
)

make(pre_impute_qc,parallelism = "clustermq",jobs = 1, console_log_file = "pre_impute_qc.out", template=list(cpus=16, partition="sgg"))

### download files and impute on Michigan server

post_impute <- drake_plan(
  #### run post-imputation quality control and processing
  download_imputation = processx::run(command="sh",c( download_imputation_script), error_on_status=F),
  run_check_imputation = if(!is.null(download_imputation)){processx::run(command="sh",c( check_imputation_script), error_on_status=F)},
  run_post_imputation = if(!is.null(download_imputation)){processx::run(command="sh",c( post_imputation_script), error_on_status=F)},
  run_final_processing = if(!is.null(run_post_imputation)){processx::run(command="sh",c( final_processing_script), error_on_status=F)},
  mk_reprt_folder = if(!is.null(run_check_imputation)){processx::run(command="mkdir",c( "docs/generated_reports"), error_on_status=F)},
  cp_qc_report = if(!is.null(mk_reprt_folder)){processx::run(command="cp",c( "analysis/QC/06_imputation_get/qcreport.html", "docs/generated_reports/"), error_on_status=F)},
  cp_qc_check = if(!is.null(mk_reprt_folder)){processx::run("/bin/sh", c("-c","cp analysis/QC/07_imputation_check/summaryOutput/*html docs/generated_reports/"), error_on_status=F)}
)

make(post_impute,parallelism = "clustermq",jobs = 2, console_log_file = "post_impute_qc.out", template=list(cpus=16, partition="sgg"))



analysis_prep <- drake_plan(
  #### prepare phenotype files for analysis in GWAS/GRS etc.
  create_pheno_folder = dir.create("data/processed/phenotype_data/GWAS_input",showWarnings=F),
  pheno_raw = readr::read_delim(file_in(pheno_file), col_types = cols(.default = col_character()), delim=";") %>% type_convert(),
  bgen_sample_file = readr::read_delim(file_in("analysis/QC/15_final_processing/FULL/PSYMETAB_GWAS.FULL.sample"), col_types = cols(.default = col_character()), delim=" ") %>% type_convert(),
  bgen_out = write.table(bgen_sample_file[,1:3],file_out("analysis/QC/15_final_processing/FULL/PSYMETAB_GWAS.FULL_nosex.sample"),row.names=F, quote=F, col.names=T),

  pc_raw = read_pcs(pc_dir) %>% as_tibble(),
  pheno_munge = munge_pheno(pheno_raw), # pheno_munge %>% count(GEN) %>% filter(n!=1) ## NO DUPLICATES !
  pheno_baseline = inner_join(pc_raw %>% mutate_at("GPCR", as.character) , pheno_munge %>% mutate_at("GEN", as.character), by = c("GPCR" = "GEN")) %>%
    replace_na(list(sex='NONE')),
  pheno_followup = munge_pheno_follow(pheno_baseline), #names(pheno_followup) is the names defined in `test_drugs`: tibble
  ## linear model GWAS:
  linear_pheno = pheno_baseline %>% dplyr::select(matches('^FID$|^IID$|_start_Drug_1')),
  linear_covar = pheno_baseline %>% dplyr::select(matches('^FID$|^IID$|^sex$|^Age_Drug_1$|Age_sq_Drug_1$|^PC[0-9]$|^PC[1][0-9]$|^PC20')),
  pheno_followup_split = split_followup(pheno_followup),
  followup_data_write = if(!is.null(create_pheno_folder)){ write_followup_data(pheno_followup_split)},

  out_linear_pheno =  if(!is.null(create_pheno_folder)){write.table(linear_pheno, file_out(!! paste0("data/processed/phenotype_data/GWAS_input/","linear_pheno_input.txt")), row.names=F, quote=F, col.names=T)},
  out_linear_covar =  if(!is.null(create_pheno_folder)){write.table(linear_covar, file_out(!! paste0("data/processed/phenotype_data/GWAS_input/","linear_covar_input.txt")), row.names=F, quote=F, col.names=T)},

  ## to be grabbed from `pheno_clean.r`


)

make(analysis_prep)
#vis_drake_graph(drake_config(analysis_prep))

init_analysis <- drake_plan(
  make_dirs_out = create_analysis_dirs(),
  gwas_interaction_track = file_in(!!gwas_interaction_script),
  gwas_linear_track = file_in(!!gwas_linear_script),

  run_gwas_interaction = target(
    if(make_dirs_out){processx::run(command="sh",c( gwas_interaction_track, variable), error_on_status=F)},
    transform = map(variable = !!unlist(test_drugs %>% dplyr::select(class)))
  ),
  run_gwas_linear = target(
    if(make_dirs_out){processx::run(command="sh",c( gwas_linear_track, variable), error_on_status=F)},
    transform = map(variable = !!paste0(baseline_vars,"_start_Drug_1"))
  ),
  run_meta_interaction = target(
    meta_interaction(run_gwas_interaction),
    transform = map(run_gwas_interaction)
  ),
  run_meta_linear = target(
    meta_linear(run_gwas_linear),
    transform = map(run_gwas_linear)
  ),
  meta_jobs = target(
    combine_targets(run_meta_linear,run_meta_interaction),
    transform = combine(run_meta_linear,run_meta_interaction)
  ),

  #run_grs = if(make_dirs_out){processx::run(command="sh",c( compute_grs_script), error_on_status=F)}
)

c(init_analysis %>% dplyr::select(target))
outdated(drake_config(init_analysis))

make(init_analysis,parallelism = "clustermq",jobs = 4, console_log_file = "init_analysis.out", template=list(cpus=16, partition="cluster"))

process_init <- plan(

  gwas_figures = target(
    process_meta(outcome_variable,interaction_variable,model),
    transform = map(.data =!!!GWAS_models)
  )


)
#### run gwas and GRS computation
### to be grabble from `GWAS.sh` and `PRSice.sh`
)

####
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

    ### dimensions of PC data
    print(eth)
    print(dim(PC_temp))

    drug_list <- unique(unlist(full_pheno %>% dplyr::select(starts_with("AP1_Drug_"))))
    ### list of drugs included in phenofile



  create_report(pheno_baseline, config = configure_report(add_plot_prcomp = FALSE))


)
