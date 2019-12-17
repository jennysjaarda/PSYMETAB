
# Workflow Plans
#

qc_prep <- drake_plan (
  # Load and clean data to prepare for QC ----------------------
  qc_pheno_raw = read_excel(file_in(!!qc_pheno_file),sheet = 1),
  id_code = read.csv(file_in("data/raw/ID_key.csv"), header=T),
  rs_conversion = fread(file_in(!!rsconv_raw_file), data.table = F),
  fam_raw = read.table(file_in(!!paste0(plink_bed_out,".fam"))), # read the .fam file
  ped_folder = file_in(!!dirname(plink_ped_raw)),
  #bed_folder = dir.create(file_out(!!dirname(plink_bed_out)),showWarnings = F),
  #create_dir = dir.create(bed_folder, showWarnings = F),
  create_bed_out = {
    in_file <- paste0(ped_folder,"/",basename(!!plink_ped_raw))
    outfile <- paste0(bed_folder,"/",basename(!!plink_bed_out))
    run(command = "plink", c( "--file", in_file, "--make-bed", "--out", outfile), error_on_status = F)
    file_out(!!paste0(plink_bed_out,".fam"))},

    fam_munge = munge_fam(fam_raw, id_code),
    qc_pheno_munge = munge_qc_pheno(qc_pheno_raw,fam_munge),
    dups = find_dups(qc_pheno_munge, fam_munge),
    sex_info = qc_pheno_munge %>% dplyr::select(c(FID,IID,Sexe)),
    eth_info = qc_pheno_munge %>% dplyr::select(c(FID,IID,Ethnie)),
    sex_format = format_sex_file(sex_info),
    eth_format = format_eth_file(eth_info),
    rs_conversion_munge = rs_conversion %>% mutate(RsID = case_when(RsID=="." ~ Name, TRUE ~ RsID)),
    out_sex = write.table(sex_format,  file_out("data/processed/phenotype_data/PSYMETAB_GWAS_sex.txt"),row.names = F, quote = F, col.names = F),
    out_eth = write.table(eth_format,  file_out("data/processed/phenotype_data/PSYMETAB_GWAS_eth.txt"),row.names = F, quote = F, col.names = F),
    out_dups = write.table(dups$dups,file_out("data/processed/phenotype_data/PSYMETAB_GWAS_dupIDs.txt"),row.names = F, quote = F, col.names = F),
    out_dups_set = write.table(dups$dups_set, file_out("data/processed/phenotype_data/PSYMETAB_GWAS_dupIDs_set.txt"),row.names = F, quote = F, col.names = F),
    out_rs_conversion = write.table(rs_conversion_munge, file_out("data/processed/reference_files/rsid_conversion.txt"), row.names = F, quote = F, col.names = F)


)

vis_drake_graph(drake_config(qc_prep))

#vis_drake_graph(drake_config(qc_prep))

pre_impute_qc <- drake_plan(
  #### run pre-imupation quality control
  run_pre_imputation = target(
    {file_in("data/processed/phenotype_data/PSYMETAB_GWAS_sex.txt"); file_in("data/processed/phenotype_data/PSYMETAB_GWAS_eth.txt");
     file_in("data/processed/phenotype_data/PSYMETAB_GWAS_dupIDs.txt"); file_in("data/processed/phenotype_data/PSYMETAB_GWAS_dupIDs_set.txt");
     file_in("data/processed/reference_files/rsid_conversion.txt");processx::run(command = "sh", c( pre_imputation_script), error_on_status = F)},
     resources = list(core = 16)
  )
)

pre_impute <- bind_plans(qc_prep,pre_impute_qc)
vis_drake_graph(drake_config(pre_impute))

## here download output and upload to Michigan server
post_impute <- drake_plan(

  # run post-imputation quality control and processing ------------
  download_imputation = target(
    processx::run(command = "sh", c(file_in(!!download_imputation_script)), error_on_status = F),
    resources = list(ncpus = 16)),
  run_check_imputation = target({
      download_imputation
      processx::run(command = "sh", c( file_in(!!check_imputation_script)), error_on_status = F)},
    resources = list(ncpus = 16)),
  run_post_imputation = target({
      download_imputation
      processx::run(command = "sh", c( file_in(!!post_imputation_script)), error_on_status = F)},
    resources = list(ncpus = 16)),
  run_final_processing = target({
      run_post_imputation
      processx::run(command = "sh", c( file_in(!!final_processing_script)), error_on_status = F)},
    resources = list(ncpus = 16)),
  cp_qc_report = target({
      run_check_imputation
      processx::run(command = "cp", c( "analysis/QC/06_imputation_get/qcreport.html", "docs/generated_reports/"), error_on_status = F)},
    resources = list(ncpus = 1)),
  cp_qc_check = target({
    run_check_imputation
    processx::run("/bin/sh", c("-c","cp analysis/QC/07_imputation_check/summaryOutput/*html docs/generated_reports/"), error_on_status = F)},
    resources = list(ncpus = 1))
)

vis_drake_graph(drake_config(post_impute))

analysis_prep <- drake_plan(
  # prepare phenotype files for analysis in GWAS/GRS etc. ------------

  pheno_raw = readr::read_delim(file_in(!!pheno_file), col_types = cols(.default = col_character()), delim = ";") %>% type_convert(),
  bgen_sample_file = target({
    run_final_processing
    readr::read_delim(file_in("analysis/QC/15_final_processing/FULL/PSYMETAB_GWAS.FULL.sample"),
      col_types = cols(.default = col_character()), delim = " ") %>% type_convert()
    }),
  bgen_nosex_out = write.table(bgen_sample_file[,1:3],file_out("analysis/QC/15_final_processing/FULL/PSYMETAB_GWAS.FULL_nosex.sample"),row.names = F, quote = F, col.names = T),
  pc_raw = read_pcs(file_in(pc_dir)) %>% as_tibble(),
  pheno_munge = munge_pheno(pheno_raw), # pheno_munge %>% count(GEN) %>% filter(n!=1) ## NO DUPLICATES !
  pheno_baseline = inner_join(pc_raw %>% mutate_at("GPCR", as.character) , pheno_munge %>% mutate_at("GEN", as.character), by = c("GPCR" = "GEN")) %>%
    replace_na(list(sex='NONE')),

  pheno_followup = munge_pheno_follow(pheno_baseline), #names(pheno_followup) is the names defined in `test_drugs`: tibble
  GWAS_input = create_GWAS_pheno(pheno_baseline, pheno_followup),

  ## linear model GWAS:
  #linear_pheno = pheno_baseline %>% dplyr::select(matches('^FID$|^IID$|_start_Drug_1')),
  #linear_covar = pheno_baseline %>% dplyr::select(matches('^FID$|^IID$|^sex$|^Age_Drug_1$|Age_sq_Drug_1$|^PC[0-9]$|^PC[1][0-9]$|^PC20')),
  #pheno_followup_split = split_followup(pheno_followup),
  #followup_data_write = if(!is.null(create_pheno_folder)){ write_followup_data(pheno_followup_split)},
  out_linear_pheno =  target({
    write.table(GWAS_input$full_pheno, file_out("data/processed/phenotype_data/GWAS_input/pheno_input.txt"), row.names = F, quote = F, col.names = T)}),
  out_linear_covar =  target({
    write.table(GWAS_input$full_covar, file_out("data/processed/phenotype_data/GWAS_input/covar_input.txt"), row.names = F, quote = F, col.names = T)}),

)

vis_drake_graph(drake_config(analysis_prep))
make(analysis_prep)
#

init_analysis <- drake_plan(
  # run initial GWAS ---------------------------
  linear_out = target({
    run_final_processing
    run_gwas(pfile = ("analysis/QC/15_final_processing/FULL/PSYMETAB_GWAS.FULL"), pheno_file = file_in("data/processed/phenotype_data/GWAS_input/pheno_input.txt"),
               covar_file = file_in("data/processed/phenotype_data/GWAS_input/covar_input.txt"),
               threads = 16, pheno_name = pheno, covar_names = covars, eths = !!eths,
               eth_sample_file = "analysis/QC/12_ethnicity_admixture/pca/PSYMETAB_GWAS_ETH_samples.txt", ## this is not a real file - "ETH" gets replaced by proper "ETH" in `run_gwas`
               remove_sample_file = "analysis/QC/11_relatedness/PSYMETAB_GWAS_related_ids.txt",
               output_dir = file_in("analysis/GWAS"), output = "PSYMETAB_GWAS")},
    transform = map(.data = !!baseline_gwas_info, .id = pheno),resources = list(ncpus = 16)),

  interaction_out = target({
    run_final_processing
    run_gwas(pfile = "analysis/QC/15_final_processing/FULL/PSYMETAB_GWAS.FULL", pheno_file = file_in("data/processed/phenotype_data/GWAS_input/pheno_input.txt"),
                covar_file = file_in("data/processed/phenotype_data/GWAS_input/covar_input.txt"),
                threads = 16, pheno_name = pheno, covar_names = covars, parameters = parameters, interaction = TRUE, interaction_name = interaction_var_name,
                eths = !!eths, eth_sample_file = "analysis/QC/12_ethnicity_admixture/pca/PSYMETAB_GWAS_ETH_samples.txt",  ## this is not a real file - "ETH" gets replaced by proper "ETH" in `run_gwas`
                remove_sample_file = "analysis/QC/11_relatedness/PSYMETAB_GWAS_related_ids.txt",
                output_dir = file_out("analysis/GWAS"), output = "PSYMETAB_GWAS")},
    transform = map(.data = !!interaction_gwas_info, .id = interaction_var_name),resources = list(ncpus = 16)),
  subgroup_out = target({
    run_final_processing
    run_gwas(pfile = "analysis/QC/15_final_processing/FULL/PSYMETAB_GWAS.FULL", pheno_file = file_in("data/processed/phenotype_data/GWAS_input/pheno_input.txt"),
                covar_file = file_in("data/processed/phenotype_data/GWAS_input/covar_input.txt"),
                threads = 16, pheno_name = pheno, covar_names = covars, subgroup_var = subgroup, group = "subgroup",
                eths = !!eths, eth_sample_file = "analysis/QC/12_ethnicity_admixture/pca/PSYMETAB_GWAS_ETH_samples.txt",  ## this is not a real file - "ETH" gets replaced by proper "ETH" in `run_gwas`
                remove_sample_file = "analysis/QC/11_relatedness/PSYMETAB_GWAS_related_ids.txt",
                output_dir = file_out("analysis/GWAS"), output = "PSYMETAB_GWAS")},
    transform = map(.data = !!subgroup_gwas_info, .id = subgroup),resources = list(ncpus = 16))
)

post_impute <- bind_plans(post_impute,analysis_prep,init_analysis)

vis_drake_graph(drake_config(post_impute))
outdated(drake_config(init_analysis))

make(init_analysis,parallelism = "future",jobs = 4, console_log_file = "init.out")

make(init_analysis,parallelism = "clustermq",jobs = 4, console_log_file = "init_analysis.out", template = list(cpus = 16, partition = "sgg"))

pull(init_analysis %>% dplyr::select(command))

process_init <- drake_plan(

  gwas_figures = target(
    process_meta(outcome_variable,interaction_variable,model),
    transform = map(.data =!!!GWAS_models)
  )


)

post_impute <- bind_plans(qc_prep,pre_impute_qc, init_analysis, process_init)

make(post_impute,parallelism = "future",jobs = 2, console_log_file = "post_impute_qc.out", resources = list(partition = "sgg"))
make(init_analysis,parallelism = "clustermq",jobs = 4, console_log_file = "init_analysis.out", template = list(cpus = 16, partition = "cluster"))



### some useful code

c(init_analysis %>% dplyr::select(target))
outdated(drake_config(init_analysis))
make(init_analysis,parallelism = "clustermq",jobs = 4, console_log_file = "init_analysis.out", template = list(cpus = 16, partition = "cluster"))




#### run gwas and GRS computation
### to be grabble from `GWAS.sh` and `PRSice.sh`
)









####
visualize <- drake_plan(
  dim(dups) #30 2
  table(sex_info3$Sexe, exclude = NULL)
  #    F    M
  # 1298 1469
  table(eth_info3$Ethnie, exclude = NULL)
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
    no_sex_eth <- fam[which(!fam[,2] %in% sex_info4[,2]), c(1,2)]
    no_sex_eth
    # empty

    ### dimensions of PC data
    print(eth)
    print(dim(PC_temp))

    drug_list <- unique(unlist(full_pheno %>% dplyr::select(starts_with("AP1_Drug_"))))
    ### list of drugs included in phenofile



  create_report(pheno_baseline, config = configure_report(add_plot_prcomp = FALSE))


)
