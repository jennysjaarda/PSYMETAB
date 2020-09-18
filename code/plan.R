

# Workflow Plans
#

qc_prep <- drake_plan (
  # Load and clean data to prepare for QC ----------------------
  qc_pheno_raw = read_excel(file_in(!!qc_pheno_file),sheet = 1),
  id_code = read.csv(file_in("data/raw/ID_key.csv"), header=T),
  rs_conversion = fread(file_in(!!rsconv_raw_file), data.table = F),
  fam_raw = read.table(file_in(!!paste0(plink_bed_out,".fam"))), # read the .fam file
  ped_folder = file_in(!!dirname(plink_ped_raw)),
  #bed_folder = file_out(!!dirname(plink_bed_out)),
  #create_dir = dir.create(bed_folder, showWarnings = F),
  create_bed_out = {
    in_file <- paste0(ped_folder,"/",basename(!!plink_ped_raw))
    outfile <- paste0(!!dirname(plink_bed_out),"/",basename(!!plink_bed_out))
    processx::run(command = "plink", c( "--file", in_file, "--make-bed", "--out", outfile), error_on_status = F)
    file_out(!!paste0(plink_bed_out,".fam"))
    file_out(!!dirname(plink_bed_out))},

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
  out_rs_conversion = {write.table(rs_conversion_munge, file_out("data/processed/reference_files/rsid_conversion.txt"), row.names = F, quote = F, col.names = F)},

)

psy_miscellaneous <- drake_plan(
  custom_snp_list = read_excel(file_in("data/raw/custom_SNP_list.xlsx"),sheet = 1),
  HRC_list = fread("/data/sgg2/jenny/data/HRC.r1-1.GRCh37.wgs.mac5.sites.tab", data.table=F),
  ## check which SNPs in the custom SNP list (defined by Aurelie and others in the group at Cery) would already be imputed using the HRC panel.
  ## if SNPs are already in the HRC list above, then there is no need to genotype them directly.

  custom_snp_list_non_hrc = custom_snp_list %>% filter(!SNP %in%  HRC_list$ID),
  out_non_hrc_snps = write.table(custom_snp_list_non_hrc, file_out("data/processed/non_hrc_custom_snps.txt"),row.names = F, quote = F, col.names = T),


)
# config <- drake_config(qc_prep)
# vis_drake_graph(config)

pre_impute_qc <- drake_plan(
  # run pre-imupation quality control ----------------------
  run_pre_imputation = target(
    {file_in("data/processed/phenotype_data/PSYMETAB_GWAS_sex.txt"); file_in("data/processed/phenotype_data/PSYMETAB_GWAS_eth.txt");
     file_in("data/processed/phenotype_data/PSYMETAB_GWAS_dupIDs.txt"); file_in("data/processed/phenotype_data/PSYMETAB_GWAS_dupIDs_set.txt");
     file_in("data/processed/reference_files/rsid_conversion.txt"); processx::run(command = "sh", c( pre_imputation_script), error_on_status = F)
     file_out("analysis/QC/00_preprocessing", "analysis/QC/01_strandflip", "analysis/QC/02_maf_zero", "analysis/QC/03_missingness", "analysis/QC/04_sexcheck")}
  ),
)

# config <- drake_config(pre_impute_qc)
# vis_drake_graph(config)

################################################################
## here download output and upload to Michigan server
################################################################

download_impute <- drake_plan(

  # download imputation ------------
  download_imputation = target({
      processx::run(command = "sh", c(file_in(!!download_imputation_script)), error_on_status = F)
  }),
)

post_impute <- drake_plan(

  # run post-imputation quality control and processing ------------
  run_check_imputation = target({
      # file_in("analysis/QC/06_imputation_get")
      # this file will get zipped later so we cannot track it
      processx::run(command = "sh", c( file_in(!!check_imputation_script)), error_on_status = F)
      file_out("analysis/QC/07_imputation_check")
    }),
  run_post_imputation = target({
      file_in("code/qc/ethnicity_check.R", "code/qc/relatedness_filter.R", "code/qc/maf_check.R", "code/qc/update_pvar.R")
      processx::run(command = "sh", c( file_in(!!post_imputation_script)), error_on_status = F)
      # file_out("analysis/QC/08_plink_convert"), "analysis/QC/09_extract_typed", "analysis/QC/10_merge_imputed", "analysis/QC/11_relatedness",
      #         "analysis/QC/12_ethnicity_admixture", "analysis/QC/13_hwecheck", "analysis/QC/14_mafcheck")
      # these files are too big to track
    }),
  run_final_processing = target({
      run_post_imputation
      file_in(!!(paste0("analysis/QC/14_mafcheck/", study_name, ".mafcheck.step14.log")))
      # file_in("analysis/QC/11_relatedness", "analysis/QC/12_ethnicity_admixture", "analysis/QC/14_mafcheck")
      processx::run(command = "sh", c( file_in(!!final_processing_script)), error_on_status = F)
      # file_out("analysis/QC/15_final_processing")
      # these files are too big to track
    }),
  cp_qc_report = target({
      file_in("analysis/QC/06_imputation_get/qcreport.html")
      processx::run(command = "cp", c( "analysis/QC/06_imputation_get/qcreport.html", "docs/generated_reports/"), error_on_status = F)
      file_out("docs/generated_reports/qcreport.html")
    }, hpc = FALSE),
  cp_qc_check = target({
      file_in("analysis/QC/07_imputation_check")
      processx::run("/bin/sh", c("-c","cp analysis/QC/07_imputation_check/summaryOutput/*html docs/generated_reports/"), error_on_status = F)
      file_out("docs/generated_reports/07_imputation_check.html", "docs/generated_reports/CrossCohortReport.html")
  }, hpc = FALSE),
  snp_weights = target({
    run_post_imputation
    read.table((file_in(!!paste0("analysis/QC/12_ethnicity_admixture/snpweights/", study_name, ".NA.predpc"))))
  }, hpc = FALSE),
  related_inds = target({
    run_post_imputation
    read.table((file_in(!!paste0("analysis/QC/11_relatedness/", study_name, "_related_ids.txt"))))
  }, hpc = FALSE),
  snp_weights_munge = target(munge_snp_weights(snp_weights, related_inds), hpc = FALSE),
  snp_weights_out = target(write.table(snp_weights_munge, file_out(!!paste0("data/processed/extractions/", study_name, ".snpweights"))), hpc = FALSE)

)

qc_process <- drake_plan(

  imputed_var = count_imputed_var(),

  missingness_variants_removed = {countLines(paste0("analysis/QC/02_maf_zero/", !!study_name, ".maf_zero.step2.bim"))[1] -
    countLines(paste0("analysis/QC/03_missingness/", !!study_name, ".missingness_geno_10.bim"))[1] +
    countLines(paste0("analysis/QC/03_missingness/", !!study_name, ".missingness_mind_10.bim"))[1] -
    countLines(paste0("analysis/QC/03_missingness/", !!study_name, ".missingness_geno_05.bim"))[1] +
    countLines(paste0("analysis/QC/03_missingness/", !!study_name, ".missingness_mind_05.bim"))[1] -
    countLines(paste0("analysis/QC/03_missingness/", !!study_name, ".missingness_geno_01.bim"))[1]},

   missingness_indiv_removed = {countLines(paste0("analysis/QC/02_maf_zero/", !!study_name, ".maf_zero.step2.fam"))[1] -
    countLines(paste0("analysis/QC/03_missingness/", !!study_name, ".missingness_mind_10.fam"))[1] +
    countLines(paste0("analysis/QC/03_missingness/", !!study_name, ".missingness_geno_05.fam"))[1] -
    countLines(paste0("analysis/QC/03_missingness/", !!study_name, ".missingness_mind_05.fam"))[1] +
    countLines(paste0("analysis/QC/03_missingness/", !!study_name, ".missingness_geno_01.fam"))[1] -
    countLines(paste0("analysis/QC/03_missingness/", !!study_name, ".missingness.step3.fam"))[1]},

  initial_variants = countLines(paste0("analysis/QC/02_maf_zero/", !!study_name, ".maf_zero.step2.bim"))[1],
  initial_ind = countLines(paste0("analysis/QC/02_maf_zero/", !!study_name, ".maf_zero.step2.fam"))[1],

  preprocess_var_init = countLines(paste0(!!plink_bed_out, ".bim")),
  preprocess_ind_init = countLines(paste0(!!plink_bed_out, ".bim")),

  preprocess_non_x = countLines(paste0("analysis/QC/00_preprocessing/non_autosome_or_x_variants")),
  preprocess_non_x_remain = countLines(paste0("analysis/QC/00_preprocessing/temp1.bim")),

  preprocess_var_dups = countLines(paste0("analysis/QC/00_preprocessing/duplicate_variants")),
  preprocess_var_dups_remain = countLines(paste0("analysis/QC/00_preprocessing/temp2.bim")),


)

# config <- drake_config(post_impute)
# vis_drake_graph(config)

analysis_prep <- drake_plan(
  # prepare phenotype files for analysis in GWAS/GRS etc. ------------

  pheno_raw = readr::read_delim(file_in(!!pheno_file), na = c("#pdt", "#ccannul", "#pam", "#di", "#pnc"), col_types = cols(.default = col_character()), delim = ",") %>% type_convert(),
  caffeine_raw = target(read_excel(file_in(!!caffeine_file), sheet=1) %>% type_convert(),
    hpc = FALSE),
  # found a mistake that the age of WIFRYNSK was wrong at one instance.

  caffeine_munge = munge_caffeine(caffeine_raw),
  bgen_sample_file = target({
    readr::read_delim(file_in(!!paste0("analysis/QC/15_final_processing/FULL/", study_name, ".FULL.sample")),
      col_types = cols(.default = col_character()), delim = " ") %>% type_convert()
    }),
  bgen_nosex_out = write.table(bgen_sample_file[,1:3],file_out(!!paste0("analysis/QC/15_final_processing/FULL/", study_name, ".FULL_nosex.sample")),row.names = F, quote = F, col.names = T),
  pc_raw = read_pcs(file_in(!!pc_dir), !!study_name, !!eths) %>% as_tibble(),

  pheno_munge = munge_pheno(pheno_raw, !!baseline_vars, !!leeway_time, caffeine_munge, !!follow_up_limit, !!interaction_outcome), # pheno_munge %>% count(GEN) %>% filter(n!=1) ## NO DUPLICATES !
  #pheno_munge = munge_pheno(pheno_raw, baseline_vars, leeway_time, caffeine_munge, follow_up_limit)
  pheno_merge = merge_pheno_caffeine(pheno_raw, caffeine_munge, !!anonymization_error), # in the end, we don't use this dataset, but if you want to merge by caffeine with the appropriate date use this data.

  pheno_baseline = inner_join(pc_raw %>% mutate_at("GPCR", as.character), pheno_munge %>% mutate_at("GEN", as.character), by = c("GPCR" = "GEN")) %>%
    replace_na(list(sex='NONE')),
  pheno_eths_out = write.table(pheno_baseline %>% tidyr::separate(FID, c("COUNT", "GPCR"), "_") %>%
    dplyr::select(COUNT,GPCR,eth ), file_out(!!paste0("data/processed/phenotype_data/", study_name, "_inferred_eths.txt")), row.names = F, quote = F, col.names = T),

  pc_raw_out = write.table(pc_raw, file_out(!!paste0("data/processed/phenotype_data/", study_name, "_PCs.txt")), row.names = F, quote = F, col.names = T),

  pheno_drug_combinations = define_pheno_drug_combos(!!drug_classes, !!interaction_vars),
  pheno_followup = target(munge_pheno_follow(pheno_baseline, !!test_drugs, pheno_drug_combinations$class, pheno_drug_combinations$pheno, !!low_inducers, !!high_inducers),
    dynamic = map(pheno_drug_combinations)), #names(pheno_followup) is the names defined in `test_drugs`: tibble
  # follow-up data is smaller than baseline data because participants were filtered if who have they had taken drug but it wasn't followed (because we don't know when they would have taken the drug)
  GWAS_input = create_GWAS_pheno(pheno_baseline, pheno_followup, !!caffeine_vars, !!test_drugs,
    !!high_inducers, !!med_inducers, !!low_inducers, !!outcomes),

  ## linear model GWAS:
  #linear_pheno = pheno_baseline %>% dplyr::select(matches('^FID$|^IID$|_start_Drug_1')),
  #linear_covar = pheno_baseline %>% dplyr::select(matches('^FID$|^IID$|^sex$|^Age_Drug_1$|Age_sq_Drug_1$|^PC[0-9]$|^PC[1][0-9]$|^PC20')),
  #pheno_followup_split = split_followup(pheno_followup, !!test_drugs),
  #followup_data_write = if(!is.null(create_pheno_folder)){ write_followup_data(pheno_followup_split)},
  out_pheno_munge =  target({
    write.table(pheno_munge, file_out("data/processed/phenotype_data/GWAS_input/pheno_munge.txt"), row.names = F, quote = F, col.names = T)}),

  out_linear_pheno =  target({
    GWAS_input
    write.table(GWAS_input$full_pheno, file_out("data/processed/phenotype_data/GWAS_input/pheno_input.txt"), row.names = F, quote = F, col.names = T)}),
  out_linear_covar =  target({
    GWAS_input
    write.table(GWAS_input$full_covar, file_out("data/processed/phenotype_data/GWAS_input/covar_input.txt"), row.names = F, quote = F, col.names = T)}),

  baseline_gwas_info = target(define_baseline_models(GWAS_input, !!baseline_vars, !!drug_classes, !!caffeine_vars, !!interaction_outcome),
    hpc = FALSE),
  interaction_gwas_info = target(define_interaction_models(GWAS_input, !!drug_classes, !!interaction_outcome, !!baseline_vars),
    hpc = FALSE),
  subgroup_gwas_info = target(define_subgroup_models(GWAS_input, !!drug_classes, !!interaction_outcome, !!baseline_vars),
    hpc = FALSE),

  baseline_resid = target(
    residualize_pheno(GWAS_input, pheno = baseline_gwas_info$pheno, covars = baseline_gwas_info$covars, outcome_type = "baseline", !!eths, eth_data = pc_raw),
    dynamic = map(baseline_gwas_info)
  ),

  baseline_resid_data = reduce(baseline_resid, inner_join),

  interaction_resid = target(
    residualize_pheno(GWAS_input, pheno = interaction_gwas_info$pheno, covars = interaction_gwas_info$covars, outcome_type = "interaction", !!eths, eth_data = pc_raw),
    dynamic = map(interaction_gwas_info)
  ),

  interaction_resid_data = reduce(interaction_resid, inner_join),

  subgroup_resid = target(
    residualize_pheno(GWAS_input, pheno = subgroup_gwas_info$pheno, covars = subgroup_gwas_info$covars, outcome_type = "subgroup", !!eths, eth_data = pc_raw),
    dynamic = map(subgroup_gwas_info)
  ),
  subgroup_resid_data = reduce(subgroup_resid, inner_join),

  drug_classes_tibble = tibble(class = !!drug_classes),

  out_baseline_resid =  target({
    write.table(baseline_resid_data, file_out("data/processed/phenotype_data/GWAS_input/baseline_input_resid.txt"), row.names = F, quote = F, col.names = T)}),

  out_interacation_resid = target(
    write_interaction_resid(interaction_resid_data, drug_classes_tibble$class),
    dynamic=map(drug_classes_tibble), format = "file"),

  out_subgroup_resid = target(
    write_subgroup_resid(subgroup_resid_data, drug_classes_tibble$class),
    dynamic=map(drug_classes_tibble), format = "file"),

)

# config <- drake_config(analysis_prep)
# vis_drake_graph(config)

##### Run GWAS in PSYMETAB ---------

init_analysis <- drake_plan(

  GWAS_input_analysis = target(list(full_pheno = read_table2(file_in("data/processed/phenotype_data/GWAS_input/pheno_input.txt")),
                                    full_covar = read_table2(file_in("data/processed/phenotype_data/GWAS_input/covar_input.txt"))),
    hpc = FALSE), # this should be identical to GWAS_input above, but needs to be a different name

  # define endpoint, covars and outputs ---------------------------

  baseline_gwas_input = target(define_baseline_models(GWAS_input_analysis, !!baseline_vars, !!drug_classes, !!caffeine_vars, !!interaction_outcome),
    hpc = FALSE),
  interaction_gwas_input = target(define_interaction_inputs(!!drug_classes),
    hpc = FALSE),
  subgroup_gwas_input = target(define_subgroup_inputs(!!drug_classes),
    hpc = FALSE),


  # run initial GWAS ---------------------------

  linear_out = target({
    #loadd(baseline_gwas_info
    file_in(!!paste0("analysis/QC/15_final_processing/FULL/", study_name, ".FULL.log"))
    run_gwas(pfile = paste0("analysis/QC/15_final_processing/FULL/", !!study_name, ".FULL"),
               pheno_file = file_in("data/processed/phenotype_data/GWAS_input/baseline_input_resid.txt"),
               threads = 8, eths = !!eths,
               eth_sample_file = paste0("analysis/QC/12_ethnicity_admixture/pca/", !!study_name, "_ETH_samples.txt"), ## this is not a real file - "ETH" gets replaced by proper "ETH" in `run_gwas`
               eth_low_maf_file = paste0("analysis/QC/14_mafcheck/", !!study_name, "_ETH_low_maf_snps.txt"), ## this is not a real file - "ETH" gets replaced by proper "ETH" in `run_gwas`
               remove_sample_file = file_in(!!paste0("analysis/QC/11_relatedness/", study_name, "_related_ids.txt")),
               output_dir = (!!plink_output_dir), output = !!study_name)}
  ),

  check_interaction_input_files = target({
    c(interaction_gwas_input$pheno_file, interaction_gwas_input$covar_file)},
    dynamic=map(interaction_gwas_input), format = "file"),

  check_subgroup_input_files = target({
    subgroup_gwas_input$pheno_file},
    dynamic=map(subgroup_gwas_input), format = "file"),


  interaction_out = target({
    check_interaction_input_files
    file_in(!!paste0("analysis/QC/15_final_processing/FULL/", study_name, ".FULL.log"))
    run_gwas(pfile = paste0("analysis/QC/15_final_processing/FULL/", !!study_name, ".FULL"), pheno_file = interaction_gwas_input$pheno_file,
                covar_file = interaction_gwas_input$covar_file,type = "interaction",
                threads = 8, covar_names = interaction_gwas_input$covars, parameters = interaction_gwas_input$parameters, output_suffix = interaction_gwas_input$output_suffix,
                eths = !!eths, eth_sample_file = paste0("analysis/QC/12_ethnicity_admixture/pca/", !!study_name, "_ETH_samples.txt"),  ## this is not a real file - "ETH" gets replaced by proper "ETH" in `run_gwas`
                eth_low_maf_file = paste0("analysis/QC/14_mafcheck/", !!study_name, "_ETH_low_maf_snps.txt"), ## this is not a real file - "ETH" gets replaced by proper "ETH" in `run_gwas`
                remove_sample_file = file_in(!!paste0("analysis/QC/11_relatedness/", study_name, "_related_ids.txt")),
                output_dir = (!!plink_output_dir), output = !!study_name)},
    dynamic = map(interaction_gwas_input)
  ),

  subgroup_out = target({
    check_subgroup_input_files
    file_in(!!paste0("analysis/QC/15_final_processing/FULL/", study_name, ".FULL.log"))
    run_gwas(pfile = paste0("analysis/QC/15_final_processing/FULL/", !!study_name, ".FULL"), pheno_file = subgroup_gwas_input$pheno_file,
                type = "subgroup", subgroup_var = subgroup_gwas_input$subgroup,
                threads = 8, output_suffix = subgroup_gwas_input$output_suffix,
                eths = !!eths, eth_sample_file = paste0("analysis/QC/12_ethnicity_admixture/pca/", !!study_name, "_ETH_samples.txt"),  ## this is not a real file - "ETH" gets replaced by proper "ETH" in `run_gwas`
                eth_low_maf_file = paste0("analysis/QC/14_mafcheck/", !!study_name, "_ETH_low_maf_snps.txt"), ## this is not a real file - "ETH" gets replaced by proper "ETH" in `run_gwas`
                remove_sample_file = file_in(!!paste0("analysis/QC/11_relatedness/", study_name, "_related_ids.txt")),
                output_dir = (!!plink_output_dir), output = !!study_name)},
    dynamic = map(subgroup_gwas_input)
  ),

  # run meta analyses ---------------------------

  linear_meta_out = target({
    linear_out
    meta(output = !!study_name, output_suffix = "", eths = !!eths,
      pheno_list = baseline_gwas_input$output_suffix, threads = 8, output_dir = (!!plink_output_dir))},
    dynamic = map(baseline_gwas_input)),
  interaction_meta_out = target({
    interaction_out
    meta(output = !!study_name, output_suffix = interaction_gwas_input$output_suffix, eth = !!eths,
      pheno_list = !!interaction_outcome, type = "interaction", threads = 8,
      interaction_var = interaction_gwas_input$output_suffix, output_dir = (!!plink_output_dir))},
    dynamic = map(interaction_gwas_input)),
  subgroup_meta_out = target({
    subgroup_out
    meta(output = !!study_name, output_suffix = subgroup_gwas_input$output_suffix, eth = !!eths,
      pheno_list = !!interaction_outcome, type = "subgroup", threads = 8, output_dir = (!!plink_output_dir))},
    dynamic = map(subgroup_gwas_input)),

)

# config <- drake_config(init_analysis)
# vis_drake_graph(config)


### troubleshooting
#
# run_gwas(pfile = ("analysis/QC/15_final_processing/FULL/PSYMETAB_GWAS.FULL"), pheno_file = file_in("data/processed/phenotype_data/GWAS_input/pheno_input.txt"),
#            covar_file = file_in("data/processed/phenotype_data/GWAS_input/covar_input.txt"),
#            threads = 8, pheno_name = baseline_gwas_info$pheno[10], covar_names = baseline_gwas_info$covars[10], eths = eths, output_suffix = baseline_gwas_info$output_suffix[10],
#            eth_sample_file = "analysis/QC/12_ethnicity_admixture/pca/PSYMETAB_GWAS_ETH_samples.txt", ## this is not a real file - "ETH" gets replaced by proper "ETH" in `run_gwas`
#            eth_low_maf_file = paste0("analysis/QC/14_mafcheck/", study_name, "_ETH_low_maf_snps.txt"),
#            remove_sample_file = file_in(paste0("analysis/QC/11_relatedness/", study_name, "_related_ids.txt")),
#            output_dir = ("analysis/GWAS"), output = "PSYMETAB_GWAS")
#
#
# run_gwas(pfile = ("analysis/QC/15_final_processing/FULL/PSYMETAB_GWAS.FULL"), pheno_file = file_in("data/processed/phenotype_data/GWAS_input/pheno_input.txt"),
#           covar_file = file_in("data/processed/phenotype_data/GWAS_input/covar_input.txt"),
#           threads = 8, pheno_name = (subgroup_gwas_info$pheno[2]), covar_names = subgroup_gwas_info$covars[2], subgroup_var = subgroup_gwas_info$subgroup[2], eths = eths,  type = "subgroup", output_suffix = subgroup_gwas_info$output_suffix[2],
#           eth_sample_file = "analysis/QC/12_ethnicity_admixture/pca/PSYMETAB_GWAS_ETH_samples.txt", ## this is not a real file - "ETH" gets replaced by proper "ETH" in `run_gwas`
#           eth_low_maf_file = paste0("analysis/QC/14_mafcheck/", study_name, "_ETH_low_maf_snps.txt"),
#           remove_sample_file = file_in(paste0("analysis/QC/11_relatedness/", study_name, "_related_ids.txt")),
#           output_dir = ("analysis/GWAS"), output = "PSYMETAB_GWAS")
#
# run_gwas2(pfile = ("analysis/QC/15_final_processing/FULL/PSYMETAB_GWAS.FULL"), pheno_file = file_in("data/processed/phenotype_data/GWAS_input/pheno_input.txt"),
#           covar_file = file_in("data/processed/phenotype_data/GWAS_input/covar_input.txt"),
#           threads = 8, pheno_name = (interaction_gwas_info$pheno[4]), covar_names = interaction_gwas_info$covars[4], parameters = interaction_gwas_info$parameters[4], eths = eths,  type = "interaction", output_suffix = interaction_gwas_info$output_suffix[4],
#           eth_sample_file = "analysis/QC/12_ethnicity_admixture/pca/PSYMETAB_GWAS_ETH_samples.txt", ## this is not a real file - "ETH" gets replaced by proper "ETH" in `run_gwas`
#           eth_low_maf_file = paste0("analysis/QC/14_mafcheck/", study_name, "_ETH_low_maf_snps.txt"),
#           remove_sample_file = file_in(paste0("analysis/QC/11_relatedness/", study_name, "_related_ids.txt")),
#           output_dir = ("analysis/GWAS"), output = "PSYMETAB_GWAS")
#
#
# meta(output = study_name, output_suffix = interaction_gwas_info$output_suffix[3], eth = eths,
#   pheno = interaction_gwas_info$pheno[3], type = "interaction", threads = 8, interaction_var = interaction_gwas_info$drug[3])
#
#
# pfile = ("analysis/QC/15_final_processing/FULL/PSYMETAB_GWAS.FULL")
# pheno_file = file_in("data/processed/phenotype_data/GWAS_input/pheno_input.txt")
# covar_file = file_in("data/processed/phenotype_data/GWAS_input/covar_input.txt")
# threads = 8
# pheno_name = subgroup_gwas_info$pheno[2]
# covar_names = subgroup_gwas_info$covars[1]
# eths = eths
# output_suffix = subgroup_gwas_info$output_suffix[2]
# eth_sample_file = "analysis/QC/12_ethnicity_admixture/pca/PSYMETAB_GWAS_ETH_samples.txt" ## this is not a real file - "ETH" gets replaced by proper "ETH" in `run_gwas`
# eth_low_maf_file = "analysis/QC/14_mafcheck/PSYMETAB_GWAS_ETH_low_maf_snps.txt" ## this is not a real file - "ETH" gets replaced by proper "ETH" in `run_gwas`
# remove_sample_file = file_in(paste0("analysis/QC/11_relatedness/", study_name, "_related_ids.txt"))
# output_dir = ("analysis/GWAS")
# output = "PSYMETAB_GWAS"
# type="subgroup"
# subgroup_var = subgroup_gwas_info$subgroup[2]
#

process_init <- drake_plan(

  # define endpoint, covars and outputs ---------------------------

  GWAS_input_process = target(list(full_pheno = read_table2(file_in("data/processed/phenotype_data/GWAS_input/pheno_input.txt")),
                                    full_covar = read_table2(file_in("data/processed/phenotype_data/GWAS_input/covar_input.txt"))),
    hpc = FALSE), # this should be identical to GWAS_input above, but needs to be a different name

  baseline_gwas_process = target(define_baseline_models(GWAS_input_process, !!baseline_vars, !!drug_classes, !!caffeine_vars, !!interaction_outcome),
    hpc = FALSE),
  interaction_gwas_process = target(define_interaction_inputs(!!drug_classes),
    hpc = FALSE),
  subgroup_gwas_process = target(define_subgroup_inputs(!!drug_classes),
    hpc = FALSE),

  baseline_gwas_files = target({
    define_baseline_files(baseline_gwas_process, output = !!study_name, eths = !!eths,
    output_dir = "analysis/GWAS", type = "full")},
    hpc = FALSE),
  interaction_gwas_files = target({
    define_interaction_files(interaction_gwas_process,  output = !!study_name, eths = !!eths,
    output_dir = "analysis/GWAS", type = "interaction", pheno_list = !!interaction_outcome)},
    hpc = FALSE),
  subgroup_gwas_files = target({
    define_subgroup_files(subgroup_gwas_process, output = !!study_name, eths = !!eths,
    output_dir = "analysis/GWAS", type = "subgroup", pheno_list = !!interaction_outcome)},
    hpc = FALSE),

  check_baseline_files = target({
    baseline_gwas_files$log_file},
    dynamic=map(baseline_gwas_files), format = "file"),

  check_interaction_files = target({
    interaction_gwas_files$log_file},
    dynamic=map(interaction_gwas_files), format = "file"),

  check_subgroup_files = target({
    subgroup_gwas_files$log_file},
    dynamic=map(subgroup_gwas_files), format = "file"),

  # write formatted GWAS file ---------------------------

  process_baseline_gwas = target({
    check_baseline_files
    process_gwas(eth = baseline_gwas_files$eth, pheno=baseline_gwas_files$pheno, drug=baseline_gwas_files$drug, file=baseline_gwas_files$plink_file,
      output = !!study_name, output_dir = "analysis/GWAS", type = "full",
      info_file = paste0("analysis/QC/15_final_processing/", !!study_name, ".info"), out_file = baseline_gwas_files$processed_file)
    baseline_gwas_files$processed_file
    },
    dynamic = map(baseline_gwas_files), format = "file"
  ),

  baseline_gwas_figures_input = target({
    process_baseline_gwas
    gwas_figures_input(eth = baseline_gwas_files$eth, pheno=baseline_gwas_files$pheno, drug=baseline_gwas_files$drug, file=baseline_gwas_files$plink_file,
      output = !!study_name, output_dir = "analysis/GWAS", type = "full",
      info_file = paste0("analysis/QC/15_final_processing/", !!study_name, ".info"), out_file = baseline_gwas_files$processed_file)},
    dynamic = map(baseline_gwas_files), hpc = FALSE),

  # create manhattan and qq figures ---------------------------

  baseline_gwas_figures = target({
    process_baseline_gwas
    create_figures(
      baseline_gwas_figures_input$joint_file, baseline_gwas_figures_input$manhattan_file_name,
      baseline_gwas_figures_input$qq_file_name, baseline_gwas_figures_input$title, !!gw_sig)

      c(baseline_gwas_figures_input$manhattan_file_name, baseline_gwas_figures_input$qq_file_name)
    },
    dynamic = map(baseline_gwas_figures_input) ,format = "file"
  ),

  # write formatted GWAS file ---------------------------

  process_interaction_gwas = target({
    check_interaction_files
    process_gwas(eth = interaction_gwas_files$eth, pheno=interaction_gwas_files$pheno, drug=interaction_gwas_files$drug, file=interaction_gwas_files$plink_file,
      output = !!study_name, output_dir = "analysis/GWAS", type = "interaction",
      info_file = paste0("analysis/QC/15_final_processing/", !!study_name, ".info"), out_file = interaction_gwas_files$processed_file)
    interaction_gwas_files$processed_file
    },
    dynamic = map(interaction_gwas_files), format = "file"
  ),

  interaction_gwas_figures_input = target({
    process_interaction_gwas
    gwas_figures_input(eth = interaction_gwas_files$eth, pheno=interaction_gwas_files$pheno, drug=interaction_gwas_files$drug, file=interaction_gwas_files$plink_file,
      output = !!study_name, output_dir = "analysis/GWAS", type = "interaction",
      info_file = paste0("analysis/QC/15_final_processing/", !!study_name, ".info"), out_file = interaction_gwas_files$processed_file, drug_class = interaction_gwas_files$drug_class)},
    dynamic = map(interaction_gwas_files), hpc = FALSE),

  # create manhattan and qq figures ---------------------------

  interaction_gwas_figures = target({
    process_interaction_gwas
    create_figures(interaction_gwas_figures_input$joint_file, interaction_gwas_figures_input$manhattan_file_name,
      interaction_gwas_figures_input$qq_file_name, interaction_gwas_figures_input$title, !!gw_sig)
    c(interaction_gwas_figures_input$manhattan_file_name, interaction_gwas_figures_input$qq_file_name)
    },
    dynamic = map(interaction_gwas_figures_input),
    format = "file"
  ),

  # write formatted GWAS file ---------------------------

  process_subgroup_gwas = target({
    check_subgroup_files
    process_gwas(eth = subgroup_gwas_files$eth, pheno=subgroup_gwas_files$pheno, drug=subgroup_gwas_files$drug, file=subgroup_gwas_files$plink_file,
      output = !!study_name, output_dir = "analysis/GWAS", type = "subgroup",
      info_file = paste0("analysis/QC/15_final_processing/", !!study_name, ".info"), out_file = subgroup_gwas_files$processed_file)
    subgroup_gwas_files$processed_file
    },
    dynamic = map(subgroup_gwas_files),
    format = "file"
  ),

  subgroup_gwas_figures_input = target({
    process_subgroup_gwas
    gwas_figures_input(eth = subgroup_gwas_files$eth, pheno=subgroup_gwas_files$pheno, drug=subgroup_gwas_files$drug, file=subgroup_gwas_files$plink_file,
      output = !!study_name, output_dir = "analysis/GWAS", type = "subgroup",
      info_file = paste0("analysis/QC/15_final_processing/", !!study_name, ".info"), out_file = subgroup_gwas_files$processed_file, drug_class = subgroup_gwas_files$drug_class)},
    dynamic = map(subgroup_gwas_files), hpc = FALSE),

  # create manhattan and qq figures ---------------------------

  subgroup_gwas_figures = target({
    process_subgroup_gwas
    create_figures(subgroup_gwas_figures_input$joint_file, subgroup_gwas_figures_input$manhattan_file_name,
      subgroup_gwas_figures_input$qq_file_name, subgroup_gwas_figures_input$title, !!gw_sig)
    c(subgroup_gwas_figures_input$manhattan_file_name, subgroup_gwas_figures_input$qq_file_name)
    },
    dynamic = map(subgroup_gwas_figures_input),
    format = "file"
  ),

  # extract significant results ---------------------------

  baseline_gwas_sig = target({

    gw_sig_extract(gwas_file = baseline_gwas_files$processed_file, !!gw_sig_nominal, eth = baseline_gwas_files$eth,
      pheno = baseline_gwas_files$pheno, drug_class = NA, drug = baseline_gwas_files$drug)
    },
    dynamic = map(baseline_gwas_files)
  ),

  interaction_gwas_sig = target({

    gw_sig_extract(gwas_file = interaction_gwas_files$processed_file, !!gw_sig_nominal, eth = interaction_gwas_files$eth,
      pheno = interaction_gwas_files$pheno, drug_class = interaction_gwas_files$drug_class, drug = interaction_gwas_files$drug)
    },
    dynamic = map(interaction_gwas_files)
  ),

  subgroup_gwas_sig = target({

    gw_sig_extract(gwas_file = subgroup_gwas_files$processed_file, !!gw_sig_nominal, eth = subgroup_gwas_files$eth,
      pheno = subgroup_gwas_files$pheno, drug_class = subgroup_gwas_files$drug_class, drug = subgroup_gwas_files$drug)
    },
    dynamic = map(subgroup_gwas_files)
  ),

  baseline_gwas_nom_com = bind_rows(baseline_gwas_sig, .id = "pheno_name"),
  interaction_gwas_nom_com = bind_rows(interaction_gwas_sig, .id = "pheno_name"),
  subgroup_gwas_nom_com = bind_rows(subgroup_gwas_sig, .id = "pheno_name"),

  baseline_gwas_sig_com = baseline_gwas_nom_com %>% filter(P < !!gw_sig),
  interaction_gwas_sig_com = interaction_gwas_nom_com %>% filter(P < !!gw_sig),
  subgroup_gwas_sig_com = subgroup_gwas_nom_com %>% filter(P < !!gw_sig),


  baseline_interest = baseline_gwas_sig_com %>% filter(eth=="CEU") %>% filter(!grepl("sensitivity", pheno)),
  subgroup_interest = subgroup_gwas_sig_com %>% filter(drug=="Drug") %>% filter(eth=="CEU") %>% filter(!grepl("sensitivity", drug_class)) %>% filter(grepl("BMI", pheno)),

  write_baseline_interest = write.csv(baseline_interest,  file_out("output/PSYMETAB_GWAS_baseline_CEU_result.csv"),row.names = F),
  write_subgroup_interest = write.csv(subgroup_interest,  file_out("output/PSYMETAB_GWAS_subgroup_CEU_result.csv"),row.names = F),


  baseline_GWAS_count = count_GWAS_n(psam_file = file_in(!!paste0("analysis/QC/15_final_processing/FULL/", study_name, ".FULL.psam")),
    pheno_file = file_in("data/processed/phenotype_data/GWAS_input/baseline_input_resid.txt"), output_suffix = "baseline",
    subgroup = NA, covar_file = NA, covars = NA,
    eths = !!eths, eth_sample_file = paste0("analysis/QC/12_ethnicity_admixture/pca/", !!study_name, "_ETH_samples.txt"),  ## this is not a real file - "ETH" gets replaced by proper "ETH" in `run_gwas`
    related_ids_file = file_in(!!paste0("analysis/QC/11_relatedness/", study_name, "_related_ids.txt"))),


  interaction_GWAS_count = target({
    count_GWAS_n(psam_file = file_in(!!paste0("analysis/QC/15_final_processing/FULL/", study_name, ".FULL.psam")),
      pheno_file = interaction_gwas_process$pheno_file, output_suffix = interaction_gwas_process$output_suffix,
      subgroup = NA, covar_file = interaction_gwas_process$covar_file, covars = interaction_gwas_process$covars,
      eths = !!eths, eth_sample_file = paste0("analysis/QC/12_ethnicity_admixture/pca/", !!study_name, "_ETH_samples.txt"),  ## this is not a real file - "ETH" gets replaced by proper "ETH" in `run_gwas`
      related_ids_file = file_in(!!paste0("analysis/QC/11_relatedness/", study_name, "_related_ids.txt")))
    }, dynamic = map(interaction_gwas_process)
  ),

  subgroup_GWAS_count = target({
    count_GWAS_n(psam_file = file_in(!!paste0("analysis/QC/15_final_processing/FULL/", study_name, ".FULL.psam")),
      pheno_file = subgroup_gwas_process$pheno_file, output_suffix = subgroup_gwas_process$output_suffix,
      subgroup = subgroup_gwas_process$subgroup, covar_file = NA, covars=NA,
      eths = !!eths, eth_sample_file = paste0("analysis/QC/12_ethnicity_admixture/pca/", !!study_name, "_ETH_samples.txt"),  ## this is not a real file - "ETH" gets replaced by proper "ETH" in `run_gwas`
      related_ids_file = file_in(!!paste0("analysis/QC/11_relatedness/", study_name, "_related_ids.txt")))
  }, dynamic = map(subgroup_gwas_process)),



  ukbb_org = ukb_df_mod("ukb21067", path = !!paste0(UKBB_dir, "/org")),
  ukb_key = ukb_df_field("ukb21067", path = !!paste0(UKBB_dir, "/org")),
  sqc = fread(file_in(!!paste0(UKBB_dir,"/geno/ukb_sqc_v2.txt")), header=F, data.table=F),
  fam = fread(file_in(!!paste0(UKBB_dir,"/plink/_001_ukb_cal_chr9_v2.fam")), header=F,data.table=F),
  relatives = read.table(file_in(!!paste0(UKBB_dir,"/geno/",ukbb_relatives_file)), header=T),
  exclusion_file = fread(file_in(!!paste0(UKBB_dir, "/org/", ukbb_exclusion_list)), data.table=F),
  bgen_file = fread(file_in(!!paste0(UKBB_dir, "/", ukbb_sample_file)), header=T,data.table=F),


  #med_codes = read_tsv(file_in(!!medication_codes)),
  qc_data = ukb_gen_sqc_names(sqc),

  sqc_munge = munge_sqc(sqc,fam),
  british_subset = get_british_ids(qc_data, fam),
  ukbb_munge = munge_ukbb(ukbb_org, ukb_key, date_followup, bmi_var, sex_var, age_var),

  related_IDs_remove = ukb_gen_samples_to_remove(relatives, ukb_with_data = as.integer(ukbb_munge$eid)),
  ukbb_filter = filter_ukbb(ukbb_munge, related_IDs_remove, exclusion_file[,1], british_subset),
  ukbb_resid = resid_ukbb(ukbb_filter, ukb_key, sqc_munge, outcome = "bmi_slope", bmi_var, sex_var, age_var),
  ukbb_bgen_order = order_bgen(bgen_file, ukbb_resid, variable = "bmi_slope_res_ivt"),
  ukbb_bgen_out = write.table(ukbb_bgen_order, file_out("data/processed/ukbb_data/BMI_slope"), sep=" ", quote=F, row.names=F,col.names = T),

  chr_num = tibble(chr = 1:22),

  bgenie_out = target({
    launch_bgenie(chr_num$chr, phenofile = file_in("data/processed/ukbb_data/BMI_slope"), threads=1)
    paste0("analysis/GWAS/UKBB/chr", chr_num$chr, ".out.gz")
  }, dynamic = map(chr_num), format = "file"),


  bgenie_unzip = target({
    bgenie_out
    unzip_bgenie(chr_num$chr)
    paste0("analysis/GWAS/UKBB/chr", chr_num$chr, ".out")
    }, dynamic = map(chr_num), format = "file"),

  bgenie_merge = {
    bgenie_unzip
    merge_bgenie_output()
  },

  subgroup_files = list.files(path="analysis/GWAS/subgroup/processed", pattern=".txt", full.names=T),

  BMI_slope_sub_files = subgroup_files[grepl("CEU[.]Drug[.]BMI_slope", subgroup_files) & !grepl("weight", subgroup_files) & !grepl( "sensitivity", subgroup_files)],

  BMI_slope_nominal = target(extract_sig_files(BMI_slope_sub_files, !!gw_sig_nominal),
    dynamic = map(BMI_slope_sub_files)),

  BMI_slope_nominal_AF = add_AF(BMI_slope_nominal, AF_file = file_in("analysis/QC/14_mafcheck/CEU/PSYMETAB_GWAS.CEU.afreq")),

  psy_ukbb_merge = {
    colnames(bgenie_merge) <- paste(colnames(bgenie_merge), "UKBB", sep = "_")
    psy_ukbb_merge <- left_join(BMI_slope_nominal_AF, bgenie_merge, by = c("SNP" = "rsid_UKBB"))
    psy_ukbb_merge
  },

  psy_ukbb_merge_AF = calc_het(psy_ukbb_merge, "SNP", "BETA", "bmi_slope_res_ivt_beta_UKBB", "SE", "bmi_slope_res_ivt_se_UKBB"),

  psy_ukbb_merge_AF_prune = prune_psy_ukbb(psy_ukbb_merge_AF),

  write_ukbb_comparison = write.csv(psy_ukbb_merge_AF_prune,  file_out("output/PSYMETAB_GWAS_UKBB_comparison.csv"),row.names = F)

  ukbb_top_snps_chr = tibble(chr = psy_ukbb_merge_AF$CHR, rsid = psy_ukbb_merge_AF$SNP),

  ukbb_geno = load_geno(bgen_file, ukbb_top_snps_chr),


  # ukbb_comparison = read.csv(file_in("output/PSYMETAB_GWAS_UKBB_comparison.csv")),
  #
  # ukbb_comparison_format = sort_ukbb_comparison(ukbb_comparison, subgroup_GWAS_count, !!interaction_outcome, !!drug_classes),
  # write.csv(ukbb_comparison_format,  file_out("output/PSYMETAB_GWAS_UKBB_comparison2.csv"),row.names = F)

)

#
# ## GARBAGE
# process_gwas(eth = baseline_gwas_files$eth[1], pheno=baseline_gwas_files$pheno[1], drug=baseline_gwas_files$drug[1], file=baseline_gwas_files$file[1],
#   output = !!study_name, output_dir = "analysis/GWAS", type = "full",
#   info_file = paste0("analysis/QC/15_final_processing/", study_name, ".info"), out_file = baseline_gwas_files$write_file[1])
#
# process_gwas(eth = interaction_gwas_files$eth[1], pheno=interaction_gwas_files$pheno[1], drug=interaction_gwas_files$drug[1], file=interaction_gwas_files$file[1],
#   output = !!study_name, output_dir = "analysis/GWAS", type = "interaction",
#   info_file = paste0("analysis/QC/15_final_processing/", study_name, ".info"), out_file = interaction_gwas_files$write_file[1])
#
# process_gwas(eth = subgroup_gwas_files$eth[510], pheno=subgroup_gwas_files$pheno[510], drug=subgroup_gwas_files$drug[510], file=subgroup_gwas_files$file[510],
#   output = study_name, output_dir = "analysis/GWAS", type = "subgroup",
#   info_file = paste0("analysis/QC/15_final_processing/", study_name, ".info"), out_file = subgroup_gwas_files$write_file[510])

# create_figures(
#   baseline_gwas_figures_input$joint_file[1], baseline_gwas_figures_input$manhattan_file_name[1],
#   baseline_gwas_figures_input$qq_file_name[1], baseline_gwas_figures_input$title[1])

######


##### Run GWAS in UKBB --------

ukbb_analysis <- drake_plan(

  ukbb_org = ukb_df_mod("ukb21067", path = !!paste0(UKBB_dir, "/org")),
  ukb_key = ukb_df_field("ukb21067", path = !!paste0(UKBB_dir, "/org")),
  sqc = fread(file_in(!!paste0(UKBB_dir,"/geno/ukb_sqc_v2.txt")), header=F, data.table=F),
  fam = fread(file_in(!!paste0(UKBB_dir,"/plink/_001_ukb_cal_chr9_v2.fam")), header=F,data.table=F),
  relatives = read.table(file_in(!!paste0(UKBB_dir,"/geno/",ukbb_relatives_file)), header=T),
  exclusion_file = fread(file_in(!!paste0(UKBB_dir, "/org/", ukbb_exclusion_list)), data.table=F),
  bgen_file = fread(file_in(!!paste0(UKBB_dir, "/", ukbb_sample_file)), header=T,data.table=F),


  #med_codes = read_tsv(file_in(!!medication_codes)),
  qc_data = ukb_gen_sqc_names(sqc),

  sqc_munge = munge_sqc(sqc,fam),
  british_subset = get_british_ids(qc_data, fam),
  ukbb_munge = munge_ukbb(ukbb_org, ukb_key, date_followup, bmi_var, sex_var, age_var),

  related_IDs_remove = ukb_gen_samples_to_remove(relatives, ukb_with_data = as.integer(ukbb_munge$eid)),
  ukbb_filter = filter_ukbb(ukbb_munge, related_IDs_remove, exclusion_file[,1], british_subset),
  ukbb_resid = resid_ukbb(ukbb_filter, ukb_key, sqc_munge, outcome = "bmi_slope", bmi_var, sex_var, age_var),
  ukbb_bgen_order = order_bgen(bgen_file, ukbb_resid, variable = "bmi_slope_res_ivt"),
  ukbb_bgen_out = write.table(ukbb_bgen_order, file_out("data/processed/ukbb_data/BMI_slope"), sep=" ", quote=F, row.names=F,col.names = T),

  chr_num = tibble(chr = 1:22),

  bgenie_out = target(launch_bgenie(chr_num$chr, phenofile = file_in("data/processed/ukbb_data/BMI_slope"), threads=1),
    dynamic = map(chr_num)),

  bgenie_unzip = target({
    bgenie_out
    unzip_bgenie(chr_num$chr)
    },
    dynamic = map(chr_num)),

  bgenie_merge = {
    bgenie_unzip
    merge_bgenie_output()
  },

  subgroup_files = list.files(path="analysis/GWAS/subgroup/processed", pattern=".txt", full.names=T),

  BMI_slope_sub_files = subgroup_files[grepl("CEU[.]Drug[.]BMI_slope", subgroup_files) & !grepl("weight", subgroup_files) & !grepl( "sensitivity", subgroup_files)],

  BMI_slope_nominal = target(extract_sig_files(BMI_slope_sub_files, !!gw_sig_nominal),
    dynamic = map(BMI_slope_sub_files)),

  BMI_slope_nominal_AF = add_AF(BMI_slope_nominal, AF_file = file_in("analysis/QC/14_mafcheck/CEU/PSYMETAB_GWAS.CEU.afreq")),

  psy_ukbb_merge = {
    colnames(bgenie_merge) <- paste(colnames(bgenie_merge), "UKBB", sep = "_")
    psy_ukbb_merge <- left_join(BMI_slope_nominal_AF, bgenie_merge, by = c("SNP" = "rsid_UKBB"))
    psy_ukbb_merge
  },

  psy_ukbb_merge_AF = calc_het(psy_ukbb_merge, "SNP", "BETA", "bmi_slope_res_ivt_beta_UKBB", "SE", "bmi_slope_res_ivt_se_UKBB"),

  psy_ukbb_merge_AF_prune = prune_psy_ukbb(psy_ukbb_merge_AF),

  write_ukbb_comparison = write.csv(psy_ukbb_merge_AF_prune,  file_out("output/PSYMETAB_GWAS_UKBB_comparison.csv"),row.names = F)



)


prs <- drake_plan(
  prs_info = target(define_prs_inputs(!!consortia_dir, "analysis/PRS"), hpc = FALSE),

  prs_ukbb_info = target(define_ukbb_inputs(!!Neale_summary_dir, !!ukbb_files, "analysis/PRS"), hpc = FALSE),

  prsice_out = target(
    {
      run_prsice(base_file=prs_info$base_file,
        threads=16, memory="100000", out_file=prs_info$out_file,
        bgen_file=paste0("analysis/QC/15_final_processing/FULL/", !!study_name, ".FULL"),
        sample_file=file_in(!!paste0("analysis/QC/15_final_processing/FULL/", study_name, ".FULL_nosex.sample")))
      c(paste0(prs_info$out_file, ".log"), paste0(prs_info$out_file, ".all.score"), paste0(prs_info$out_file, ".prsice"))
    },
    dynamic = map(prs_info),
    format = "file"
  ),

  prs_format = target(
    {
      prsice_out
      format_prs(all_score_file=paste0(prs_info$out_file, ".all.score"), out_file=paste0(prs_info$out_file, ".format"))
      paste0(prs_info$out_file, ".format")
    },
    dynamic = map(prs_info), hpc = FALSE,
    format = "file"
  ),

  prsice_ukbb_out = target(
    {
      run_prsice(base_file=prs_ukbb_info$base_file,
        threads=16, memory="100000", out_file=prs_ukbb_info$out_file,
        bgen_file=paste0("analysis/QC/15_final_processing/FULL/", !!study_name, ".FULL"),
        sample_file=file_in(!!paste0("analysis/QC/15_final_processing/FULL/", study_name, ".FULL_nosex.sample")),
        snp_col="SNP", chr_col="CHR", effect_allele_col="EFFECT_ALLELE",
        other_allele_col="OTHER_ALLELE", beta_or_col="BETA", p_col="PVAL")
      c(paste0(prs_ukbb_info$out_file, ".log"), paste0(prs_ukbb_info$out_file, ".all.score"), paste0(prs_ukbb_info$out_file, ".prsice"))
    },
    dynamic = map(prs_ukbb_info),
    format = "file"
  ),

  prs_ukbb_format = target(
    {
      prsice_out
        format_prs(all_score_file=paste0(prs_ukbb_info$out_file, ".all.score"), out_file=paste0(prs_ukbb_info$out_file, ".format"))
      paste0(prs_ukbb_info$out_file, ".format")
    },
    dynamic = map(prs_ukbb_info), hpc = FALSE,
    format = "file"
  ),

  # so far this portion of the plan only applies to caffeine data

  analyze_ukbb_prs = target(
    {
      analyze_prs(prs_file = prs_ukbb_format, pheno_file = file_in("data/processed/phenotype_data/GWAS_input/pheno_munge.txt"),
        pc_file = file_in(!!paste0("data/processed/phenotype_data/", study_name, "_PCs.txt")),
        linear_pheno_columns = c("logCAF_caffeine", "logPARAX_caffeine", "logTHEOPH_caffeine", "logTHEOBR_caffeine",
          "logCAFPARAX_caffeine", "logCAFPARAXTHEOPH_caffeine"),
          glm_pheno_columns = c("Sleep_disorder_caffeine"), covars = c("sex", "Age_caffeine", "Age_caffeine_sq", "BMI_caffeine", "Hospital_stay_caffeine", "Diag_caffeine"))
    },
    dynamic = map(prs_ukbb_format)
  ),


  ukbb_prs_out = target(
    {
      out_name <- paste0("output/", gsub(".format", "_analysis.csv", basename(prs_ukbb_format)))
      write_ukbb_prs_analysis(analyze_ukbb_prs[[1]], out_name)
    },
    dynamic = map(analyze_ukbb_prs, prs_ukbb_format),
    format = "file"
  )

)

### PRS analysis

prs_analysis <- drake_plan(
  prs_info_analysis = target(define_prs_inputs(!!consortia_dir, "analysis/PRS"), hpc = FALSE),
  prs_ukbb_info_analysis = target(define_ukbb_inputs(!!Neale_summary_dir, !!ukbb_files, "analysis/PRS"), hpc = FALSE),
  prs_files = target(paste0(prs_info_analysis$out_file, ".format"),
    dynamic = map(prs_info_analysis), hpc = FALSE, format = "file"
  ),
  prs_ukbb_files = target(paste0(prs_ukbb_info_analysis$out_file, ".format"),
    dynamic = map(prs_ukbb_info_analysis), hpc = FALSE,
    format = "file"
  ),

  # merge ukbb_files with ukbb_prs_analysis which adds out_file
  # file target for out_file
  analyze_ukbb_prs = target({
    analyze_prs(prs_file = "analysis/PRS/coffee_consumed_Neale_UKBB.format", pheno_file = file_in("data/processed/phenotype_data/GWAS_input/pheno_munge.txt"),
      pc_file = file_in(!!paste0("data/processed/phenotype_data/", study_name, "_PCs.txt")),
      linear_pheno_columns = c("logCAF_caffeine", "logPARAX_caffeine", "logTHEOPH_caffeine", "logTHEOBR_caffeine",
        "logCAFPARAX_caffeine", "logCAFPARAXTHEOPH_caffeine"),
      glm_pheno_columns = c("Sleep_disorder_caffeine"), covars = c("sex", "Age_caffeine", "Age_caffeine_sq", "BMI_caffeine", "Hospital_stay_caffeine", "Diag_caffeine"))
  })#, dynamic = map(prs_ukbb_format))



)
