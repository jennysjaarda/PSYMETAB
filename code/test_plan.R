
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

  pheno_raw = readr::read_delim(file_in(!!pheno_file), col_types = cols(.default = col_character()), delim = ",") %>% type_convert(),
  caffeine_raw = target(read_excel(file_in(!!caffeine_file), sheet=1) %>% type_convert(),
    hpc = FALSE),
  caffeine_man = target(caffeine_raw %>% mutate(age=replace(age, GEN=="WIFRYNSK" & as.Date(Date)=="2009-01-07", 21)),
    hpc = FALSE),
  # found a mistake that the age of WIFRYNSK was wrong at one instance.

  caffeine_munge = munge_caffeine(caffeine_man),
  bgen_sample_file = target({
    readr::read_delim(file_in(!!paste0("analysis/QC/15_final_processing/FULL/", study_name, ".FULL.sample")),
      col_types = cols(.default = col_character()), delim = " ") %>% type_convert()
    }),
  bgen_nosex_out = write.table(bgen_sample_file[,1:3],file_out(!!paste0("analysis/QC/15_final_processing/FULL/", study_name, ".FULL_nosex.sample")),row.names = F, quote = F, col.names = T),
  pc_raw = read_pcs(file_in(!!pc_dir), !!study_name, !!eths) %>% as_tibble(),

  pheno_munge = munge_pheno(pheno_raw, !!baseline_vars, !!leeway_time, caffeine_munge, !!follow_up_limit), # pheno_munge %>% count(GEN) %>% filter(n!=1) ## NO DUPLICATES !
  #pheno_munge = munge_pheno(pheno_raw, baseline_vars, leeway_time, caffeine_munge, follow_up_limit)
  pheno_merge = merge_pheno_caffeine(pheno_raw, caffeine_munge, !!anonymization_error), # in the end, we don't use this dataset, but if you want to merge by caffeine with the appropriate date use this data.

  pheno_baseline = inner_join(pc_raw %>% mutate_at("GPCR", as.character), pheno_munge %>% mutate_at("GEN", as.character), by = c("GPCR" = "GEN")) %>%
    replace_na(list(sex='NONE')),
  pheno_eths_out = write.table(pheno_baseline %>% tidyr::separate(FID, c("COUNT", "GPCR"), "_") %>%
    dplyr::select(COUNT,GPCR,eth ), file_out(!!paste0("data/processed/phenotype_data/", study_name, "_inferred_eths.txt")), row.names = F, quote = F, col.names = T),

  pheno_followup = munge_pheno_follow(pheno_baseline, !!test_drugs), #names(pheno_followup) is the names defined in `test_drugs`: tibble
  GWAS_input = create_GWAS_pheno(pheno_baseline, pheno_followup, !!caffeine_vars),

  ## linear model GWAS:
  #linear_pheno = pheno_baseline %>% dplyr::select(matches('^FID$|^IID$|_start_Drug_1')),
  #linear_covar = pheno_baseline %>% dplyr::select(matches('^FID$|^IID$|^sex$|^Age_Drug_1$|Age_sq_Drug_1$|^PC[0-9]$|^PC[1][0-9]$|^PC20')),
  #pheno_followup_split = split_followup(pheno_followup, !!test_drugs),
  #followup_data_write = if(!is.null(create_pheno_folder)){ write_followup_data(pheno_followup_split)},
  out_pheno_munge =  target({
    write.table(pheno_munge, file_out("data/processed/phenotype_data/GWAS_input/pheno_munge.txt"), row.names = F, quote = F, col.names = T)}),

  out_linear_pheno =  target({
    write.table(GWAS_input$full_pheno, file_out("data/processed/phenotype_data/GWAS_input/pheno_input.txt"), row.names = F, quote = F, col.names = T)}),
  out_linear_covar =  target({
    write.table(GWAS_input$full_covar, file_out("data/processed/phenotype_data/GWAS_input/covar_input.txt"), row.names = F, quote = F, col.names = T)}),

)

# config <- drake_config(analysis_prep)
# vis_drake_graph(config)

init_analysis <- drake_plan(

  GWAS_input_analysis = target(list(full_pheno = read_table2(file_in("data/processed/phenotype_data/GWAS_input/pheno_input.txt")),
                                    full_covar = read_table2(file_in("data/processed/phenotype_data/GWAS_input/covar_input.txt"))),
    hpc = FALSE), # this should be identical to GWAS_input above, but needs to be a different name

  # define endpoint, covars and outputs ---------------------------
  baseline_gwas_info = target(define_baseline_inputs(GWAS_input_analysis, !!baseline_vars, !!drug_classes, !!caffeine_vars),
    hpc = FALSE),
  interaction_gwas_info = target(define_interaction_inputs(GWAS_input_analysis, !!drug_classes),
    hpc = FALSE),
  subgroup_gwas_info = target(define_subgroup_inputs(GWAS_input_analysis, !!drug_classes),
    hpc = FALSE),

  # run initial GWAS ---------------------------
  linear_out = target({
    #loadd(baseline_gwas_info
    file_in(!!paste0("analysis/QC/15_final_processing/FULL/", study_name, ".FULL.log"))
    run_gwas(pfile = paste0("analysis/QC/15_final_processing/FULL/", !!study_name, ".FULL"), pheno_file = file_in("data/processed/phenotype_data/GWAS_input/pheno_input.txt"),
               covar_file = file_in("data/processed/phenotype_data/GWAS_input/covar_input.txt"),
               threads = 8, pheno_name = baseline_gwas_info$pheno, covar_names = baseline_gwas_info$covars, eths = !!eths, output_suffix = baseline_gwas_info$output_suffix,
               eth_sample_file = paste0("analysis/QC/12_ethnicity_admixture/pca/", !!study_name, "_ETH_samples.txt"), ## this is not a real file - "ETH" gets replaced by proper "ETH" in `run_gwas`
               eth_low_maf_file = paste0("analysis/QC/14_mafcheck/", !!study_name, "_ETH_low_maf_snps.txt"), ## this is not a real file - "ETH" gets replaced by proper "ETH" in `run_gwas`
               remove_sample_file = file_in(!!paste0("analysis/QC/11_relatedness/", study_name, "_related_ids.txt")),
               output_dir = ("analysis/GWAS"), output = !!study_name)},
    dynamic = map(baseline_gwas_info)
  ),

  interaction_out = target({
    #loadd(interaction_gwas_info)
    file_in(!!paste0("analysis/QC/15_final_processing/FULL/", study_name, ".FULL.log"))
    run_gwas(pfile = paste0("analysis/QC/15_final_processing/FULL/", !!study_name, ".FULL"), pheno_file = file_in("data/processed/phenotype_data/GWAS_input/pheno_input.txt"),
                covar_file = file_in("data/processed/phenotype_data/GWAS_input/covar_input.txt"),type = "interaction",
                threads = 8, pheno_name = interaction_gwas_info$pheno, covar_names = interaction_gwas_info$covars, parameters = interaction_gwas_info$parameters, output_suffix = interaction_gwas_info$output_suffix,
                eths = !!eths, eth_sample_file = paste0("analysis/QC/12_ethnicity_admixture/pca/", !!study_name, "_ETH_samples.txt"),  ## this is not a real file - "ETH" gets replaced by proper "ETH" in `run_gwas`
                eth_low_maf_file = paste0("analysis/QC/14_mafcheck/", !!study_name, "_ETH_low_maf_snps.txt"), ## this is not a real file - "ETH" gets replaced by proper "ETH" in `run_gwas`
                remove_sample_file = file_in(!!paste0("analysis/QC/11_relatedness/", study_name, "_related_ids.txt")),
                output_dir = ("analysis/GWAS"), output = !!study_name)},
    dynamic = map(interaction_gwas_info)
  ),

  subgroup_out = target({
    #loadd(subgroup_gwas_info)
    file_in(!!paste0("analysis/QC/15_final_processing/FULL/", study_name, ".FULL.log"))
    run_gwas(pfile = paste0("analysis/QC/15_final_processing/FULL/", !!study_name, ".FULL"), pheno_file = file_in("data/processed/phenotype_data/GWAS_input/pheno_input.txt"),
                covar_file = file_in("data/processed/phenotype_data/GWAS_input/covar_input.txt"),
                threads = 8, pheno_name = subgroup_gwas_info$pheno, covar_names = subgroup_gwas_info$covars, subgroup_var = subgroup_gwas_info$subgroup, type = "subgroup", output_suffix = subgroup_gwas_info$output_suffix,
                eths = !!eths, eth_sample_file = paste0("analysis/QC/12_ethnicity_admixture/pca/", !!study_name, "_ETH_samples.txt"),  ## this is not a real file - "ETH" gets replaced by proper "ETH" in `run_gwas`
                eth_low_maf_file = paste0("analysis/QC/14_mafcheck/", !!study_name, "_ETH_low_maf_snps.txt"), ## this is not a real file - "ETH" gets replaced by proper "ETH" in `run_gwas`
                remove_sample_file = file_in(!!paste0("analysis/QC/11_relatedness/", study_name, "_related_ids.txt")),
                output_dir = ("analysis/GWAS"), output = !!study_name)},

    dynamic = map(subgroup_gwas_info)
  ),
  linear_meta_out = target({
    linear_out
    meta(output = !!study_name, output_suffix = baseline_gwas_info$output_suffix, eths = !!eths,
      pheno = baseline_gwas_info$pheno, threads = 8)},
    dynamic = map(baseline_gwas_info)),
  interaction_meta_out = target({
    interaction_out
    meta(output = !!study_name, output_suffix = interaction_gwas_info$output_suffix, eth = !!eths,
      pheno = interaction_gwas_info$pheno, type = "interaction", threads = 8)},
    dynamic = map(interaction_gwas_info)),
  subgroup_meta_out = target({
    subgroup_out
    meta(output = !!study_name, output_suffix = subgroup_gwas_info$output_suffix, eth = !!eths,
      pheno = subgroup_gwas_info$pheno, type = "subgroup", threads = 8)},
    dynamic = map(subgroup_gwas_info)),

)



# config <- drake_config(init_analysis)
# vis_drake_graph(config)



 ### GARBAGE
#
# run_gwas(pfile = ("analysis/QC/15_final_processing/FULL/PSYMETAB_GWAS.FULL"), pheno_file = file_in("data/processed/phenotype_data/GWAS_input/pheno_input.txt"),
#            covar_file = file_in("data/processed/phenotype_data/GWAS_input/covar_input.txt"),
#            threads = 8, pheno_name = baseline_gwas_info$pheno[1], covar_names = baseline_gwas_info$covars[1], eths = eths, output_suffix = baseline_gwas_info$output_suffix[1],
#            eth_sample_file = "analysis/QC/12_ethnicity_admixture/pca/PSYMETAB_GWAS_ETH_samples.txt", ## this is not a real file - "ETH" gets replaced by proper "ETH" in `run_gwas`
#            eth_low_maf_file = "analysis/QC/14_mafcheck/pca/PSYMETAB_GWAS_ETH_low_maf_snps.txt", ## this is not a real file - "ETH" gets replaced by proper "ETH" in `run_gwas`
#            remove_sample_file = file_in(paste0("analysis/QC/11_relatedness/", study_name, "_related_ids.txt")),
#            output_dir = ("analysis/GWAS"), output = "PSYMETAB_GWAS")
#
# pfile = ("analysis/QC/15_final_processing/FULL/PSYMETAB_GWAS.FULL")
# pheno_file = file_in("data/processed/phenotype_data/GWAS_input/pheno_input.txt")
# covar_file = file_in("data/processed/phenotype_data/GWAS_input/covar_input.txt")
# threads = 8
# pheno_name = baseline_gwas_info$pheno[1]
# covar_names = baseline_gwas_info$covars[1]
# eths = eths
# output_suffix = baseline_gwas_info$output_suffix[1]
# eth_sample_file = "analysis/QC/12_ethnicity_admixture/pca/PSYMETAB_GWAS_ETH_samples.txt" ## this is not a real file - "ETH" gets replaced by proper "ETH" in `run_gwas`
# eth_low_maf_file = "analysis/QC/14_mafcheck/PSYMETAB_GWAS_ETH_low_maf_snps.txt" ## this is not a real file - "ETH" gets replaced by proper "ETH" in `run_gwas`
# remove_sample_file = file_in(paste0("analysis/QC/11_relatedness/", study_name, "_related_ids.txt"))
# output_dir = ("analysis/GWAS")
# output = "PSYMETAB_GWAS"


process_init <- drake_plan(
  # define endpoint, covars and outputs ---------------------------

  GWAS_input_process = target(list(full_pheno = read_table2(file_in("data/processed/phenotype_data/GWAS_input/pheno_input.txt")),
                                    full_covar = read_table2(file_in("data/processed/phenotype_data/GWAS_input/covar_input.txt"))),
    hpc = FALSE), # this should be identical to GWAS_input above, but needs to be a different name

  # define endpoint, covars and outputs ---------------------------
  baseline_gwas_process = target(define_baseline_inputs(GWAS_input_process, !!baseline_vars, !!drug_classes, !!caffeine_vars),
    hpc = FALSE),
  interaction_gwas_process = target(define_interaction_inputs(GWAS_input_process, !!drug_classes),
    hpc = FALSE),
  subgroup_gwas_process = target(define_subgroup_inputs(GWAS_input_process, !!drug_classes),
    hpc = FALSE),

  baseline_gwas_files = target({
    define_baseline_files(baseline_gwas_process, output = !!study_name, eths,
    output_dir = "analysis/GWAS", type = "full")},
    hpc = FALSE),
  interaction_gwas_files = target({
    define_interaction_files(interaction_gwas_process,  output = !!study_name, eths,
    output_dir = "analysis/GWAS", type = "interaction")},
    hpc = FALSE),
  subgroup_gwas_files = target({
    define_subgroup_files(subgroup_gwas_process, output = !!study_name, eths,
    output_dir = "analysis/GWAS", type = "subgroup")},
    hpc = FALSE),

  check_baseline_files = target({
    baseline_gwas_files$file},
    dynamic=map(baseline_gwas_files), format = "file"),

  check_interaction_files = target({
    interaction_gwas_files$file},
    dynamic=map(interaction_gwas_files), format = "file"),

  check_subgroup_files = target({
    subgroup_gwas_files$file},
    dynamic=map(subgroup_gwas_files), format = "file"),


  process_baseline_gwas = target({ # write formatted GWAS file
    check_baseline_files
    process_gwas(eth = baseline_gwas_files$eth, pheno=baseline_gwas_files$pheno, drug=baseline_gwas_files$drug, file=baseline_gwas_files$file,
    output = !!study_name, output_dir = "analysis/GWAS", type = "full",
    info_file = paste0("analysis/QC/15_final_processing/", !!study_name, ".info"), out_file = baseline_gwas_files$write_file)
    baseline_gwas_files$write_file
    },
    dynamic = map(baseline_gwas_files), format = "file"
  ),
  baseline_gwas_figures_input = target({
    process_baseline_gwas
    gwas_figures_input(eth = baseline_gwas_files$eth, pheno=baseline_gwas_files$pheno, drug=baseline_gwas_files$drug, file=baseline_gwas_files$file,
      output = !!study_name, output_dir = "analysis/GWAS", type = "full",
      info_file = paste0("analysis/QC/15_final_processing/", !!study_name, ".info"), out_file = baseline_gwas_files$write_file)},
    dynamic = map(baseline_gwas_files), hpc = FALSE),

  baseline_gwas_figures = target({
    process_baseline_gwas
    create_figures(
      baseline_gwas_figures_input$joint_file, baseline_gwas_figures_input$manhattan_file_name,
      baseline_gwas_figures_input$qq_file_name, baseline_gwas_figures_input$title)

    #eth = baseline_gwas_files$eth, pheno=baseline_gwas_files$pheno, drug=baseline_gwas_files$drug, file=baseline_gwas_files$file,
    #output = !!study_name, output_dir = "analysis/GWAS", type = "full",
    #info_file = "analysis/QC/15_final_processing/PSYMETAB_GWAS.info", out_file = baseline_gwas_files$write_file)
    c(baseline_gwas_figures_input$manhattan_file_name, baseline_gwas_figures_input$qq_file_name)
    },
    #dynamic = map(baseline_gwas_figures_input)
    dynamic = map(baseline_gwas_figures_input) ,format = "file"
  ),


  process_interaction_gwas = target({
    check_interaction_files
    process_gwas(eth = interaction_gwas_files$eth, pheno=interaction_gwas_files$pheno, drug=interaction_gwas_files$drug, file=interaction_gwas_files$file,
    output = !!study_name, output_dir = "analysis/GWAS", type = "interaction",
    info_file = paste0("analysis/QC/15_final_processing/", !!study_name, ".info"), out_file = interaction_gwas_files$write_file)
    interaction_gwas_files$write_file
    },
    dynamic = map(interaction_gwas_files), format = "file"
  ),
  interaction_gwas_figures_input = target({
    process_interaction_gwas
    gwas_figures_input(eth = interaction_gwas_files$eth, pheno=interaction_gwas_files$pheno, drug=interaction_gwas_files$drug, file=interaction_gwas_files$file,
      output = !!study_name, output_dir = "analysis/GWAS", type = "interaction",
      info_file = paste0("analysis/QC/15_final_processing/", !!study_name, ".info"), out_file = interaction_gwas_files$write_file)},
    dynamic = map(interaction_gwas_files), hpc = FALSE),

  interaction_gwas_figures = target({
    process_interaction_gwas
    create_figures(interaction_gwas_figures_input$joint_file, interaction_gwas_figures_input$manhattan_file_name,
    interaction_gwas_figures_input$qq_file_name, interaction_gwas_figures_input$title)


      #process_interaction_gwas$joint, process_interaction_gwas$manhattan_file_name,
      #process_interaction_gwas$qq_file_name, process_interaction_gwas$title)
    c(interaction_gwas_figures_input$manhattan_file_name, interaction_gwas_figures_input$qq_file_name)
    },
    dynamic = map(interaction_gwas_figures_input),
    format = "file"
  ),

  process_subgroup_gwas = target({
    check_subgroup_files
    process_gwas(eth = subgroup_gwas_files$eth, pheno=subgroup_gwas_files$pheno, drug=subgroup_gwas_files$drug, file=subgroup_gwas_files$file,
    output = !!study_name, output_dir = "analysis/GWAS", type = "subgroup",
    info_file = paste0("analysis/QC/15_final_processing/", !!study_name, ".info"), out_file = subgroup_gwas_files$write_file)
    subgroup_gwas_files$write_file
    },
    dynamic = map(subgroup_gwas_files), format = "file"
  ),
  subgroup_gwas_figures_input = target({
    process_subgroup_gwas
    gwas_figures_input(eth = subgroup_gwas_files$eth, pheno=subgroup_gwas_files$pheno, drug=subgroup_gwas_files$drug, file=subgroup_gwas_files$file,
      output = !!study_name, output_dir = "analysis/GWAS", type = "subgroup",
      info_file = paste0("analysis/QC/15_final_processing/", !!study_name, ".info"), out_file = subgroup_gwas_files$write_file)},
    dynamic = map(subgroup_gwas_files), hpc = FALSE),
  subgroup_gwas_figures = target({
    process_subgroup_gwas
    create_figures(subgroup_gwas_figures_input$joint_file, subgroup_gwas_figures_input$manhattan_file_name,
    subgroup_gwas_figures_input$qq_file_name, subgroup_gwas_figures_input$title)
      #process_subgroup_gwas$joint, process_subgroup_gwas$manhattan_file_name,
      #process_subgroup_gwas$qq_file_name, process_subgroup_gwas$title)
    c(subgroup_gwas_figures_input$manhattan_file_name, subgroup_gwas_figures_input$qq_file_name)
    },
    dynamic = map(subgroup_gwas_figures_input),
    format = "file"
  ),

  calc_gw_sig = target({
    process_subgroup_gwas
    gw_sig(interaction_gwas_figures_input$joint_file)

    gwas <- fread(interaction_gwas_figures_input[[17]])
    sig <- gwas %>% filter(P < 5e-08)

  })


  ## create meta directories
  ## read from sort_gwas.r

  # gwas_figures = target(
  #   process_meta(outcome_variable,interaction_variable,model="subgroup", output_dir = ("analysis/GWAS")),
  #   dynamic = map(subgroup_gwas_info)
  # )

)



#
# ## GARBAGE
# process_gwas(eth = baseline_gwas_files$eth[19], pheno=baseline_gwas_files$pheno[19], drug=baseline_gwas_files$drug[19], file=baseline_gwas_files$file[19],
# output = "PSYMETAB_GWAS", output_dir = "analysis/GWAS", type = "full",
# info_file = "analysis/QC/15_final_processing/PSYMETAB_GWAS.info", out_file = baseline_gwas_files$write_file[19])
#
# create_figures(
#   baseline_gwas_figures_input$joint_file[1], baseline_gwas_figures_input$manhattan_file_name[1],
#   baseline_gwas_figures_input$qq_file_name[1], baseline_gwas_figures_input$title[1])

######

prs <- drake_plan(
  prs_info = target(define_prs_inputs(!!consortia_dir, "analysis/PRS"), hpc = FALSE),
  prs_ukbb_info = target(define_ukbb_inputs(!!Neale_summary_dir, !!ukbb_files, "analysis/PRS"), hpc = FALSE),
  prsice_out = target({
    run_prsice(base_file=prs_info$base_file,
      threads=16, memory="100000", out_file=prs_info$out_file,
      bgen_file=paste0("analysis/QC/15_final_processing/FULL/", !!study_name, ".FULL"),
      sample_file=file_in(!!paste0("analysis/QC/15_final_processing/FULL/", !!study_name, ".FULL_nosex.sample")))
    c(paste0(prs_info$out_file, ".log"), paste0(prs_info$out_file, ".all.score"), paste0(prs_info$out_file, ".prsice"))
  },
    dynamic = map(prs_info),
    format = "file"
  ),
  prs_format = target({
    prsice_out
    format_prs(all_score_file=paste0(prs_info$out_file, ".all.score"), out_file=paste0(prs_info$out_file, ".format"))
    paste0(prs_info$out_file, ".format")
  },
  dynamic = map(prs_info), hpc = FALSE,
  format = "file"
  ),

  prsice_ukbb_out = target({
    run_prsice(base_file=prs_ukbb_info$base_file,
      threads=16, memory="100000", out_file=prs_ukbb_info$out_file,
      bgen_file=paste0("analysis/QC/15_final_processing/FULL/", !!study_name, ".FULL"),
      sample_file=file_in(!!paste0("analysis/QC/15_final_processing/FULL/", !!study_name, ".FULL_nosex.sample")),
      snp_col="SNP", chr_col="CHR", effect_allele_col="EFFECT_ALLELE",
      other_allele_col="OTHER_ALLELE", beta_or_col="BETA", p_col="PVAL")
    c(paste0(prs_ukbb_info$out_file, ".log"), paste0(prs_ukbb_info$out_file, ".all.score"), paste0(prs_ukbb_info$out_file, ".prsice"))
  },
    dynamic = map(prs_ukbb_info),
    format = "file"
  ),
  prs_ukbb_format = target({
    prsice_out
    format_prs(all_score_file=paste0(prs_ukbb_info$out_file, ".all.score"), out_file=paste0(prs_ukbb_info$out_file, ".format"))
    paste0(prs_ukbb_info$out_file, ".format")
  },
  dynamic = map(prs_ukbb_info), hpc = FALSE,
  format = "file"
),

  analyze_ukbb_prs = target({
    analyze_prs(prs_file = prs_ukbb_format, pheno_file = file_in("data/processed/phenotype_data/GWAS_input/pheno_munge.txt"),
      pc_eth_data = pc_raw,
      linear_pheno_columns = c("logCAF", "logPARAX", "logTHEOPH", "logTHEOBR", "logCAFPARAX", "logCAFPARAXTHEOPH"),
      glm_pheno_columns = c("Sleep_disorder"), covars = c("Age_caffeine", "Age_caffeine_sq", "sex"))
  }, dynamic = map(prs_ukbb_format))


)

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
    analys_prs(prs_file = prs_ukbb_format, pheno_file = file_in("data/processed/phenotype_data/GWAS_input/pheno_munge.txt"),
      pc_raw,
      linear_pheno_columns = c("logCAF", "logPARAX", "logTHEOPH", "logTHEOBR", "logCAFPARAX", "logCAFPARAXTHEOPH"),
      glm_pheno_columns = c("Sleep_disorder"), covars = c("Age_caffeine", "Age_caffeine_sq"))
  }, dynamic = map(prs_ukbb_format))


)

#prs_analysis <- bind_plans(analysis_prep, prs)


### A MESS
#
# pull(init_analysis %>% dplyr::select(command))
#
#
#
# post_impute <- bind_plans(pre_impute_qc, init_analysis, process_init)
#
# make(post_impute,parallelism = "future",jobs = 2, console_log_file = "post_impute_qc.out", resources = list(partition = "sgg"))
# make(init_analysis,parallelism = "clustermq",jobs = 4, template = list(cpus = 16, partition = "cluster"))
#
#
#
# ### some useful code
#
# c(init_analysis %>% dplyr::select(target))
# outdated(drake_config(init_analysis))
# make(init_analysis,parallelism = "clustermq",jobs = 4, console_log_file = "init_analysis.out", template = list(cpus = 16, partition = "cluster"))
#
#
#
#
# #### run gwas and GRS computation
# ### to be grabble from `GWAS.sh` and `PRSice.sh`
# )
#
#
#
#
#
#
#
#
#
# ####
# visualize <- drake_plan(
#   dim(dups) #30 2
#   table(sex_info3$Sexe, exclude = NULL)
#   #    F    M
#   # 1298 1469
#   table(eth_info3$Ethnie, exclude = NULL)
#   # africain  africain + caucasien       amerique du sud
#   #      101                     1                     1
#   # Antilles                 arabe     arabe + caucasien
#   #        1                    60                     3
#   # asiatique asiatique + caucasien                 autre
#   #       25                     1                   178
#   # caucasien               inconnu                  <NA>
#   #     1667                   458                   271
#
#   ## double check no duplicates
#   sex_info4 <- unique(sex_info)
#   dim(sex_info)==dim(sex_info4)
#   sex_info4[which(duplicated(sex_info4[,2])),]
#   duplicate_IDs <- sex_info4[which(duplicated(sex_info4[,2])),2]
#   sex_info4[which(sex_info4[,2]==duplicate_IDs),]
#   #TRUE
#
#   eth_info4 <- unique(eth_info)
#   dim(eth_info3)==dim(eth_info4)
#   #TRUE
#   eth_info4[which(duplicated(eth_info4[,2])),]
#   #NONE
#
#     ## anyone missign?
#     no_sex_eth <- fam[which(!fam[,2] %in% sex_info4[,2]), c(1,2)]
#     no_sex_eth
#     # empty
#
#     ### dimensions of PC data
#     print(eth)
#     print(dim(PC_temp))
#
#     drug_list <- unique(unlist(full_pheno %>% dplyr::select(starts_with("AP1_Drug_"))))
#     ### list of drugs included in phenofile
#
#
#
#   create_report(pheno_baseline, config = configure_report(add_plot_prcomp = FALSE))
#
#
# )
