

# Workflow Plans
#

qc_prep <- drake_plan (
  # Load and clean data to prepare for QC ----------------------
  qc_pheno_raw = read_excel(file_in(!!qc_pheno_file), sheet = 1),
  id_code = read.csv(file_in("data/raw/ID_key.csv"), header = T),
  rs_conversion = fread(file_in(!!rsconv_raw_file), data.table = F),
  fam_raw = read.table(file_in(!!paste0(
    plink_bed_out, ".fam"
  ))),
  # read the .fam file
  ped_folder = file_in(!!dirname(plink_ped_raw)),
  #bed_folder = file_out(!!dirname(plink_bed_out)),
  #create_dir = dir.create(bed_folder, showWarnings = F),
  create_bed_out = {
    in_file <- paste0(ped_folder, "/", basename(!!plink_ped_raw))
    outfile <-
      paste0(!!dirname(plink_bed_out), "/", basename(!!plink_bed_out))
    processx::run(
      command = "plink",
      c("--file", in_file, "--make-bed", "--out", outfile),
      error_on_status = F
    )
    file_out(!!paste0(plink_bed_out, ".fam"))
    file_out(!!dirname(plink_bed_out))
  },

  fam_munge = munge_fam(fam_raw, id_code),
  qc_pheno_munge = munge_qc_pheno(qc_pheno_raw, fam_munge),
  dups = find_dups(qc_pheno_munge, fam_munge),
  sex_info = qc_pheno_munge %>% dplyr::select(c(FID, IID, Sexe)),
  eth_info = qc_pheno_munge %>% dplyr::select(c(FID, IID, Ethnie)),
  sex_format = format_sex_file(sex_info),
  eth_format = format_eth_file(eth_info),
  rs_conversion_munge = rs_conversion %>% mutate(RsID = case_when(RsID ==
                                                                    "." ~ Name, TRUE ~ RsID)),
  out_sex = write.table(
    sex_format,
    file_out("data/processed/phenotype_data/PSYMETAB_GWAS_sex.txt"),
    row.names = F,
    quote = F,
    col.names = F
  ),
  out_eth = write.table(
    eth_format,
    file_out("data/processed/phenotype_data/PSYMETAB_GWAS_eth.txt"),
    row.names = F,
    quote = F,
    col.names = F
  ),
  out_dups = write.table(
    dups$dups,
    file_out("data/processed/phenotype_data/PSYMETAB_GWAS_dupIDs.txt"),
    row.names = F,
    quote = F,
    col.names = F
  ),
  out_dups_set = write.table(
    dups$dups_set,
    file_out(
      "data/processed/phenotype_data/PSYMETAB_GWAS_dupIDs_set.txt"
    ),
    row.names = F,
    quote = F,
    col.names = F
  ),
  out_rs_conversion = {
    write.table(
      rs_conversion_munge,
      file_out("data/processed/reference_files/rsid_conversion.txt"),
      row.names = F,
      quote = F,
      col.names = F
    )
  },

)


# vis_drake_graph(qc_prep)

psy_miscellaneous <- drake_plan(
  custom_snp_list = read_excel(file_in("data/raw/custom_SNP_list.xlsx"), sheet = 1),
  HRC_list = fread(
    "/data/sgg2/jenny/data/HRC.r1-1.GRCh37.wgs.mac5.sites.tab",
    data.table = F
  ),
  ## check which SNPs in the custom SNP list (defined by Aurelie and others in the group at Cery) would already be imputed using the HRC panel.
  ## if SNPs are already in the HRC list above, then there is no need to genotype them directly.

  custom_snp_list_non_hrc = custom_snp_list %>% filter(!SNP %in%  HRC_list$ID),
  out_non_hrc_snps = write.table(
    custom_snp_list_non_hrc,
    file_out("data/processed/non_hrc_custom_snps.txt"),
    row.names = F,
    quote = F,
    col.names = T
  ),
)

# vis_drake_graph(psy_miscellaneous)


# run pre-imupation quality control ----------------------



pre_impute_qc <- drake_plan(run_pre_imputation = target({
  file_in("data/processed/phenotype_data/PSYMETAB_GWAS_sex.txt")
  file_in("data/processed/phenotype_data/PSYMETAB_GWAS_eth.txt")

  file_in("data/processed/phenotype_data/PSYMETAB_GWAS_dupIDs.txt")
  file_in("data/processed/phenotype_data/PSYMETAB_GWAS_dupIDs_set.txt")

  file_in("data/processed/reference_files/rsid_conversion.txt")
  processx::run(command = "sh",
                c(pre_imputation_script),
                error_on_status = F)
  file_out(
    "analysis/QC/00_preprocessing",
    "analysis/QC/01_strandflip",
    "analysis/QC/02_maf_zero",
    "analysis/QC/03_missingness",
    "analysis/QC/04_sexcheck"
  )
}), )

# vis_drake_graph(pre_impute_qc)

################################################################
## here download output and upload to Michigan server
################################################################

# download imputation ------------
download_impute <- drake_plan(download_imputation = target({
  processx::run(command = "sh",
                c(file_in(!!download_imputation_script)),
                error_on_status = F)
}),)

# vis_drake_graph(download_impute)


post_impute <- drake_plan(
  # run post-imputation quality control and processing ------------
  run_check_imputation = target({
    # file_in("analysis/QC/06_imputation_get")
    # this file will get zipped later so we cannot track it
    processx::run(command = "sh",
                  c(file_in(!!check_imputation_script)),
                  error_on_status = F)
    file_out("analysis/QC/07_imputation_check")
  }),
  run_post_imputation = target({
    file_in(
      "code/qc/ethnicity_check.R",
      "code/qc/relatedness_filter.R",
      "code/qc/maf_check.R",
      "code/qc/update_pvar.R"
    )
    processx::run(command = "sh",
                  c(file_in(!!post_imputation_script)),
                  error_on_status = F)
    # file_out("analysis/QC/08_plink_convert"), "analysis/QC/09_extract_typed", "analysis/QC/10_merge_imputed", "analysis/QC/11_relatedness",
    #         "analysis/QC/12_ethnicity_admixture", "analysis/QC/13_hwecheck", "analysis/QC/14_mafcheck")
    # these files are too big to track
  }),
  run_final_processing = target({
    run_post_imputation
    file_in(!!(
      paste0(
        "analysis/QC/14_mafcheck/",
        study_name,
        ".mafcheck.step14.log"
      )
    ))
    # file_in("analysis/QC/11_relatedness", "analysis/QC/12_ethnicity_admixture", "analysis/QC/14_mafcheck")
    processx::run(command = "sh",
                  c(file_in(!!final_processing_script)),
                  error_on_status = F)
    # file_out("analysis/QC/15_final_processing")
    # these files are too big to track
  }),
  cp_qc_report = target({
    file_in("analysis/QC/06_imputation_get/qcreport.html")
    processx::run(
      command = "cp",
      c(
        "analysis/QC/06_imputation_get/qcreport.html",
        "docs/generated_reports/"
      ),
      error_on_status = F
    )
    file_out("docs/generated_reports/qcreport.html")
  }, hpc = FALSE),
  cp_qc_check = target({
    file_in("analysis/QC/07_imputation_check")
    processx::run(
      "/bin/sh",
      c(
        "-c",
        "cp analysis/QC/07_imputation_check/summaryOutput/*html docs/generated_reports/"
      ),
      error_on_status = F
    )
    file_out(
      "docs/generated_reports/07_imputation_check.html",
      "docs/generated_reports/CrossCohortReport.html"
    )
  }, hpc = FALSE),
  snp_weights = target({
    run_post_imputation
    read.table((file_in(
      !!paste0(
        "analysis/QC/12_ethnicity_admixture/snpweights/",
        study_name,
        ".NA.predpc"
      )
    )))
  }, hpc = FALSE),
  related_inds = target({
    run_post_imputation
    read.table((file_in(
      !!paste0(
        "analysis/QC/11_relatedness/",
        study_name,
        "_related_ids.txt"
      )
    )))
  }, hpc = FALSE),
  snp_weights_munge = target(munge_snp_weights(snp_weights, related_inds), hpc = FALSE),
  snp_weights_out = target(write.table(snp_weights_munge, file_out(
    !!paste0("data/processed/extractions/", study_name, ".snpweights")
  )), hpc = FALSE),

  ## an info file was also created in the "final_processing" script but it iddin't extract all relevant fields.
  # info_file = target({
  #   run_final_processing
  #   system(!!paste0("bcftools query -f '%CHROM %POS %ID %REF %ALT %INFO/R2{0} %INFO/ER2 %INFO/IMPUTED %INFO/TYPED %INFO/TYPED_ONLY\n' analysis/QC/15_final_processing/FULL/", study_name, ".FULL.vcf > analysis/QC/", study_name, ".impute_info"))
  # }
  # ),

)

# vis_drake_graph(post_impute)


qc_process <- drake_plan(
  imputed_var = count_imputed_var(),

  missingness_variants_removed = {
    countLines(paste0(
      "analysis/QC/02_maf_zero/",
      !!study_name,
      ".maf_zero.step2.bim"
    ))[1] -
      countLines(
        paste0(
          "analysis/QC/03_missingness/",
          !!study_name,
          ".missingness_geno_10.bim"
        )
      )[1] +
      countLines(
        paste0(
          "analysis/QC/03_missingness/",
          !!study_name,
          ".missingness_mind_10.bim"
        )
      )[1] -
      countLines(
        paste0(
          "analysis/QC/03_missingness/",
          !!study_name,
          ".missingness_geno_05.bim"
        )
      )[1] +
      countLines(
        paste0(
          "analysis/QC/03_missingness/",
          !!study_name,
          ".missingness_mind_05.bim"
        )
      )[1] -
      countLines(
        paste0(
          "analysis/QC/03_missingness/",
          !!study_name,
          ".missingness_geno_01.bim"
        )
      )[1]
  },

  missingness_indiv_removed = {
    countLines(paste0(
      "analysis/QC/02_maf_zero/",
      !!study_name,
      ".maf_zero.step2.fam"
    ))[1] -
      countLines(
        paste0(
          "analysis/QC/03_missingness/",
          !!study_name,
          ".missingness_mind_10.fam"
        )
      )[1] +
      countLines(
        paste0(
          "analysis/QC/03_missingness/",
          !!study_name,
          ".missingness_geno_05.fam"
        )
      )[1] -
      countLines(
        paste0(
          "analysis/QC/03_missingness/",
          !!study_name,
          ".missingness_mind_05.fam"
        )
      )[1] +
      countLines(
        paste0(
          "analysis/QC/03_missingness/",
          !!study_name,
          ".missingness_geno_01.fam"
        )
      )[1] -
      countLines(paste0(
        "analysis/QC/03_missingness/",
        !!study_name,
        ".missingness.step3.fam"
      ))[1]
  },

  initial_variants = countLines(
    paste0(
      "analysis/QC/02_maf_zero/",
      !!study_name,
      ".maf_zero.step2.bim"
    )
  )[1],
  initial_ind = countLines(
    paste0(
      "analysis/QC/02_maf_zero/",
      !!study_name,
      ".maf_zero.step2.fam"
    )
  )[1],

  preprocess_var_init = countLines(paste0(!!plink_bed_out, ".bim")),
  preprocess_ind_init = countLines(paste0(!!plink_bed_out, ".bim")),

  preprocess_non_x = countLines(
    paste0("analysis/QC/00_preprocessing/non_autosome_or_x_variants")
  ),
  preprocess_non_x_remain = countLines(paste0(
    "analysis/QC/00_preprocessing/temp1.bim"
  )),

  preprocess_var_dups = countLines(paste0(
    "analysis/QC/00_preprocessing/duplicate_variants"
  )),
  preprocess_var_dups_remain = countLines(paste0(
    "analysis/QC/00_preprocessing/temp2.bim"
  )),


)

analysis_prep <- drake_plan(
  # prepare phenotype files for analysis in GWAS/GRS etc. ------------

  pheno_raw = readr::read_delim(
    file_in(!!pheno_file),
    na = c("#pdt", "#ccannul", "#pam", "#di", "#pnc"),
    col_types = cols(.default = col_character()),
    delim = ","
  ) %>% type_convert(),

  ## create a file with start and end dates for each particpant for extracting ECG data
  pheno_dates = pheno_raw %>% group_by(GEN) %>% mutate(Date = as.Date(Date, format = '%d.%m.%y')) %>% mutate(min_date = min(Date, na.rm =
                                                                                                                              T)) %>% mutate(max_date = max(Date, na.rm = T)) %>% dplyr::select(GEN, min_date, max_date) %>% unique(),
  pheno_dates_out =  write.table(
    pheno_dates,
    file_out("output/PSYMETAB_inclusion_dates.txt"),
    row.names = F,
    quote = F,
    col.names = T
  ),

  caffeine_raw = target(read_excel(file_in(!!caffeine_file), sheet = 1) %>% type_convert(),
                        hpc = FALSE),
  # found a mistake that the age of WIFRYNSK was wrong at one instance.

  caffeine_munge = munge_caffeine(caffeine_raw),
  psymet_bgen_sample = target({
    readr::read_delim(
      file_in(
        !!paste0(
          "analysis/QC/15_final_processing/FULL/",
          study_name,
          ".FULL.sample"
        )
      ),
      col_types = cols(.default = col_character()),
      delim = " "
    ) %>% type_convert()
  }),
  psymet_bgen_nosex_out = write.table(
    psymet_bgen_sample[, 1:3],
    file_out(
      !!paste0(
        "analysis/QC/15_final_processing/FULL/",
        study_name,
        ".FULL_nosex.sample"
      )
    ),
    row.names = F,
    quote = F,
    col.names = T
  ),
  pc_raw = read_pcs(file_in(!!pc_dir),!!study_name,!!eths) %>% as_tibble(),

  pheno_munge = munge_pheno(
    pheno_raw,
    !!baseline_vars,
    !!leeway_time,
    caffeine_munge,
    !!follow_up_6mo,
    !!follow_up_3mo,
    !!follow_up_1mo,
    !!interaction_outcome
  ),
  # pheno_munge %>% count(GEN) %>% filter(n!=1) ## NO DUPLICATES !
  #pheno_munge = munge_pheno(pheno_raw, baseline_vars, leeway_time, caffeine_munge, follow_up_limit)
  pheno_merge = merge_pheno_caffeine(pheno_raw, caffeine_munge,!!anonymization_error),
  # in the end, we don't use this dataset, but if you want to merge by caffeine with the appropriate date use this data.

  pheno_baseline = inner_join(
    pc_raw %>% mutate_at("GPCR", as.character),
    pheno_munge %>% mutate_at("GEN", as.character),
    by = c("GPCR" = "GEN")
  ) %>%
    replace_na(list(sex = 'NONE')),
  pheno_eths_out = write.table(
    pheno_baseline %>% tidyr::separate(FID, c("COUNT", "GPCR"), "_") %>%
      dplyr::select(COUNT, GPCR, eth),
    file_out(
      !!paste0(
        "data/processed/phenotype_data/",
        study_name,
        "_inferred_eths.txt"
      )
    ),
    row.names = F,
    quote = F,
    col.names = T
  ),

  pc_raw_out = write.table(
    pc_raw,
    file_out(
      !!paste0("data/processed/phenotype_data/", study_name, "_PCs.txt")
    ),
    row.names = F,
    quote = F,
    col.names = T
  ),

  pheno_drug_combinations = define_pheno_drug_combos(!!drug_classes,!!interaction_vars),
  pheno_followup = target(
    munge_pheno_follow(
      pheno_baseline,
      !!test_drugs,
      pheno_drug_combinations$class,
      pheno_drug_combinations$pheno,
      !!low_inducers,
      !!high_inducers
    ),
    dynamic = map(pheno_drug_combinations)
  ),
  #names(pheno_followup) is the names defined in `test_drugs`: tibble
  # follow-up data is smaller than baseline data because participants were filtered if who have they had taken drug but it wasn't followed (because we don't know when they would have taken the drug)

  GWAS_input = create_GWAS_pheno(
    pheno_baseline,
    pheno_followup,
    !!caffeine_vars,
    !!test_drugs,!!high_inducers,
    !!med_inducers,
    !!low_inducers,
    !!outcomes
  ),


  ## linear model GWAS:
  #linear_pheno = pheno_baseline %>% dplyr::select(matches('^FID$|^IID$|_start_Drug_1')),
  #linear_covar = pheno_baseline %>% dplyr::select(matches('^FID$|^IID$|^sex$|^Age_Drug_1$|Age_sq_Drug_1$|^PC[0-9]$|^PC[1][0-9]$|^PC20')),
  #pheno_followup_split = split_followup(pheno_followup, !!test_drugs),
  #followup_data_write = if(!is.null(create_pheno_folder)){ write_followup_data(pheno_followup_split)},
  out_pheno_munge =  target({
    write.table(
      pheno_munge,
      file_out(
        "data/processed/phenotype_data/GWAS_input/pheno_munge.txt"
      ),
      row.names = F,
      quote = F,
      col.names = T
    )
  }),

  out_linear_pheno =  target({
    GWAS_input
    write.table(
      GWAS_input$full_pheno,
      file_out(
        "data/processed/phenotype_data/GWAS_input/pheno_input.txt"
      ),
      row.names = F,
      quote = F,
      col.names = T
    )
  }),
  out_linear_covar =  target({
    GWAS_input
    write.table(
      GWAS_input$full_covar,
      file_out(
        "data/processed/phenotype_data/GWAS_input/covar_input.txt"
      ),
      row.names = F,
      quote = F,
      col.names = T
    )
  }),

  baseline_gwas_info = target(
    define_baseline_models(
      GWAS_input,
      !!baseline_vars,
      !!drug_classes,
      !!caffeine_vars,
      !!interaction_outcome
    ),
    hpc = FALSE
  ),
  interaction_gwas_info = target(
    define_interaction_models(
      GWAS_input,
      !!drug_classes,
      !!interaction_outcome,
      !!baseline_vars
    ),
    hpc = FALSE
  ),
  subgroup_gwas_info = target(
    define_subgroup_models(
      GWAS_input,
      !!drug_classes,
      !!interaction_outcome,
      !!baseline_vars
    ),
    hpc = FALSE
  ),

  baseline_resid = target(
    residualize_pheno(
      GWAS_input,
      pheno = baseline_gwas_info$pheno,
      covars = baseline_gwas_info$covars,
      outcome_type = "baseline",
      !!eths,
      eth_data = pc_raw
    ),
    dynamic = map(baseline_gwas_info)
  ),

  baseline_resid_data = reduce(baseline_resid, inner_join),

  interaction_resid = target(
    residualize_pheno(
      GWAS_input,
      pheno = interaction_gwas_info$pheno,
      covars = interaction_gwas_info$covars,
      outcome_type = "interaction",
      !!eths,
      eth_data = pc_raw
    ),
    dynamic = map(interaction_gwas_info)
  ),

  interaction_resid_data = reduce(interaction_resid, inner_join),

  subgroup_resid = target(
    residualize_pheno(
      GWAS_input,
      pheno = subgroup_gwas_info$pheno,
      covars = subgroup_gwas_info$covars,
      outcome_type = "subgroup",
      !!eths,
      eth_data = pc_raw
    ),
    dynamic = map(subgroup_gwas_info)
  ),
  subgroup_resid_data = reduce(subgroup_resid, inner_join),

  drug_classes_tibble = tibble(class = !!drug_classes),

  out_baseline_resid =  target({
    write.table(
      baseline_resid_data,
      file_out(
        "data/processed/phenotype_data/GWAS_input/baseline_input_resid.txt"
      ),
      row.names = F,
      quote = F,
      col.names = T
    )
  }),

  out_interacation_resid = target(
    write_interaction_resid(interaction_resid_data, drug_classes_tibble$class),
    dynamic = map(drug_classes_tibble),
    format = "file"
  ),

  out_subgroup_resid = target(
    write_subgroup_resid(subgroup_resid_data, drug_classes_tibble$class),
    dynamic = map(drug_classes_tibble),
    format = "file"
  ),

)

# vis_drake_graph(analysis_prep)

# based on Zoltan's suggestion, divide PSYMETAB into subgroups based on different drugs such that each participant is only in one subgroup.


analysis_prep_case_only <- drake_plan(
  interaction_vars_tibble = tibble(phenotype = !!interaction_vars),

  pheno_case_only = target(
    munge_pheno_case_only(
      pheno_baseline,
      !!drug_prioritization,
      interaction_vars_tibble$phenotype,
      !!low_inducers,
      !!high_inducers
    ),
    dynamic = map(interaction_vars_tibble)
  ),

  GWAS_input_case_only = create_GWAS_pheno_case_only(
    pheno_baseline,
    pheno_case_only,
    !!caffeine_vars,
    !!test_drugs,!!high_inducers,
    !!med_inducers,
    !!low_inducers,
    !!outcomes
  ),

  out_case_only_pheno =  target({
    GWAS_input
    write.table(
      GWAS_input_case_only$full_pheno,
      file_out(
        "data/processed/phenotype_data/GWAS_input/case_only_pheno_input.txt"
      ),
      row.names = F,
      quote = F,
      col.names = T
    )
  }),
  out_case_only_covar =  target({
    GWAS_input
    write.table(
      GWAS_input_case_only$full_covar,
      file_out(
        "data/processed/phenotype_data/GWAS_input/case_only_covar_input.txt"
      ),
      row.names = F,
      quote = F,
      col.names = T
    )
  }),


  case_only_gwas_info = target(
    define_case_only_models(
      GWAS_input_case_only,
      !!drug_prioritization,
      !!interaction_outcome,
      !!baseline_vars
    ),
    hpc = FALSE
  ),

  case_only_resid = target(
    residualize_pheno_case_only(
      GWAS_input_case_only,
      pheno = case_only_gwas_info$pheno,
      covars = case_only_gwas_info$covars,
      outcome_type = "case_only",
      !!eths,
      eth_data = pc_raw
    ),
    dynamic = map(case_only_gwas_info)
  ),

  case_only_resid_data = reduce(case_only_resid, inner_join),

  out_case_only_resid = target({
    write.table(
      case_only_resid_data,
      file_out(
        "data/processed/phenotype_data/GWAS_input/case_only_input_resid.txt"
      ),
      row.names = F,
      quote = F,
      col.names = T
    )
  }),


)

#analysis_prep_merge <- bind_plans(analysis_prep, analysis_prep_case_only)

## Case only resid data has number of columns equal to length(drug_prioritization) * length(outcomes),
## where each drug in `drug_prioritization` corresponds to length(outcomes) columns (i.e. BMI slope, 6mo, 3mo, weighted slope, etc.).
## The database is split into length(drug_prioritization) independent groups.

## Test the independence with the following code:

# test <- reduce(case_only_resid[seq(1, 144, 16)], inner_join)
#
#
# result <- numeric()
# for (i in 3:dim(test)[2]){
#
#   non_na_i <- which(!is.na(test[,i]))
#   result <- c(result, non_na_i)
# }
#
#
# length(unique(result))
# length(result)
# dim(test)[2]

# Test if they all equal each other:
# all(sapply(list(length(unique(result)), length(result), dim(test)[1]), function(x) x == dim(test)[1]))



##### Run GWAS in PSYMETAB ---------

init_analysis <- drake_plan(

  GWAS_input_analysis = target(list(full_pheno = read_table2(file_in("data/processed/phenotype_data/GWAS_input/pheno_input.txt")),
                                    full_covar = read_table2(file_in("data/processed/phenotype_data/GWAS_input/covar_input.txt"))),
    hpc = FALSE), # this should be identical to GWAS_input above, but needs to be a different name

  GWAS_input_case_only_analysis = target(list(full_pheno = read_table2(file_in("data/processed/phenotype_data/GWAS_input/case_only_pheno_input.txt")),
                                    full_covar = read_table2(file_in("data/processed/phenotype_data/GWAS_input/case_only_covar_input.txt"))),
    hpc = FALSE), # this should be identical to GWAS_input_case_only above, but needs to be a different name


  # define endpoint, covars and outputs ---------------------------

  baseline_gwas_input = target(define_baseline_models(GWAS_input_analysis, !!baseline_vars, !!drug_classes, !!caffeine_vars, !!interaction_outcome),
    hpc = FALSE),
  interaction_gwas_input = target(define_interaction_inputs(!!drug_classes),
    hpc = FALSE),
  subgroup_gwas_input = target(define_subgroup_inputs(!!drug_classes),
    hpc = FALSE),
  case_only_gwas_input = target(define_case_only_models(GWAS_input_case_only_analysis, !!drug_prioritization, !!interaction_outcome, !!baseline_vars),
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

  case_only_out = target({
    #loadd(baseline_gwas_info
    file_in(!!paste0("analysis/QC/15_final_processing/FULL/", study_name, ".FULL.log"))
    run_gwas(pfile = paste0("analysis/QC/15_final_processing/FULL/", !!study_name, ".FULL"),
               pheno_file = file_in("data/processed/phenotype_data/GWAS_input/case_only_input_resid.txt"),
               threads = 8, eths = !!eths, type="case_only",
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

  # get frequency counts for each GWAS


  linear_freq_out = target({
    #loadd(baseline_gwas_info
    file_in(!!paste0("analysis/QC/15_final_processing/FULL/", study_name, ".FULL.log"))
    run_gwas_freq_counts(pfile = paste0("analysis/QC/15_final_processing/FULL/", !!study_name, ".FULL"),
             pheno_file = file_in("data/processed/phenotype_data/GWAS_input/baseline_input_resid.txt"),
             threads = 8, eths = !!eths,
             eth_sample_file = paste0("analysis/QC/12_ethnicity_admixture/pca/", !!study_name, "_ETH_samples.txt"), ## this is not a real file - "ETH" gets replaced by proper "ETH" in `run_gwas`
             eth_low_maf_file = paste0("analysis/QC/14_mafcheck/", !!study_name, "_ETH_low_maf_snps.txt"), ## this is not a real file - "ETH" gets replaced by proper "ETH" in `run_gwas`
             remove_sample_file = file_in(!!paste0("analysis/QC/11_relatedness/", study_name, "_related_ids.txt")),
             output_dir = (!!plink_output_dir), output = paste0(!!study_name, "_FREQ"))}
  ),

  interaction_freq_out = target({
    check_interaction_input_files
    file_in(!!paste0("analysis/QC/15_final_processing/FULL/", study_name, ".FULL.log"))
    run_gwas_freq_counts(pfile = paste0("analysis/QC/15_final_processing/FULL/", !!study_name, ".FULL"), pheno_file = interaction_gwas_input$pheno_file,
             covar_file = interaction_gwas_input$covar_file,type = "interaction",
             threads = 8, covar_names = interaction_gwas_input$covars, parameters = interaction_gwas_input$parameters, output_suffix = interaction_gwas_input$output_suffix,
             eths = !!eths, eth_sample_file = paste0("analysis/QC/12_ethnicity_admixture/pca/", !!study_name, "_ETH_samples.txt"),  ## this is not a real file - "ETH" gets replaced by proper "ETH" in `run_gwas`
             eth_low_maf_file = paste0("analysis/QC/14_mafcheck/", !!study_name, "_ETH_low_maf_snps.txt"), ## this is not a real file - "ETH" gets replaced by proper "ETH" in `run_gwas`
             remove_sample_file = file_in(!!paste0("analysis/QC/11_relatedness/", study_name, "_related_ids.txt")),
             output_dir = (!!plink_output_dir), output = paste0(!!study_name, "_FREQ"))},
    dynamic = map(interaction_gwas_input)
  ),

  subgroup_freq_out = target({
    check_subgroup_input_files
    file_in(!!paste0("analysis/QC/15_final_processing/FULL/", study_name, ".FULL.log"))
    run_gwas_freq_counts(pfile = paste0("analysis/QC/15_final_processing/FULL/", !!study_name, ".FULL"), pheno_file = subgroup_gwas_input$pheno_file,
                type = "subgroup", subgroup_var = subgroup_gwas_input$subgroup,
                threads = 8, output_suffix = subgroup_gwas_input$output_suffix,
                eths = !!eths, eth_sample_file = paste0("analysis/QC/12_ethnicity_admixture/pca/", !!study_name, "_ETH_samples.txt"),  ## this is not a real file - "ETH" gets replaced by proper "ETH" in `run_gwas`
                eth_low_maf_file = paste0("analysis/QC/14_mafcheck/", !!study_name, "_ETH_low_maf_snps.txt"), ## this is not a real file - "ETH" gets replaced by proper "ETH" in `run_gwas`
                remove_sample_file = file_in(!!paste0("analysis/QC/11_relatedness/", study_name, "_related_ids.txt")),
                output_dir = (!!plink_output_dir), output = paste0(!!study_name, "_FREQ"))},
    dynamic = map(subgroup_gwas_input)
  ),

  case_only_freq_out = target({
    #loadd(baseline_gwas_info
    file_in(!!paste0("analysis/QC/15_final_processing/FULL/", study_name, ".FULL.log"))
    run_gwas_freq_counts(pfile = paste0("analysis/QC/15_final_processing/FULL/", !!study_name, ".FULL"),
               pheno_file = file_in("data/processed/phenotype_data/GWAS_input/case_only_input_resid.txt"),
               threads = 8, eths = !!eths, type="case_only",
               eth_sample_file = paste0("analysis/QC/12_ethnicity_admixture/pca/", !!study_name, "_ETH_samples.txt"), ## this is not a real file - "ETH" gets replaced by proper "ETH" in `run_gwas`
               eth_low_maf_file = paste0("analysis/QC/14_mafcheck/", !!study_name, "_ETH_low_maf_snps.txt"), ## this is not a real file - "ETH" gets replaced by proper "ETH" in `run_gwas`
               remove_sample_file = file_in(!!paste0("analysis/QC/11_relatedness/", study_name, "_related_ids.txt")),
               output_dir = (!!plink_output_dir), output = paste0(!!study_name, "_FREQ"))}
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

  case_only_meta_eths_out = target({
    case_only_out
    meta_case_only_eths(output = !!study_name, output_suffix = "", eths = !!eths,
      pheno = case_only_gwas_input$output_suffix, threads = 8, output_dir = (!!plink_output_dir), type = "case_only")},
    dynamic = map(case_only_gwas_input)),

  interaction_outcome_tibble = tibble(outcome = !!interaction_outcome),

  case_only_meta_drugs_out = target({ ##zip during this run
    case_only_out
    case_only_meta_eths_out
    meta_case_only_drugs(output = !!study_name, output_suffix = "", eths = !!eths, drug_groups = !!drug_prioritization,
      outcome = interaction_outcome_tibble$outcome, threads = 8, output_dir = (!!plink_output_dir), type = "case_only")},
    dynamic = map(interaction_outcome_tibble)),

)

# vis_drake_graph(init_analysis)


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

  GWAS_input_case_only_process = target(list(full_pheno = read_table2(file_in("data/processed/phenotype_data/GWAS_input/case_only_pheno_input.txt")),
                                    full_covar = read_table2(file_in("data/processed/phenotype_data/GWAS_input/case_only_covar_input.txt"))),
    hpc = FALSE), # this should be identical to GWAS_input_case_only above, but needs to be a different name

  baseline_gwas_process = target(define_baseline_models(GWAS_input_process, !!baseline_vars, !!drug_classes, !!caffeine_vars, !!interaction_outcome),
    hpc = FALSE),
  interaction_gwas_process = target(define_interaction_inputs(!!drug_classes),
    hpc = FALSE),
  subgroup_gwas_process = target(define_subgroup_inputs(!!drug_classes),
    hpc = FALSE),
  case_only_gwas_process = target(define_case_only_models(GWAS_input_case_only_process, !!drug_prioritization, !!interaction_outcome, !!baseline_vars),
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
  case_only_gwas_files = target({
    define_case_only_files(case_only_gwas_process, output = !!study_name, eths = !!eths,
    output_dir = "analysis/GWAS", type = "case_only", pheno_list = !!interaction_outcome)},
    hpc = FALSE),

  subgroup_freq_files = subgroup_gwas_files %>% mutate_if(is.character, ~str_replace(., "PSYMETAB_GWAS", "PSYMETAB_GWAS_FREQ")) %>% filter(eth!="META"),
  interaction_freq_files = interaction_gwas_files %>% mutate_if(is.character, ~str_replace(., "PSYMETAB_GWAS", "PSYMETAB_GWAS_FREQ")) %>% filter(eth!="META"),
  baseline_freq_files = baseline_gwas_files %>% mutate_if(is.character, ~str_replace(., "PSYMETAB_GWAS", "PSYMETAB_GWAS_FREQ")) %>% filter(eth!="META"),

  case_only_freq_files = case_only_gwas_files %>% mutate_if(is.character, ~str_replace(., "PSYMETAB_GWAS", "PSYMETAB_GWAS_FREQ")) %>% filter(eth!="META_DRUGS") %>% filter(eth!="META"),

  check_baseline_files = target({
    baseline_gwas_files$log_file},
    dynamic=map(baseline_gwas_files), format = "file"),

  check_interaction_files = target({
    interaction_gwas_files$log_file},
    dynamic=map(interaction_gwas_files), format = "file"),

  check_subgroup_files = target({
    subgroup_gwas_files$log_file},
    dynamic=map(subgroup_gwas_files), format = "file"),

  check_case_only_files = target({
    case_only_gwas_files$log_file},
    dynamic=map(case_only_gwas_files), format = "file"),

  ## save frequency statistics for subgroup and case only analyses

  check_subgroup_freq_files = target({
    subgroup_freq_files$log_file},
    dynamic=map(subgroup_freq_files), format = "file"),

  check_case_only_freq_files = target({
    case_only_freq_files$log_file},
    dynamic=map(case_only_freq_files), format = "file"),

  check_interaction_freq_files = target({
    interaction_freq_files$log_file},
    dynamic=map(interaction_freq_files), format = "file"),

  check_baseline_freq_files = target({
    baseline_freq_files$log_file},
    dynamic=map(baseline_freq_files), format = "file"),

  # write formatted GWAS file ---------------------------

  gwas_reference_file_out = create_gwas_reference_file(freq_file = "analysis/QC/15_final_processing/CEU/PSYMETAB_GWAS.CEU.afreq", info_file = "analysis/QC/15_final_processing/PSYMETAB_GWAS.info",
                        output_file = file_out("analysis/GWAS/GWAS_reference_file.txt")),

  process_baseline_gwas = target({
    check_baseline_files
    gwas_reference_file_out
    process_gwas(eth = baseline_gwas_files$eth, pheno=baseline_gwas_files$pheno, drug=baseline_gwas_files$drug, file=baseline_gwas_files$plink_file,
      output = !!study_name, output_dir = "analysis/GWAS", type = "full",
      reference_file = file_in("analysis/GWAS/GWAS_reference_file.txt"), out_file = baseline_gwas_files$processed_file)
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
    gwas_reference_file_out
    process_gwas(eth = interaction_gwas_files$eth, pheno=interaction_gwas_files$pheno, drug=interaction_gwas_files$drug, file=interaction_gwas_files$plink_file,
      output = !!study_name, output_dir = "analysis/GWAS", type = "interaction",
      reference_file = file_in("analysis/GWAS/GWAS_reference_file.txt"), out_file = interaction_gwas_files$processed_file)
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
    gwas_reference_file_out
    process_gwas(eth = subgroup_gwas_files$eth, pheno=subgroup_gwas_files$pheno, drug=subgroup_gwas_files$drug, file=subgroup_gwas_files$plink_file,
      output = !!study_name, output_dir = "analysis/GWAS", type = "subgroup",
      reference_file = file_in("analysis/GWAS/GWAS_reference_file.txt"), out_file = subgroup_gwas_files$processed_file)
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

  # write formatted GWAS file ---------------------------

  case_only_gwas_files_mod = case_only_gwas_files %>% mutate(eth_mod = case_when(eth == "META_DRUGS" ~ "META", TRUE ~ eth)),

  process_case_only_gwas = target({
    check_case_only_files
    gwas_reference_file_out
    process_gwas(eth = case_only_gwas_files_mod$eth_mod, pheno=case_only_gwas_files_mod$pheno, drug=case_only_gwas_files_mod$drug, file=case_only_gwas_files_mod$plink_file,
      output = !!study_name, output_dir = "analysis/GWAS", type = "case_only",
      reference_file = file_in("analysis/GWAS/GWAS_reference_file.txt"), out_file = case_only_gwas_files_mod$processed_file)
    case_only_gwas_files_mod$processed_file
    },
    dynamic = map(case_only_gwas_files_mod), format = "file"
  ),

  case_only_gwas_figures_input = target({
    process_case_only_gwas
    gwas_figures_input(eth = case_only_gwas_files_mod$eth_mod, pheno=case_only_gwas_files_mod$pheno, drug=case_only_gwas_files_mod$drug, file=case_only_gwas_files_mod$plink_file,
      output = !!study_name, output_dir = "analysis/GWAS", type = "case_only",
      info_file = paste0("analysis/QC/15_final_processing/", !!study_name, ".info"), out_file = case_only_gwas_files_mod$processed_file)},
    dynamic = map(case_only_gwas_files_mod), hpc = FALSE),

  # create manhattan and qq figures ---------------------------

  case_only_gwas_figures = target({
    process_case_only_gwas
    create_figures(
      case_only_gwas_figures_input$joint_file, case_only_gwas_figures_input$manhattan_file_name,
      case_only_gwas_figures_input$qq_file_name, case_only_gwas_figures_input$title, !!gw_sig)

      c(case_only_gwas_figures_input$manhattan_file_name, case_only_gwas_figures_input$qq_file_name)
    },
    dynamic = map(case_only_gwas_figures_input) ,format = "file"
  ),


  # extract nominally significant results (results in a list of nominally significant tibbles)---------------------------

  baseline_gwas_sig = target({
    check_baseline_freq_files
    process_baseline_gwas
    gw_sig_extract(gwas_file = baseline_gwas_files$processed_file, gz_plink_file = baseline_gwas_files$plink_file, !!gw_sig_nominal, eth = baseline_gwas_files$eth,
      pheno = baseline_gwas_files$pheno, drug_class = NA, drug = baseline_gwas_files$drug)
    },
    dynamic = map(baseline_gwas_files)
  ),

  interaction_gwas_sig = target({
    check_interaction_freq_files
    process_interaction_gwas
    gw_sig_extract(gwas_file = interaction_gwas_files$processed_file, gz_plink_file = interaction_gwas_files$plink_file, !!gw_sig_nominal, eth = interaction_gwas_files$eth,
      pheno = interaction_gwas_files$pheno, drug_class = interaction_gwas_files$drug_class, drug = interaction_gwas_files$drug, interaction = T)
    },
    dynamic = map(interaction_gwas_files)
  ),

  subgroup_gwas_sig = target({
    check_subgroup_freq_files
    process_subgroup_gwas
    gw_sig_extract(gwas_file = subgroup_gwas_files$processed_file, gz_plink_file = subgroup_gwas_files$plink_file, !!gw_sig_nominal, eth = subgroup_gwas_files$eth,
      pheno = subgroup_gwas_files$pheno, drug_class = subgroup_gwas_files$drug_class, drug = subgroup_gwas_files$drug)
    },
    dynamic = map(subgroup_gwas_files)
  ),

  case_only_gwas_sig = target({
    check_case_only_freq_files
    process_case_only_gwas
    gw_sig_extract(gwas_file = case_only_gwas_files_mod$processed_file, gz_plink_file = case_only_gwas_files_mod$plink_file, !!gw_sig_nominal, eth = case_only_gwas_files_mod$eth,
      pheno = case_only_gwas_files_mod$pheno, drug_class = NA, drug = case_only_gwas_files_mod$drug)
    },
    dynamic = map(case_only_gwas_files_mod)
  ),

  case_only_celine_extract = target({
    #chr = 10 #BP_start=61815632 # BP_stop=61820501
    #BP+- 500K around 61834898
    check_case_only_freq_files
    process_case_only_gwas
    region_extract(gwas_file = case_only_gwas_files_mod$processed_file, gz_plink_file = case_only_gwas_files_mod$plink_file, chr = 10, bp_start = 61834898-500000, bp_stop = 61834898 + 500000, eth = case_only_gwas_files_mod$eth,
                   pheno = case_only_gwas_files_mod$pheno, drug_class = NA, drug = case_only_gwas_files_mod$drug)
  },
  dynamic = map(case_only_gwas_files_mod)
  ),


  baseline_gwas_sig_clean = target(clean_gwas_summary(baseline_gwas_sig), dynamic = map(baseline_gwas_sig)),
  interaction_gwas_sig_clean = target(clean_gwas_summary(interaction_gwas_sig), dynamic = map(interaction_gwas_sig)),
  subgroup_gwas_sig_clean = target(clean_gwas_summary(subgroup_gwas_sig), dynamic = map(subgroup_gwas_sig)),

  case_only_gwas_sig_clean = target(clean_case_only_summary(case_only_gwas_sig),
    dynamic = map(case_only_gwas_sig)),

  case_only_celine_extract_clean = target(clean_case_only_summary(case_only_celine_extract),
                                    dynamic = map(case_only_celine_extract)),

  # combined nominally significant results into one tibble ---------------------------

  baseline_gwas_nom_com = bind_rows(baseline_gwas_sig_clean, .id = "pheno_name"),
  interaction_gwas_nom_com = bind_rows(interaction_gwas_sig_clean, .id = "pheno_name"),
  subgroup_gwas_nom_com = bind_rows(subgroup_gwas_sig_clean, .id = "pheno_name"),
  case_only_gwas_nom_com = bind_rows(case_only_gwas_sig_clean, .id = "pheno_name"),

  case_only_celine_extract_com = bind_rows(case_only_celine_extract_clean, .id = "pheno_name"),
  # filter above nominally significant results into GW-significant results  ---------------------------

  baseline_gwas_sig_com = baseline_gwas_nom_com %>% filter(P < !!gw_sig),
  interaction_gwas_sig_com = interaction_gwas_nom_com %>% filter(P < !!gw_sig),
  subgroup_gwas_sig_com = subgroup_gwas_nom_com %>% filter(P < !!gw_sig),
  case_only_gwas_sig_com = case_only_gwas_nom_com %>% filter(P < !!gw_sig),

  # count significant results ---------------------------

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

  # filter into only categories of interest ---------------------------

  baseline_interest = baseline_gwas_sig_com %>% filter(eth=="CEU") %>% filter(!grepl("sensitivity", pheno)) %>% dplyr::filter(outcome!="slope_weight_6mo") %>% dplyr::filter(outcome!="slope_weight"),
  subgroup_interest = subgroup_gwas_sig_com %>% filter(drug=="Drug") %>% filter(eth=="CEU") %>% filter(!grepl("sensitivity", drug_class)) %>% filter(grepl("BMI", pheno)) %>% dplyr::filter(outcome!="slope_weight_6mo") %>% dplyr::filter(outcome!="slope_weight"),
  case_only_interest = case_only_gwas_sig_com %>% filter(group=="META_DRUGS") %>% filter(eth=="CEU") %>% filter(variable=="BMI") %>% dplyr::filter(outcome!="slope_weight_6mo") %>% dplyr::filter(outcome!="slope_weight"),

  case_only_celine_interest = case_only_celine_extract_com %>% filter(group=="META_DRUGS") %>% filter(eth=="CEU") %>% filter(variable=="BMI") %>% dplyr::filter(outcome!="slope_weight_6mo") %>% dplyr::filter(outcome!="slope_weight"),

  # extract model specific results for case-only and pull Q-statistic

  case_only_interest_summarize = target(summarize_case_only_meta(case_only_interest$SNP, case_only_interest$eth, case_only_interest$variable, case_only_interest$outcome,
                                                                 !!study_name, !!drug_prioritization, directory = "analysis/GWAS/case_only"),
                                        dynamic=map(case_only_interest)
  ),

  #case_only_celine_interest_summarize = target(summarize_case_only_meta(case_only_celine_interest$SNP, case_only_celine_interest$eth, case_only_celine_interest$variable, case_only_celine_interest$outcome,
  #                                                                      !!study_name, !!drug_prioritization, directory = "analysis/GWAS/case_only"),
  #                                             dynamic=map(case_only_celine_interest)
  #),


  subgroup_interest_prune = crude_prune(subgroup_interest, grouping_var = "pheno_name"), # prune by position

  ## Define list of SNPs that are significant in case_only and subgroup models and extract these SNPs from files of interest

  SNPs_interest = c(case_only_interest$SNP, subgroup_interest$SNP),

  SNPs_interest_chr = tibble(chr = c(case_only_interest$CHR, subgroup_interest$CHR), rsid = SNPs_interest),
  consort_files = list.files(!!paste0(consortia_dir, "/formatted"), full.names = TRUE),

  check_consort_files = target({
    consort_files},
    dynamic=map(consort_files), format = "file"),

  consort_extract = target({
    check_consort_files
    consort_lookup(SNPs_interest, consort_files)
  }, dynamic = map(consort_files)),

  consort_extract_agg = rbindlist(consort_extract),

  ## have the lists be sorted SNPs instead of by consort

  consort_extract_snpwise = target({
    temp <- consort_extract_agg  %>% group_by(SNP)
    out <- temp %>% group_split() %>% set_names(unlist(group_keys(temp)))
  }
  ),

  ## check if significant snps (`SNPs_interest`) have consitent findings in similar GWAS (i.e. 3 mo, 6 mo, full)

  case_only_extract = target({
    check_case_only_freq_files
    process_case_only_gwas
    snp_psymetab_extract(gwas_file = case_only_gwas_files_mod$processed_file, gz_plink_file = case_only_gwas_files_mod$plink_file, SNPs_interest, eth = case_only_gwas_files_mod$eth,
                         pheno = case_only_gwas_files_mod$pheno, drug_class = NA, drug = case_only_gwas_files_mod$drug)
  },
  dynamic = map(case_only_gwas_files_mod)
  ),

  subgroup_extract = target({
    check_subgroup_freq_files
    process_subgroup_gwas
    snp_psymetab_extract(gwas_file = subgroup_gwas_files$processed_file, gz_plink_file = subgroup_gwas_files$plink_file, SNPs_interest, eth = subgroup_gwas_files$eth,
                         pheno = subgroup_gwas_files$pheno, drug_class = subgroup_gwas_files$drug_class, drug = subgroup_gwas_files$drug)
  },
  dynamic = map(subgroup_gwas_files)
  ),

  case_only_extract_agg = rbindlist(case_only_extract, fill=T),
  case_only_extract_snpwise = target({
    temp <- case_only_extract_agg  %>% group_by(SNP)
    out <- temp %>% group_split() %>% set_names(unlist(group_keys(temp)))
  }
  ),

  case_only_extract_snpwise_clean = target({
    clean_case_only_snpwise(case_only_extract_snpwise)
  }, dynamic = map(case_only_extract_snpwise)
  ),

  subgroup_extract_agg = rbindlist(subgroup_extract,  fill=T),

  ## have the lists be sorted SNPs instead of by consort
  subgroup_extract_snpwise = target({
    temp <- subgroup_extract_agg  %>% group_by(SNP)
    out <- temp %>% group_split() %>% set_names(unlist(group_keys(temp)))
  }
  ),

  write_baseline_interest = write.csv(baseline_interest,  file_out("output/PSYMETAB_GWAS_baseline_CEU_result.csv"),row.names = F),
  write_subgroup_interest = write.csv(subgroup_interest,  file_out("output/PSYMETAB_GWAS_subgroup_CEU_result.csv"),row.names = F),
  write_case_only_interest = write.csv(case_only_interest_summarize,  file_out("output/PSYMETAB_GWAS_case_only_CEU_meta_result.csv"),row.names = F),
  #write_case_only_celine_interest = write.csv(case_only_celine_interest_summarize,  file_out("output/PSYMETAB_GWAS_celine_extract.csv"),row.names = F),

)

# vis_drake_graph(process_init)

ukbb_analysis <- drake_plan(
    ukbb_org = target(ukb_df_mod("ukb21067", path = !!paste0(UKBB_dir, "/org")), hpc = FALSE),
    ukb_key = target(ukb_df_field("ukb21067", path = !!paste0(UKBB_dir, "/org")), hpc = FALSE),
    sqc = target(fread(file_in(!!paste0(UKBB_dir, ukbb_sqc_file)), header=F, data.table=F), hpc = FALSE),
    fam = target(fread(file_in(!!paste0(UKBB_dir, ukbb_fam_file)), header=F,data.table=F), hpc = FALSE),
    relatives = target(read.table(file_in(!!paste0(UKBB_dir, ukbb_relatives_file)), header=T), hpc = FALSE),
    exclusion_list = target(fread(file_in(!!paste0(UKBB_dir, ukbb_exclusion_file)), data.table=F), hpc = FALSE),
    ukbb_bgen_sample = target(fread(file_in(!!paste0(UKBB_dir, ukbb_sample_file)), header=T, data.table=F), hpc = FALSE),

    qc_data = target(ukb_gen_sqc_names(sqc), hpc = FALSE),
    sqc_munge = target(munge_sqc(sqc,fam), hpc = FALSE),
    british_subset = target(get_british_ids(qc_data, fam), hpc = FALSE),

    ukbb_med_codes = read_tsv(file_in(!!ukbb_med_codes_file)),
    medication_columns = ukb_key %>% filter(field.showcase == "20003") %>% filter(str_extract(field.html, "(?<=.[-]).")==0) %>% pull(col.name),

    psy_ukbb_drugs = read_excel(file_in(!!psychiatric_ukbb_drugs_file), sheet=1),
    psy_ukbb_drugs2 = read_excel(file_in(!!psychiatric_ukbb_drugs_file), sheet=2) %>% # this list of drugs were not in the original analysis but can be included to increase sample size
      mutate_at(vars(Coding), ~as.character(.)),
    psy_ukbb_drugs3 = read_excel(file_in(!!psychiatric_ukbb_drugs_file), sheet=3) %>% # this list of drugs were not in the original analysis but can be included to increase sample size
      mutate(risk = 0) %>% rename(molecule = Molecule) %>%  mutate_at(vars(Coding), ~as.character(.)),
    psy_ukbb_drugs_all = dplyr::bind_rows(psy_ukbb_drugs, psy_ukbb_drugs2, psy_ukbb_drugs3) %>%
      mutate(molecule = case_when(molecule == "Sulpiride (N05AL01)" ~ "amisulpride",
                                TRUE ~ molecule
                             )),

    drug_codes_summary = psy_ukbb_drugs_all %>% group_by(molecule) %>% summarize(codes= list(Coding)), #group codes together by molecule
    drug_codes_interest = drug_codes_summary %>%
      filter(molecule %in% c("olanzapine", "clozapine", "valproate", "quetiapine", "risperidone", "aripiprazole", "mirtazapine", "amisulpride")),

    ukbb_munge_bmi_slope = munge_ukbb_bmi_slope(ukbb_org, ukb_key, !!date_followup, !!bmi_var),
    ukbb_filter_bmi_slope = target(filter_ukbb(ukbb_munge_bmi_slope, relatives, exclusion_list, british_subset),
      dynamic = map(ukbb_munge_bmi_slope)),

    ukbb_munge_drug_users = munge_ukbb_drug_users(ukbb_org, drug_codes_interest, medication_columns),
    ## get counts in UKBB for each molecule
    ukbb_drug_counts = count_ukbb_drug_users(ukbb_munge_drug_users, drug_codes_interest),

    ukbb_munge_drug_users_bmi = munge_ukbb_drug_users_bmi(ukbb_org, ukb_key, ukbb_munge_drug_users, !!bmi_var),
    ukbb_filter_drug_users_bmi = target(filter_ukbb(ukbb_munge_drug_users_bmi, relatives, exclusion_list, british_subset),
      dynamic = map(ukbb_munge_drug_users_bmi)),


    ukbb_filter_all = c(ukbb_filter_bmi_slope, ukbb_filter_drug_users_bmi),

    # ### 28/05/2021: UKBB GWAS - NEEDS TO BE RE-RUN ADJUSTING FOR BASELINE BMI BUT TAKES SOME TIME

    # residualize and ivt
    # ukbb_resid_list = target(resid_ukbb(ukbb_filter_all, ukbb_org, ukb_key, sqc_munge, !!bmi_var, !!sex_var, !!age_var),
    #   dynamic = map(ukbb_filter_all)),
    #
    # ukbb_resid_join = reduce(ukbb_resid_list, full_join),
    # #ukbb_resid_bmi_slope = resid_ukbb(ukbb_filter_bmi_slope, ukb_key, sqc_munge, outcome = "bmi_slope", !!bmi_var, !!sex_var, !!age_var),
    # ukbb_bgen_order = order_bgen(ukbb_bgen_sample, ukbb_resid_join),
    #
    #
    # ukbb_bgen_out = write.table(ukbb_bgen_order, file_out("data/processed/ukbb_data/ukbb_GWAS"), sep=" ", quote=F, row.names=F, col.names = T),
    #
    # #chr_num = tibble(chr = 1:22),
    #
    # v2_snp_list_files = create_UKBB_v2_snp_list(!!UKBB_processed),
    #
    # check_v2_snp_list_files = target({
    #   c(v2_snp_list_files$v2_snp_list_files)},
    #   dynamic = map(v2_snp_list_files), format = "file"),
    #
    # chunk_ukbb = target({
    #   check_v2_snp_list_files
    #   make_ukbb_chunks(v2_snp_list_files$v2_snp_list_files, chunk_size=1e6)
    # }, dynamic = map(v2_snp_list_files)),
    #
    # chunk_ukbb_run = as_tibble(chunk_ukbb),
    #
    #

    # bgenie_out = target({
    #   ukbb_bgen_out
    #   if(!file.exists(paste0("analysis/GWAS/UKBB/chr", chunk_ukbb_run$chr, "_chunk", chunk_ukbb_run$chunk_num, ".out.gz"))){
    #     launch_bgenie(chunk_ukbb_run$chr, phenofile = file_in("data/processed/ukbb_data/ukbb_GWAS"), !!UKBB_dir, chunk_ukbb_run$chr_char, chunk_ukbb_run$start, chunk_ukbb_run$end, chunk_ukbb_run$chunk_num)
    #   }
    #   paste0("analysis/GWAS/UKBB/chr", chunk_ukbb_run$chr, "_chunk", chunk_ukbb_run$chunk_num, ".out.gz")
    # }, dynamic = map(chunk_ukbb_run), format = "file"),
    #
    # # Not all runs will produce a file, because there are no variants in the selected chunk, for example:
    # # launch_bgenie(chunk_ukbb_run$chr[133], phenofile = ("data/processed/ukbb_data/ukbb_GWAS"), UKBB_dir, chunk_ukbb_run$chr_char[133], chunk_ukbb_run$start[133], chunk_ukbb_run$end[133], chunk_ukbb_run$chunk_num[133])
    #
    # bgenie_unzip = target({
    #   bgenie_out
    #   if(!file.exists(paste0("analysis/GWAS/UKBB/chr", chunk_ukbb_run$chr, "_chunk", chunk_ukbb_run$chunk_num, ".out"))){
    #     unzip_bgenie(chunk_ukbb_run$chr, chunk_ukbb_run$chunk_num)
    #   }
    #   paste0("analysis/GWAS/UKBB/chr", chunk_ukbb_run$chr, "_chunk", chunk_ukbb_run$chunk_num, ".out")
    # }, dynamic = map(chunk_ukbb_run), format = "file"),
    #
    # ukbb_gwas = {
    #   bgenie_unzip
    #   process_bgenie(directory="analysis/GWAS/UKBB/", extension=".out", HRC_panel=file_in(!!HRC_panel))
    # },
    # bgenie_process_out = fwrite(ukbb_gwas, file_out("analysis/GWAS/UKBB/UKBB_bgenie_HRC_filtered.txt"), row.names=F, quote=F, sep = "\t"),


)

# vis_drake_graph(ukbb_analysis)

ukbb_control <- drake_plan(

  # bgenie_read = target(
  #   fread(paste0("analysis/GWAS/UKBB/chr", chr_num$chr, ".out"), data.table=F),
  #   #read_bgenie_output(paste0("analysis/GWAS/UKBB/chr", chr_num$chr, ".out")),
  #   dynamic = map(chr_num)
  # ),

  subgroup_files = list.files(path="analysis/GWAS/subgroup/processed", pattern=".txt", full.names=T),

  BMI_slope_sub_files = subgroup_files[grepl("CEU[.]Drug[.]BMI_slope", subgroup_files) & !grepl("weight", subgroup_files) & !grepl( "sensitivity", subgroup_files)],

  BMI_slope_nominal = target(extract_sig_results(BMI_slope_sub_files, !!gw_sig_nominal, eth="CEU"),
    dynamic = map(BMI_slope_sub_files)),


  # BMI_slope_nominal_AF = add_AF(BMI_slope_nominal, AF_file = file_in("analysis/QC/14_mafcheck/CEU/PSYMETAB_GWAS.CEU.afreq")),

  # psy_ukbb_merge = target({
  #   bgenie_unzip
  #   bgenie_read <- fread(paste0("analysis/GWAS/UKBB/chr", chr_num$chr, ".out"), data.table=F)
  #   colnames(bgenie_read) <- paste(colnames(bgenie_read), "UKBB", sep = "_")
  #   psy_ukbb_merge <- left_join(BMI_slope_nominal_AF, bgenie_read, by = c("SNP" = "rsid_UKBB"))
  #   psy_ukbb_merge},
  #   dynamic = map(chr_num)),

  psy_UKBB_merge = target(merge_psy_UKBB(BMI_slope_nominal, UKBB_GWAS_file=file_in("analysis/GWAS/UKBB/UKBB_bgenie_HRC_filtered.txt")),
    dynamic = map(BMI_slope_nominal)),

  psy_UKBB_harmonize = target(harmonize_ukbb_plink_data(psy_UKBB_merge, SNP_col="SNP", REF1_col="A1", ALT1_col="A2", REF2_col="a_1_UKBB", ALT2_col="a_0_UKBB"),
    dynamic = map(psy_UKBB_merge)),

  #psy_ukbb_merge = left_join(BMI_slope_nominal, ukbb_gwas %>% dplyr::select(), by = c("SNP" = "rsid")),
  psy_ukbb_het = target(calc_het(psy_UKBB_harmonize, "SNP", "BETA", "bmi_slope_resid_beta_UKBB", "SE", "bmi_slope_resid_se_UKBB", "match_description"),
    dynamic = map(psy_UKBB_harmonize)),

  psy_ukbb_het_prune = target(prune_psy_ukbb(psy_ukbb_het),
    dynamic = map(psy_ukbb_het)),

  psy_ukbb_het_prune_interest = {
    psy_ukbb_het_prune %>% dplyr::select(-contains("users")) %>% rename(bmi_slope_het_pval = het_pval)},

  write_ukbb_comparison = write.csv(psy_ukbb_het_prune_interest,  file_out("output/PSYMETAB_GWAS_UKBB_controls.csv"),row.names = F),

  # ukbb_comparison = read.csv(file_in("output/PSYMETAB_GWAS_UKBB_comparison.csv")),
  #
  # ukbb_comparison_format = sort_ukbb_comparison(ukbb_comparison, subgroup_GWAS_count, !!interaction_outcome, !!drug_classes),
  # write.csv(ukbb_comparison_format,  file_out("output/PSYMETAB_GWAS_UKBB_comparison2.csv"),row.names = F)

)

ukbb_replication <- drake_plan(

  ukbb_top_snps_chr = target(tibble(chr = psy_ukbb_het$CHR, rsid = psy_ukbb_het$SNP),
    dynamic = map(psy_ukbb_het)),

  #ukbb_geno = target(load_geno(ukbb_bgen_sample, ukbb_top_snps_chr, !!UKBB_dir),
  #  dynamic = map(ukbb_top_snps_chr)),

  psy_UKBB_replicate = target(replicate_psy_UKBB(psy_UKBB_merge, eth="CEU", !!test_drugs),
    dynamic = map(psy_UKBB_merge)),

  psy_UKBB_replicate_sig = target(extract_ukbb_rep_sig(psy_UKBB_replicate, eth="CEU"),
    dynamic = map(psy_UKBB_replicate)),

  # check if BMI correlates with BMI_slope to justify the analysis above
  bmi_slope_corr_data = target(corr_bmi_slope(ukbb_filter_bmi_slope, ukbb_org, ukb_key, !!bmi_var),
    dynamic = map(ukbb_filter_bmi_slope)),

)

ukbb_icd_vs_drug <- drake_plan(

  ukbb_diagnoses = ukb_icd_diagnosis(ukbb_org, ukbb_org$eid, icd.version = 10),

  # ICD10 diagnoses: https://biobank.ctsu.ox.ac.uk/showcase/field.cgi?id=41270
  ukbb_psychiatric_icd10 = ukbb_diagnoses %>% filter(grepl('^F', code)),
  ukbb_psychiatric_icd10_samples = ukbb_psychiatric_icd10 %>% pull(sample) %>% unique(),

  ukbb_psy_subgroups = create_ukbb_psy_subgroups(ukbb_munge_drug_users, ukbb_psychiatric_icd10_samples),
  ukbb_psy_subgroups_bmi_slope = ukbb_subgroups_add_var(ukbb_psy_subgroups, ukbb_filter_bmi_slope$bmi_slope),

  ukbb_psy_subgroups_filter_bmi_slope = target(filter_ukbb(ukbb_psy_subgroups_bmi_slope, relatives, exclusion_list, british_subset),
                                      dynamic = map(ukbb_psy_subgroups_bmi_slope)),


  # residualize and ivt
  ukbb_psy_subgroups_resid_bmi_slope = target(resid_ukbb(ukbb_psy_subgroups_filter_bmi_slope, ukbb_org, ukb_key, sqc_munge, !!bmi_var, !!sex_var, !!age_var),
                           dynamic = map(ukbb_psy_subgroups_filter_bmi_slope)),

  ukbb_geno = load_geno(ukbb_bgen_sample, SNPs_interest_chr, !!UKBB_dir),

  ukbb_psy_subgroups_assoc = target(ukbb_resid_geno_assoc(ukbb_psy_subgroups_resid_bmi_slope, ukbb_geno),
         dynamic = map(ukbb_psy_subgroups_resid_bmi_slope)),

  ukbb_psy_subgroups_assoc_agg = rbindlist(ukbb_psy_subgroups_assoc),


  ukbb_psy_subgroups_assoc_snpwise = target({
    temp <- ukbb_psy_subgroups_assoc_agg  %>% group_by(snp)
    out <- temp %>% group_split() %>% set_names(unlist(group_keys(temp)))
  }
  ),

)





# vis_drake_graph(process_init)

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
