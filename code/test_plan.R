
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
