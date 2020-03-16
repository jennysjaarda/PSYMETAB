
px <- processx:::get_tool("px")

read_fam <- function(create_bed_out)
{
  read.table(paste0(plink_bed_out,".fam"))
}

munge_fam <- function(fam_raw, id_code){
  colnames(fam_raw) <- c("FID", "IID", "fatherID", "motherID", "sex", "pheno")
  fam_id <- merge(fam_raw, id_code, by.x="IID", by.y="randomizedID")
  return(fam_id)
}

munge_qc_pheno <- function(qc_pheno_raw,fam_munge){  ## output of munge_fam
  qc_pheno_raw %>%
    mutate(RecodeSample_ID = str_replace(RecodeSample_ID, "_", "")) %>% mutate_at("RecodeSample_ID", as.character) %>%
    right_join(. , fam_munge %>% mutate_at("sampleID", as.character), by=c("RecodeSample_ID"="sampleID"))
}

find_dups <- function(qc_pheno_munge, fam_munge){

  dups <- qc_pheno_munge[grep(".{8}002", qc_pheno_munge$IID),c("FID", "IID")]
  #dim(dups)  # 15  2
  ## find both dups
  dups2 <- dups %>%
    rowwise() %>%
    mutate(org_ID = unlist(str_split(IID, "(?<=.{8})002"))[1])  %>%
    mutate_at("org_ID", as.character) %>%
    right_join(fam_munge %>% dplyr::select("FID","IID") %>% mutate_at("IID", as.character) , .,by = c("IID" = "org_ID")) %>%
    mutate_all(as.character) %>%
    dplyr::rename(IID.x = IID) %>%
    mutate(id = seq_len(n())) %>%
    gather(v, value, FID.x:IID.y) %>%
    separate(v, c("var", "col"))  %>%
    spread(var, value)  %>%
    dplyr::select("FID","IID")
    return(list(dups=dups, dups_set=dups2))
}

format_sex_file <- function(sex_info){
  sex_info[which(is.na(sex_info$Sexe)),"Sexe"] <- 0
  colnames(sex_info) <- c("FID", "IID", "SEX")
  return(sex_info)
}

format_eth_file <- function(eth_info){
  colnames(eth_info) <- c("FID", "IID", "ETH")
  eth_info$ETH <- as.factor(eth_info$ETH)
  eth_info$ETH2 <- vector(length=nrow(eth_info))
  eth_info$ETH2 <- replace(eth_info$ETH2, eth_info$ETH=="caucasien", 1)
  eth_info$ETH2 <- replace(eth_info$ETH2, eth_info$ETH=="arabe + caucasien", 1)
  eth_info$ETH2 <- replace(eth_info$ETH2, eth_info$ETH=="africain", 2)
  eth_info$ETH2 <- replace(eth_info$ETH2, eth_info$ETH=="amerique du sud", 3)
  ### South America -> Latin
  eth_info$ETH2 <- replace(eth_info$ETH2, eth_info$ETH=="Antilles", 3)
  ### Antilles=West Indes. Should be African?
  eth_info$ETH2 <- replace(eth_info$ETH2, eth_info$ETH=="arabe", 4)
  ### Should be European?
  eth_info$ETH2 <- replace(eth_info$ETH2, eth_info$ETH=="asiatique", 5)
  ### Which part of Asia? South or east?
  eth_info$ETH2 <- replace(eth_info$ETH2, eth_info$ETH=="asiatique + caucasien", 6)
  eth_info$ETH2 <- replace(eth_info$ETH2, eth_info$ETH=="autre", 6)
  eth_info$ETH2 <- replace(eth_info$ETH2, eth_info$ETH=="africain + caucasien", 6)
  eth_info$ETH2 <- replace(eth_info$ETH2, eth_info$ETH=="inconnu", 7)
  eth_info$ETH2 <- replace(eth_info$ETH2, is.na(eth_info$ETH), 7)
  eth_info$ETH2 <- factor(eth_info$ETH2, levels=c(1:7), labels=c("EUROPEAN","AFRICAN", "LATIN", "SOUTH_ASIAN", "EAST_ASIAN", "OTHER", "UNKNOWN"))
  eth_info <- eth_info[,c(1,2,4)]
  colnames(eth_info) <- c("FID", "IID", "ETH")
  return(eth_info)
}

count_imputed_var <- function(){

  imputed_var <- numeric()
  for (chr in 20:22){

    temp <- countLines(paste0("analysis/QC/06_imputation_get/chr", chr, ".info.gz"))[1] - 1
    temp <- cbind(chr, temp)
    imputed_var <- rbind(imputed_var, temp)

  }
  total <- cbind("all", sum(imputed_var[,2]))
  imputed_var <- rbind(imputed_var, total)
  colnames(imputed_var) <- c("Chr", "Num imputed variants")
  imputed_var <- as.data.frame(imputed_var)
  return(imputed_var)

}

rename_meds <- function(x)
{
  out <- rep(1, length(x))
  regimen_num <- 1
  for(i in 1:length(x))
  {
    if(is.na(x[i])){out[i] <- 1}
    if(x[i]=="sensible"){out[i] <- regimen_num}
    if(x[i]=="new_regimen"){
      regimen_num <- regimen_num + 1
      out[i] <- regimen_num
    }
  }
  return(out)
}

bmi_diff <- function(x){n_obs <- length(x); x[n_obs]-x[1]}

check_sex <- function(x){
  out <- NA
  if(("M" %in% x) & ("F" %in% x)) {
    out <- NA
  } else {
    if("M" %in% x) {out <- "M"}
    if("F" %in% x) {out <- "F"}
  }
  return(out)
}

check_drug <- function(PatientRec, Drug){
  if(any(is.na(PatientRec)) |
     any(is.na(Drug))) {out <- "missing Drug or PatientRec"}
  else if(length(unique(na.omit(PatientRec)))==1 &
     length(unique(na.omit(Drug)))==1) {out <- "sensible"} else{
    out <- "non-match"
  }
  return(out)

}
check_height <- function(x){
  if(all(is.na(x))) {out <- NA} else{
    out <- mean(x, na.rm=T)
  }
  return(out)
}

ever_drug <- function(x) {
  if(any(x==1)){out <- 1} else{
    out <- 0
  }
  return(out)
}

munge_pheno <- function(pheno_raw, baseline_vars){
  pheno_raw %>%
    mutate(Date = as.Date(Date, format = '%d.%m.%Y'))  %>% filter(!is.na(Date)) %>% arrange(Date)  %>%
    mutate(AP1 = gsub(" ", "_",AP1)) %>% mutate_at("AP1",as.factor) %>% mutate(AP1 = gsub("_.*$","", AP1)) %>% mutate(AP1 = na_if(AP1, "")) %>% ## merge retard/depot with original
    group_by(GEN) %>%  mutate(sex = check_sex(Sexe)) %>%  filter(!is.na(Sexe)) %>% ## if any sex is missing take sex from other entries
    mutate_at("PatientsTaille", as.numeric) %>% mutate(height = check_height(PatientsTaille)) %>% ### take average of all heights
    mutate_at(vars(Quetiapine:Doxepine), list(ever_drug = ever_drug)) %>% ungroup() %>%  ### create ever on any drug
    mutate(BMI = Poids/(PatientsTaille/100)^2)%>% filter(!is.na(BMI)) %>% ## create BMI
    group_by(GEN,AP1) %>% mutate(drug_instance = row_number()) %>%
    mutate(date_difference = as.numeric(difftime(lag(Date),Date, units = "days"))) %>%
    mutate(follow_up = case_when(abs(date_difference) >= (Mois-lag(Mois))*30-leeway_time &  abs(date_difference) <= (Mois-lag(Mois))*30+leeway_time ~ "sensible",
                                  is.na(date_difference) ~ "NA",
                                  Mois == 0 ~ "new_regimen",
                                  date_difference < 0 ~ "check",
                                  TRUE ~ "check")) %>%
    filter(!follow_up == "check") %>%
    mutate(AP1_mod = paste0(AP1, "_round",rename_meds(follow_up))) %>%
    ungroup() %>% group_by(GEN,AP1_mod) %>%
    mutate(AP1_duration= as.numeric(max(Date)-min(Date)))%>% mutate(AP1_duration = na_if(AP1_duration, 0)) %>%
    mutate(Num_followups = n()) %>%
    ungroup() %>% group_by(GEN,AP1) %>%
    filter(Num_followups == max(Num_followups)) %>%
    mutate(BMI_change = bmi_diff(BMI)) %>%
    mutate(time_between_visits = as.numeric(Date-lag(Date))) %>% replace_na(list(time_between_visits=0)) %>%
    dplyr::distinct(AP1, .keep_all=T) %>% ungroup() %>%
    rename_at(baseline_vars, function(x) paste0( x, "_start")) %>%
    mutate_at(paste0( baseline_vars, "_start"), destring) %>%
    group_by(GEN) %>%
    mutate(Drug_Number=paste0("Drug_",row_number())) %>%
    pivot_wider(id_cols=c(GEN,sex, ends_with("_ever_drug")), names_from=Drug_Number, values_from=c("AP1", "Age","Date", paste0(baseline_vars, "_start"),
        "AP1_duration", "BMI_change","Num_followups")) %>%
    dplyr::select(matches('GEN|sex|AP1|Age|Date|BMI|_start_Drug_1|Num_followups|_ever_drug')) %>%
    ungroup() %>%
    mutate(Age_sq_Drug_1 = Age_Drug_1^2) %>%
    mutate_at(vars(starts_with("AP1_Drug_")) , as.factor)
}

read_pcs <- function(pc_dir){

  pclist = list()
  for (eth in eths)
  {
    PC_temp<- fread(file.path(pc_dir, eth, paste0(study_name,"_",eth,"_projections.txt")))
    PC_temp <- PC_temp %>%
      separate(FID, c("counter", "GPCR"), sep="_")  %>%
      mutate(eth = eth)  %>%
      mutate(FID = IID)  %>%
      dplyr::select(FID, IID, counter, GPCR, eth, everything())
    pclist[[eth]] <- PC_temp # add it to your list

  }
  PC_data <- dplyr::bind_rows(pclist)
  return(PC_data)
}

munge_snp_weights <- function(snp_weights, related_inds){
  colnames(snp_weights) <- c("sample_ID", "population_label", "n_SNPs",
    "PC1", "PC2", "PC3", "YRI%", "CEU%", "EA%", "NA%")
  t <- snp_weights  %>%
    filter(!sample_ID %in% related_inds[,1]) %>%
    separate(sample_ID, into = c("ID", "GPCR"), sep = "_") %>%
    dplyr::select(-population_label)
  return(t)

}

clean_drugs <- function(x){
  drugs <- as.character(unlist(x %>% dplyr::select(starts_with("AP1_Drug_"))))
  time <- as.numeric(unlist(x %>% dplyr::select(starts_with("AP1_duration_Drug_"))))
  bmi_change <-  as.numeric(unlist(x %>% dplyr::select(starts_with("BMI_change_Drug_"))))
  num_followups <- as.numeric(unlist(x %>% dplyr::select(starts_with("Num_followups_Drug_"))))
  age_started <- as.numeric(unlist(x %>% dplyr::select(starts_with("Age_Drug_"))))
  drug_num <- c(1:length(drugs))
  out <- tibble(drug=drugs,time=time,bmi_change=bmi_change, drug_num=drug_num, age_started=age_started, num_followups=num_followups)
  return(out)
}

classify_drugs <- function(x, case_categories, preferential_control_categories) {
  drug_info <- clean_drugs(x)
  case_match <- unlist(drug_info %>% filter(time >= min_follow_up) %>%
    filter(drug %in% case_categories) %>% dplyr::select(drug_num))[1]
  control_match1 <- unlist(drug_info %>% filter(time >= min_follow_up) %>%
    filter(!drug %in% case_categories) %>% filter(drug %in% preferential_control_categories) %>%
    filter(num_followups==max(num_followups)) %>% dplyr::select(drug_num))[1]
  control_match2 <- unlist(drug_info %>% filter(time >= min_follow_up) %>%
    filter(!drug %in% case_categories) %>% filter(!drug %in% preferential_control_categories) %>%
    filter(num_followups==max(num_followups)) %>% dplyr::select(drug_num))[1]
  control_match <- ifelse(is.na(control_match1), control_match2,control_match1 )
  case <- unname(ifelse (is.na(case_match), 0, 1))
  if(is.na(control_match) & is.na(case_match)) {case <- NA}
  best_match <- unname(ifelse (is.na(case_match), control_match, case_match))
  time_temp <- dplyr::pull(drug_info[best_match,"time"])
  drug_temp <- dplyr::pull(drug_info[best_match,"drug"])
  bmi_change <- dplyr::pull(drug_info[best_match,"bmi_change"])
  age_started <- dplyr::pull(drug_info[best_match,"age_started"])
  return(list(case=case, drug_num=best_match, drug_name=drug_temp,age_started=age_started ,bmi_diff=bmi_change, duration=time_temp))
}


munge_pheno_follow <-  function(pheno_baseline, test_drugs) {
  out <- list()
  for(i in 1:dim(test_drugs)[1])
  {
    drug_list <- unlist(test_drugs %>% dplyr::select(drugs) %>% dplyr::slice(i))
    drug_class <- unlist(test_drugs %>% dplyr::select(class) %>% dplyr::slice(i))
    col_match <- paste(paste0(drug_list,"_ever_drug"), collapse = "|")
    followup_data <- pheno_baseline %>% rowwise() %>%
      dplyr::do({
            result = as_tibble(.)
            x =  classify_drugs(result,drug_list,low_inducers)
            result$high_inducer=x$case
            result$high_inducer_drug_num=x$drug_num
            result$high_inducer_drug_name=x$drug_name
            result$bmi_change=x$bmi_diff
            result$follow_up_time=x$duration
            result$age_started=x$age_started
            result
        }) %>% ungroup() %>%
      mutate(ever_drug_match = rowSums(dplyr::select(., matches(col_match)) == 1) > 0) %>%
      filter(!(ever_drug_match & high_inducer==0)) %>% ##filter out individuals who have taken this drug but it wasn't followed
      mutate(follow_up_time_sq = follow_up_time^2) %>%
      mutate(age_sq=age_started^2)

    out[[drug_class]] <- followup_data
  }
  return(out)
}

create_GWAS_pheno <- function(pheno_baseline, pheno_followup){
  na_to_none <- function(x) ifelse(is.na(x),'NONE',x)
  linear_pheno = pheno_baseline %>% dplyr::select(matches('^FID$|^IID$|_start_Drug_1'))
  linear_covar = pheno_baseline %>% dplyr::select(matches('^FID$|^IID$|^sex$|^Age_Drug_1$|Age_sq_Drug_1$|^PC[0-9]$|^PC[1][0-9]$|^PC20'))

  list_covar <- list(linear=linear_covar)
  list_pheno <- list(linear=linear_pheno)
  for(sub in names(pheno_followup))
  {
    data_drug <- pheno_followup[[sub]]
    all_var <- data_drug %>% dplyr::select(matches('^FID$|^IID$|bmi_change$|^high_inducer$|^age_|^follow_up_time',ignore.case = F)) %>%
      mutate_at("high_inducer", as.character) %>%
      mutate(high_inducer=dplyr::recode(high_inducer, "0" = "NoDrug", "1" = "Drug")) %>%
      rename_at(vars(-FID,-IID),function(x) paste0(x,"_",sub)) %>%
      mutate_if(.predicate = is.character,.funs = na_to_none)

    pheno_var <- all_var %>% dplyr::select(matches('^FID$|^IID$|^bmi_change|^high_inducer',ignore.case = F))
    covar_var <- all_var %>% dplyr::select(matches('^FID$|^IID$|^age_|^high_inducer|^follow_up_time',ignore.case = F))
    list_pheno[[sub]] <- pheno_var
    list_covar[[sub]] <- covar_var
  }

  full_pheno <- reduce(list_pheno, left_join, by = c("FID","IID")) %>%
    mutate_if(.predicate = is.character,.funs = na_to_none)
  full_covar <- reduce(list_covar, left_join, by = c("FID","IID")) %>%
    mutate_if(.predicate = is.character,.funs = na_to_none) %>%
    mutate_at(vars(starts_with("high_inducer_")), function(x) case_when(x == "NONE" ~ NA_real_, x == "Drug" ~ 1, x == "NoDrug" ~ 0))  %>%
    mutate_at("sex", function(x) case_when(x == "NONE" ~ NA_real_, x == "F" ~ 1, x == "M" ~ 0))
  out <- list(full_pheno = full_pheno, full_covar = full_covar)
}

split_followup <- function(pheno_followup, test_drugs){
  out <- list()
  for(i in 1:dim(test_drugs)[1])
  {
    drug_list <- unlist(test_drugs %>% dplyr::select(drugs) %>% dplyr::slice(i))
    drug_class <- unlist(test_drugs %>% dplyr::select(class) %>% dplyr::slice(i))
    data_drug <- pheno_followup[[drug_class]]
    temp_pheno = data_drug %>% dplyr::select(matches('^FID$|^IID$|bmi_change$'))
    temp_covar = data_drug %>% dplyr::select(matches('^FID$|^IID$|^sex$|^high_inducer$|^age_|^follow_up_time|^PC[0-9]$|^PC[1][0-9]$|^PC20',ignore.case = F)) %>%
      mutate(sex = ifelse(sex == "M",0,1))
    temp_list <- list("pheno"=temp_pheno,"covar"=temp_covar )
    out[[drug_class]] <- temp_list
  }
  return(out)
}

write_followup_data <- function(pheno_followup_split, test_drugs){
  for(i in 1:dim(test_drugs)[1])
  {
    drug_list <- unlist(test_drugs %>% dplyr::select(drugs) %>% dplyr::slice(i))
    drug_class <- unlist(test_drugs %>% dplyr::select(class) %>% dplyr::slice(i))
    data_drug <- pheno_followup_split[[drug_class]]
    out_interaction_pheno = write.table(data_drug$pheno, file_out(paste0("data/processed/phenotype_data/GWAS_input/",drug_class,"_interaction_pheno_input.txt")), row.names=F, quote=F, col.names=T)
    out_interaction_covar = write.table(data_drug$covar, file_out( paste0("data/processed/phenotype_data/GWAS_input/",drug_class,"_interaction_covar_input.txt")), row.names=F, quote=F, col.names=T)

  }
}

create_analysis_dirs <- function(top_level){
  dir.create(top_level,showWarnings=F)
  to_create <- c("analysis/GWAS/interaction","analysis/GWAS/full","analysis/GWAS/subgroup")
  for(dir in to_create){
    dir.create(dir,showWarnings=F)
    for (eth in eths){
      dir.create(file.path(dir, eth),showWarnings=F)
    }
  }
  dir.create("analysis/PRS",showWarnings=F)
  return(TRUE)
}

define_baseline_inputs <- function(GWAS_input, baseline_vars, drug_classes){
  col_match<- paste(paste0(baseline_vars,"_start_Drug_1"), collapse = "|")
  linear_vars <- colnames(dplyr::select(GWAS_input$full_pheno,matches(col_match)))
  linear_covars <- rep(list(c(standard_covars,baseline_covars)),length(baseline_vars))
  for(i in 1:length(drug_classes)){
    linear_vars <- c(linear_vars,
      colnames(GWAS_input$full_pheno)[which(colnames(GWAS_input$full_pheno)==paste0("bmi_change","_",drug_classes[i]))])
    covars <- c(standard_covars)
    covars <- c(covars, colnames((dplyr::select(GWAS_input$full_covar, ends_with(drug_classes[i])) %>%
        dplyr::select(starts_with("age")))))

    linear_covars[[length(linear_covars)+1]]  <- covars
  }

  baseline_gwas_info <- tibble ( pheno = linear_vars,
                                 covars = linear_covars,
                                 output_suffix = linear_vars)
  return(baseline_gwas_info)
}

define_interaction_inputs <- function(GWAS_input, drug_classes){
  interaction_vars <- numeric()
  interaction_covars <- list()
  interaction_pams <- numeric()
  for(i in 1:length(drug_classes)){

   interaction_vars <- c(interaction_vars,
     colnames(GWAS_input$full_pheno)[which(colnames(GWAS_input$full_pheno)==paste0("bmi_change","_",drug_classes[i]))])
   covars <- c(standard_covars)
   covars <- c(covars,
     colnames(dplyr::select(GWAS_input$full_covar, ends_with(drug_classes[i]))))
   interaction_covars[[length(interaction_covars)+1]]  <- covars
  }

  interaction_pams <- rep("1-27,49", length(drug_classes)) ## to check
  interaction_gwas_info <- tibble ( pheno = interaction_vars, ## this is y
                                    covars = interaction_covars, ## these are the x'x
                                    parameters = interaction_pams, ## this is passed to PLINK
                                    output_suffix = drug_classes) ## this will be the output name
  return(interaction_gwas_info)
}

define_subgroup_inputs <- function(GWAS_input, drug_classes){
  subgroup_vars <- numeric()
  subgroup_covars <- list()
  subgroup_pheno <- list()
  for(i in 1:length(drug_classes)){
    pheno <- c(paste0("bmi_change_",drug_classes[i]),paste0("high_inducer_",drug_classes[i]))
    subgroup_vars <- c(subgroup_vars,
      colnames(GWAS_input$full_pheno)[which(colnames(GWAS_input$full_pheno)==paste0("high_inducer","_",drug_classes[i]))])
    covars <- c(standard_covars)
    covars <- c(covars,
      colnames(dplyr::select(GWAS_input$full_covar, ends_with(drug_classes[i])) %>%
      dplyr::select(-starts_with("high_inducer"))))
    subgroup_covars[[length(subgroup_covars)+1]]  <- covars
    subgroup_pheno[[length(subgroup_pheno)+1]]  <- pheno
  }
  subgroup_gwas_info <- tibble( pheno = subgroup_pheno,
                                subgroup = subgroup_vars,
                                covars = subgroup_covars,
                                output_suffix = drug_classes)
  return(subgroup_gwas_info)
}


run_gwas <- function(pfile, pheno_file, pheno_name, covar_file, covar_names, eths,
  eth_sample_file, output_dir, output,
  remove_sample_file=NULL,threads=1,type="full",
  subgroup=NULL, subgroup_var="", interaction=FALSE,parameters=NULL, output_suffix=""){

  # eth_sample_file : in the form of "keep_this_ETH.txt"
  # default option: --variance-standardize
  file_name <- output
  file_name <- paste0(file_name,"_",output_suffix)
  if(type=="subgroup")
  {
    subgroup_commands <- c("--loop-cats", subgroup_var)
  } else subgroup_commands <- NULL

  if(!is.null(remove_sample_file)){
    remove_commands <- c("--remove",remove_sample_file)
  } else remove_commands <- NULL

  if(type=="interaction")
  {
    analysis_commands <- c("--glm", "interaction", "--parameters", parameters)
    file_name <- paste0(file_name,"_int")
  } else analysis_commands <- c("--glm", "hide-covar")

  general_commands <- c("--pfile", pfile,"--pheno", pheno_file, "--pheno-name", pheno_name, "--covar",covar_file,
          "--covar-name", covar_names,
          "--threads", threads, "--variance-standardize")

  out <- list()
  for (eth in eths)
  {
    keep_file <- str_replace(eth_sample_file, "ETH", eth)
    eth_count <- dim(fread(keep_file))[1]
    file_name <- paste0(file_name,"_",eth)
    write_dir <- file.path(output_dir, type)
    full_output <- file.path(write_dir,eth,file_name)
    eth_commands <- c("--keep", keep_file, "--out", full_output)

    #if(eth_count > 100)
    #{
      plink_input <- c(general_commands, analysis_commands, remove_commands, subgroup_commands,eth_commands)
      temp_out <- processx::run(command="plink2",plink_input, error_on_status=F)
      out[[eth]] <- temp_out
    #}
  }
  return(out)
}

meta_interaction <- function(run_gwas_interaction){
    input_char <- deparse(substitute(run_gwas_interaction))
    inter_variable <- sub("run_meta_interaction_", "", input_char)
    folder <- paste0("analysis/GWAS/interaction/",inter_variable)
    subfolders <- list.files(path=folder,full.names = T)
    out <- list()
    for(result_dir in subfolders)
    {
      output_var <- basename(result_dir)
      gwas_results <- list.files(path=result_dir, pattern=".filter", full.names=T)

      temp_out <- processx::run(command="plink",c( "--meta-analysis", gwas_results,  "+", "qt", "no-map",
        "--meta-analysis-snp-field", "ID", "--out", file.path(result_dir, paste0(inter_variable, ".", output_var)), "--threads", "16"), error_on_status=F)
      out[[output_var]] <- temp_out
    }

}

meta_linear <- function(run_gwas_linear){
  input_char <- deparse(substitute(run_gwas_linear))
  variable <- sub("run_gwas_linear_", "", input_char)

  folder <- paste0("analysis/GWAS/linear/",variable)
  result_dir <-folder
  gwas_results <- list.files(path=folder, pattern=".linear", full.names=T)

  out <-  processx::run(command="plink",c( "--meta-analysis", gwas_results,  "+", "qt", "no-map",
  "--meta-analysis-snp-field", "ID", "--out", file.path(result_dir, paste0(variable)), "--threads", "16"), error_on_status=F)

  return(out)
}

combine_targets <- function(...){
  temp <- substitute(list(...))[-1]
  c(sapply(temp, deparse))
}

process_subgroup <- function(nodrug, pheno_list, output = "PSYMETAB_GWAS", output_suffix = ""){

  # pheno_list <- subgroup_gwas_info$pheno
  # output_suffix <- subgroup_gwas_info$output_suffix
  in_folder <- file.path("analysis/GWAS/subgroup/CEU")
  gwas_files <- list.files(in_folder, ".linear")
  eth <- "CEU"

  info <- fread("analysis/QC/15_final_processing/PSYMETAB_GWAS.info")
  eth="CEU"
  freq <- fread(paste0("analysis/QC/15_final_processing/",eth, "/PSYMETAB_GWAS.CEU.afreq"), header=T) %>%
    rename(CHROM = "#CHROM")


  snp_metrics <- full_join(info, freq, by = c("CHROM", "ID", "REF", "ALT"))

  for(k in 1:length(pheno_list)){

    outcome_var <- pheno_list[[k]][1]
    suffix <- output_suffix[k]
    file_drug <- paste0(output,"_", suffix,"_", eth, ".Drug.", outcome_var, ".glm.linear" )
    file_nodrug <- paste0(output,"_", suffix,"_", eth, ".NoDrug.", outcome_var, ".glm.linear" )

    if(file_drug %in% gwas_files & file_nodrug %in% gwas_files){

    nodrug <- fread(file.path(in_folder,file_nodrug), data.table=F, stringsAsFactors=F) %>%
      rename(BP = POS) %>% filter(!is.na(P))
    drug <- fread(file.path(in_folder,file_drug), data.table=F, stringsAsFactors=F) %>%
      rename(BP = POS) %>% filter(!is.na(P))

    joint <- inner_join(nodrug, drug, by = c("#CHROM", "BP", "ID", "TEST"), suffix = c("_nodrug", "_drug")) %>%
      rename(CHROM = "#CHROM") %>% left_join(snp_metrics, by = c("CHROM", "ID"))

    filter <- joint %>% filter(R2 > info_threshold) %>%
      filter(ALT_FREQS > maf_threshold & ALT_FREQS < (1 - maf_threshold))

    het_out <- numeric()
    for(snp in 1:dim(filter)[1] ){
      rsid <- as.character(filter[snp,"ID"])
      beta_D <-  as.numeric(as.character(filter[snp,"BETA_drug"]))
      beta_ND <-   as.numeric(as.character(filter[snp,"BETA_nodrug"]))
      SE_D <-  as.numeric(as.character(filter[snp,"SE_drug"]))
      SE_ND <-  as.numeric(as.character(filter[snp,"SE_nodrug"]))
      se <- sqrt( (SE_D^2) + (SE_ND^2) )
      t <- (beta_D-beta_ND)/se
      p_het <- 2*pnorm(-abs(t))
      temp <- cbind(rsid, p_het)
      #het_out <- rbind(het_out,temp)
      write.table(temp, paste0("analysis/GWAS/subgroup/", eth, "/", output, "_", suffix,"_", eth, ".", outcome_var, ".het"), append=T, col.names=F, row.names=F)
      if((snp %% 50000)==0){print(snp)}

    }

    }
  }

}

pca_plot <- function(data_clean, col, colname, title){
  ggplot(data_clean) +
  # add scatter points
  geom_point(aes_string(x = "PC1", y = "PC2", col = col),
             alpha = 0.5) +

  # change labels
  labs(title = title,
       x = "PC1",
       y = "PC2",
       #size = "Percent European",
       col = colname) +
  # change the size scale
  scale_size(range = c(0.1, 10)) +
  # add a nicer theme
  theme_classic() +
  # place legend at top and grey axis lines
  theme(legend.position = "top")
}

munge_snpweights <- function(study_name, pc_data, snp_weights, eth_file){

  colnames(snp_weights) <- c("FID", "IID", "num_SNPs",  paste0("snpweights_PC",1:3), "YRI_perc", "CEU_perc", "EA_perc", "NA_perc")

  pheno_eth <- read.table(eth_file,header=F) # reported ethnicity
  colnames(pheno_eth) <- c(  "FID", "IID", "ETH")
  pheno_eth <- pheno_eth %>%
    unite("ID2", FID:IID, remove = FALSE, sep="_")
    #mutate_at("ID2", as.factor)

  data <- merge(snp_weights, pc_data, by=c("FID", "IID"))

  factor_eths <- function(x) {ifelse(x< 0.8, FALSE, TRUE)}
  factor_yri <- function(x) {ifelse(x< 0.7, FALSE, TRUE)}

  data_clean <- data %>%
      mutate_at(c("CEU_perc","EA_perc","NA_perc"), list(factor=factor_eths)) %>%
      mutate_at(c("YRI_perc"), list(YRI_perc_factor=factor_yri)) %>%

      rename_at(vars(ends_with("_factor")),
        list( ~ str_replace(., "_perc", ""))) %>%
      dplyr::select(ends_with("_factor"), FID)  %>%
      tidyr::gather(eth, value, -FID) %>%
      filter(value==TRUE) %>%
      dplyr::select(-value) %>%
      left_join(data %>% dplyr::select(-ends_with("_factor")), .) %>%
      arrange(IID)  %>%
      mutate_at(c("YRI_perc","CEU_perc","EA_perc","NA_perc"), list(factor=factor_eths)) %>%
      rename_at(vars(ends_with("_factor")),
        list( ~ str_replace(., "_perc", "")))  %>%
      replace_na(list(eth = "MIXED"))  %>%
      mutate(genetic_eth = recode(eth, CEU_factor = "CEU",
           EA_factor = "EA",
           NA_factor = "NA",
           YRI_factor = "YRI"
       )) %>%
       mutate_at("IID", as.character) %>%
       left_join(., pheno_eth %>% dplyr::select(ID2,ETH), by = c("IID" = "ID2") ) %>%
       dplyr::rename(reported_eth = ETH)

  return(data_clean)

}

define_prs_inputs <- function(consortia_dir, output_folder){

  prs_base_files <- list.files(paste0(consortia_dir, "/formatted"), full.names=T)
  prs_names <- paste0(output_folder, "/", gsub(".txt", "", list.files(paste0(consortia_dir, "/formatted"))))

  prs_info <- tibble( base_file = prs_base_files,
                      out_file = prs_names)
  return(prs_info)

}

run_prsice <- function(base_file, threads=1, memory="7900", out_file, PRSice_dir="/data/sgg2/jenny/bin/PRSice/",
  snp_col="SNP", chr_col="CHR", effect_allele_col="EFFECT_ALLELE",
  other_allele_col="OTHER_ALLELE", beta_or_col="BETA", p_col="PVAL",
  data_type="bgen", bgen_file="", sample_file="", plink_file="",
  pheno_file="", pheno_col="",
  bar_levels="0.00000005,0.0000005,0.000005,0.00005,0.0005,0.005,0.05,0.1", no_regress=TRUE){
  if(data_type=="plink"){
    data_spec_commands <- c("--target", plink_file)
  }
  if(data_type=="bgen"){
    data_spec_commands <- c("--target", paste0(bgen_file, ",", sample_file), "--type", "bgen")
  }
  general_commands <- c(paste0(PRSice_dir, "PRSice.R"), "--dir", ".", "--prsice", paste0(PRSice_dir, "PRSice_linux"),
    "--base", base_file, "--thread", threads, "--snp", snp_col, "--chr", chr_col, "--A1", effect_allele_col, "--A2",
    other_allele_col, "--stat", beta_or_col, "--pvalue", p_col, "--out", out_file, "--bar-levels", bar_levels)
  if(no_regress==TRUE){
    regress_commands <- c("--no-regress")
  }
  if(no_regress==FALSE){
    regress_commands <- c("--pheno-file", pheno_file, "--pheno-col", pheno_col,
    "--ignore-fid")
  }
  prsice_input <- c(general_commands, data_spec_commands, regress_commands)
  out <- processx::run(command="Rscript", prsice_input, error_on_status=F)
  return(out)
}

#### TO BE REVISED
process_gwas <- function(outcome_variable,interaction_variable,model){

  if(model=="subgroup"){
    out_folder <- file.path("analysis/GWAS/interaction",model)


  }
  if(model=="interaction"){
    out_folder <- file.path("analysis/GWAS/interaction",interaction_variable, outcome_variable)
    eth_gwas_files <- list.files(path=out_folder, pattern = "\\.filter$",full.names=T)
  }
  if(model=="linear"){
    out_folder <- file.path("analysis/GWAS/interaction", outcome_variable)
    eth_gwas_files <- list.files(path=out_folder, pattern = "\\.linear$",full.names=T)

  }


  eths <- sapply(eth_gwas_files,function(x){ gsub(".*PSYMETAB_GWAS_(.*?)\\_.*", "\\1", x)})

  meta_file <- file_in(list.files(path=out_folder, pattern=".meta", full.names=T))
  gwas_files <- tibble(eth=c(eths,"META"), file=c(eth_gwas_files,meta_file))

  info <- fread("analysis/QC/15_final_processing/PSYMETAB_GWAS.info")
  for(i in 1:dim(gwas_files)[1])
  {
    result <- dplyr::pull(gwas_files %>% dplyr::select(file))[i]
    eth <- dplyr::pull(gwas_files %>% dplyr::select(eth))[i]
    gwas_result <- fread(result, data.table=F, stringsAsFactors=F)

    if(eth!="META")
    {
      gwas_munge <- gwas_result %>% rename(CHR = "#CHROM") %>% rename(BP = POS) %>% filter(!is.na(P))
      freq <- fread(paste0("analysis/QC/15_final_processing/",eth, "/PSYMETAB_GWAS.CEU.afreq"), header=T)
      joint <- reduce(list(gwas_munge,freq,info), full_join, by = "ID")
      sig <- joint %>%
              mutate_at("P", as.numeric) %>%
              filter(P < gw_sig) %>%
              filter(ALT_FREQS > maf_threshold) %>%
              filter(R2 > info_threshold)

    }


    joint_maf <- joint %>% filter(ALT_FREQS > maf_threshold & ALT_FREQS < (1- maf_threshold))%>% mutate_at("P", as.numeric)
    sig <- joint_maf  %>% filter(P < gw_sig)
    png("man_interaction.png", width=2000, height=1000, pointsize=18)
    manhattan(joint_maf)
    dev.off()

    png("interaction_qq2.png", width=2000, height=1000, pointsize=18)
    qq(joint_maf$P)
    dev.off()
    qq(gwas_result2$P)
  }

}


#### GWAS
# sig_nodrug <-nodrug[ which(nodrug$P < 5e-08),]
# sig_drug <- drug[which(drug$P < 5e-08),]
#
# info <- fread("/data/sgg2/jenny/projects/PSYMETAB/analysis/QC/15_final_processing/PSYMETAB_GWAS.info")
# eth="CEU"
# freq <- fread(paste0("/data/sgg2/jenny/projects/PSYMETAB/analysis/QC/15_final_processing/",eth, "/PSYMETAB_GWAS.CEU.afreq"), header=T)

# for(data in c("nodrug", "drug"))
# {
# gwas_result <- get(data)
# gwas_munge <- gwas_result %>% rename(CHR = "#CHROM") %>% rename(BP = POS) %>% filter(!is.na(P))
# joint <- reduce(list(gwas_munge,freq,info), full_join, by = "ID")
# sig <- joint %>%
#         mutate_at("P", as.numeric) %>%
#         filter(P < 5e-06) %>%
#         filter(ALT_FREQS > maf_threshold & ALT_FREQS < (1 - maf_threshold)) %>%
#         filter(R2 > info_threshold)
# # sig_nodrug <- sig
# }

#
# for(pam in pams )
# {
#   beta_F <-  as.numeric(as.character(trait_result[pam,"female"]))
#   beta_M <-   as.numeric(as.character(trait_result[pam,"male"]))
#   SE_F <-  as.numeric(as.character(trait_result[paste0(pam,"_se"),"female"]))
#   SE_M <-  as.numeric(as.character(trait_result[paste0(pam,"_se"),"male"]))
#   se <- sqrt( (SE_F^2) + (SE_M^2) )
#   t <- (beta_F-beta_M)/se
#   p_het <- 2*pnorm(-abs(t))
#   het_out <- cbind(het_out, p_het)
# }
