
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

check_sex <- function(x){out <- NA
  if(("M" %in% x) & ("F" %in% x)) {out <- NA} else {
                            if("M" %in% x) {out <- "M"}
                            if("F" %in% x) {out <- "F"}}
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


munge_pheno <- function(pheno_raw){
  pheno_raw %>%
    mutate(Date = as.Date(Date, format = '%d.%m.%Y'))  %>% filter(!is.na(Date)) %>% arrange(Date)  %>%
    mutate(AP1 = gsub(" ", "_",AP1)) %>% mutate_at("AP1",as.factor) %>% mutate(AP1 = gsub("_.*$","", AP1)) %>% mutate(AP1 = na_if(AP1, "")) %>%## merge retard/depot with original
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
    pivot_wider(id_cols=c(GEN,sex, ends_with("_ever_drug")), names_from=Drug_Number, values_from=c("AP1", "Age","Date", "BMI_start", "LDL_start","Glucose_start","Creatinine_start",
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


munge_pheno_follow <-  function(pheno_baseline) {

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

split_followup <- function(pheno_followup){
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

write_followup_data <- function(pheno_followup_split){

  for(i in 1:dim(test_drugs)[1])
  {
    drug_list <- unlist(test_drugs %>% dplyr::select(drugs) %>% dplyr::slice(i))
    drug_class <- unlist(test_drugs %>% dplyr::select(class) %>% dplyr::slice(i))
    data_drug <- pheno_followup_split[[drug_class]]

    out_interaction_pheno = write.table(data_drug$pheno, file_out(paste0("data/processed/phenotype_data/GWAS_input/",drug_class,"_interaction_pheno_input.txt")), row.names=F, quote=F, col.names=T)
    out_interaction_covar = write.table(data_drug$covar, file_out( paste0("data/processed/phenotype_data/GWAS_input/",drug_class,"_interaction_covar_input.txt")), row.names=F, quote=F, col.names=T)

  }
}

create_analysis_dirs <- function(){
  dir.create("analysis/GWAS",showWarnings=F)
  dir.create("analysis/GWAS/interaction",showWarnings=F)
  dir.create("analysis/GWAS/linear",showWarnings=F)
  dir.create("analysis/GRS",showWarnings=F)
  return(TRUE)
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

combine_targets <- function(...)
{
  temp <- substitute(list(...))[-1]
  c(sapply(temp, deparse))

}



#### TO BE REVISED
process_gwas <- function(outcome_variable,interaction_variable,model){

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
      freq <- fread(paste0("analysis/QC/15_final_processing/,",eth, "/PSYMETAB_GWAS.CEU.afreq"), header=T)
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
