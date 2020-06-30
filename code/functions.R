
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

get_group_number = function(){
    i = 0
    function(){
        i <<- i+1
        i
    }
}


rename_meds <- function(AP1, PatientsRecNum, Date)
{

  out_index <- numeric()
  data <- tibble(AP1 = AP1, PatientsRecNum = as.factor(PatientsRecNum), Date = as.Date(Date))

  PatientsRecNum_previous <- data[["PatientsRecNum"]][1]
  count_previous <- 1
  drugs_seen <- c()
  for(i in 1:dim(data)[1]){

    drug_i <- data[["AP1"]][i]
    PatientsRecNum_i <- data[["PatientsRecNum"]][i]
    if(drug_i %in% drugs_seen){

      if(PatientsRecNum_previous==PatientsRecNum_i){count <- count_previous}
      if(PatientsRecNum_previous!=PatientsRecNum_i){count <- count_previous + 1}
    }

    if(!drug_i %in% drugs_seen){
      count <- 1
    }
    out_index <- c(out_index, count)
    PatientsRecNum_previous <- PatientsRecNum_i
    count_previous <- count
    drugs_seen <- c(drugs_seen, drug_i)
  }

  return(paste0(data[["AP1"]], "_round", out_index))
}



bmi_diff <- function(x){
  n_obs <- length(x)
  if(length(x) > 1) {out <- x[n_obs]-x[1]}
  if(length(x) ==1) {out <- NA}
  return(out)
}

bmi_diff_time <- function(x, date_difference_first, time_limit){
  data <- tibble(BMI = x, date_difference_first = date_difference_first)
  bmi_diff(data %>% filter(date_difference_first <= time_limit) %>% pull(BMI))
}

calc_slope <- function(x, y, time_limit){
  x <- x[which(x <= time_limit)]
  y <- y[which(x <= time_limit)]
  beta <- NA
  if(length(x) > 1 & (length(unique(x)) >1 | length(unique(y)) > 1)){
    if(length(unique(y))==1){
      beta <- 0
    } else beta <- summary(lm(y~x))$coefficients[2,1]

  }

  return(beta)

}

calc_duration_short <- function(Date, date_difference_first, time_limit){
  data <- tibble(Date = Date, date_difference_first = date_difference_first)

  Date_6mo_max <- as.Date(data %>% filter(date_difference_first <= time_limit) %>% filter(Date==max(Date)) %>% pull(Date), format='%d.%m.%y')
  temp <- data %>% mutate(AP1_duration= case_when(date_difference_first <= time_limit ~ as.numeric(Date_6mo_max-min(Date)),
                                                  date_difference_first > time_limit ~ NA_real_)) %>%
                   mutate(AP1_duration = na_if(AP1_duration, 0)) %>% pull(AP1_duration)
  return(temp)
}

calc_weighted_slope <- function(x, y, time_limit){

  x <- x[which(x <= time_limit)]
  y <- y[which(x <= time_limit)]
  std_beta <- NA
  se <- 1
  if(length(x) > 1 & (length(unique(x)) >1 | length(unique(y)) > 1)){
    if(length(unique(y))==1){
      beta <- 0
    } else {
      beta <- summary(lm(y~x))$coefficients[2,1]
      se <- summary(lm(y~x))$coefficients[2,2]
    }
      std_beta <- beta/se
  }

  return(std_beta)

}

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
     length(unique(na.omit(Drug)))==1) {out <- "sensible"}
  else if (length(unique(na.omit(PatientRec)))==1 &
     length(unique(na.omit(Drug)))==2 & "Paliperidone" %in% Drug & "Risperidone" %in% Drug) {
       out <- "sensible"
     }
 else if (length(unique(na.omit(PatientRec)))==1 &
    length(unique(na.omit(Drug)))==2 & "Risperdal" %in% Drug & "Risperidone" %in% Drug) {
      out <- "sensible"
    } else{
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
  if(any(x==1, na.rm=TRUE)){out <- 1} else{
    out <- 0
  }
  return(out)
}

check_sleep_disorder <- function(x, caf){
  max_caf <- which.max(caf)
  out <- NA
  if((0 %in% x) & (1 %in% x)) {
    out <- x[max_caf]
  } else {
    if(0 %in% x) {out <- 0}
    if(1 %in% x) {out <- 1}
  }
  return(out)
}

munge_caffeine <- function(caffeine_raw){
  caffeine_raw %>%
    mutate_all(~replace(., . == 999, NA)) %>%
    mutate(Date = as.Date(Date)) %>%
    rename(Sleep_disorder = Sleep_disorders) %>%
    group_by(GEN) %>%
    filter(Date == min(Date)) %>%
    mutate_at( vars(-GEN, -Sexe, -Age, -Diag, -Hospital_stay, -starts_with("Sleep_disorder")), mean) %>%
    mutate(Sexe = case_when(Sexe == 1 ~ "F",
                            Sexe == 0 ~ "M")) %>%
    rename_at(vars(-GEN),function(x) paste0(x,"_caffeine")) %>%
    mutate(Diag_caffeine = str_replace_all(Diag_caffeine, "& ", "")) %>%
    mutate(Diag_caffeine = str_replace_all(Diag_caffeine, " ", "_")) %>%
    unique()

}

merge_pheno_caffeine <- function(pheno, caffeine, anonymization_error){

  pheno_temp <- pheno %>% mutate(Date = as.Date(Date, format = '%d.%m.%y')) %>%
    mutate(Date_temp_start = Date - anonymization_error) %>%
    mutate(Date_temp_end = Date + anonymization_error)

  caffeine_temp <- caffeine %>%
    mutate(Date_temp_start = Date_caffeine - anonymization_error) %>%
    mutate(Date_temp_end = Date_caffeine + anonymization_error)

  out <- fuzzy_full_join(
      pheno_temp, caffeine_temp,
      by = c(
        "GEN" = "GEN",
        "Sexe" = "Sexe_caffeine",
        "Date_temp_start" = "Date_temp_end",
        "Date_temp_end" = "Date_temp_start"
        ),
      match_fun = list(`==`, `==`, `<=`, `>=`)
    ) %>%
    filter(!(Sexe != Sexe_caffeine) | any(is.na(Sexe), is.na(Sexe_caffeine))) %>% #remove sex's that don't match between the two databases
    dplyr::select(-GEN.y, -Sexe_caffeine, -starts_with("Date_temp")) %>%
    rename(GEN = GEN.x)

  return(out)
}

remove_outliers <- function(x){
  out <- (x > mean(x, na.rm=T)-sd(x, na.rm=T)*3) & (x < mean(x, na.rm=T)+sd(x, na.rm=T)*3)
}

munge_pheno <- function(pheno_raw, baseline_vars, leeway_time, caffeine_munge, follow_up_limit){
  pheno_raw %>%
    mutate_all(~replace(., . == 999, NA)) %>% filter(!is.na(PatientsTaille) & !is.na(Poids)) %>%
    filter(remove_outliers(PatientsTaille)) %>%
    filter(remove_outliers(Poids)) %>%
    mutate(Date = as.Date(Date, format = '%d.%m.%y'))  %>%
    filter(!is.na(Date)) %>% arrange(PatientsRecNum, Date)  %>%
    mutate(AP1 = gsub(" ", "_",AP1)) %>% mutate_at("AP1",as.factor) %>% mutate(AP1 = gsub("_.*$","", AP1)) %>% mutate(AP1 = na_if(AP1, "")) %>% ## merge retard/depot with original
    filter(!is.na(AP1)) %>%
    group_by(GEN) %>%  mutate(sex = check_sex(Sexe)) %>%  filter(!is.na(Sexe)) %>% ## if any sex is missing take sex from other entries
    mutate_at("PatientsTaille", as.numeric) %>% mutate(height = check_height(PatientsTaille)) %>% ### take average of all heights
    mutate_at(vars(Quetiapine:Doxepine), list(ever_drug = ever_drug)) %>% ungroup() %>%  ### create ever on any drug
    rename(weight = Poids) %>%
    mutate(BMI = weight/(height/100)^2) %>% filter(!is.na(BMI)) %>% ## create BMI
    group_by(GEN,PatientsRecNum) %>% mutate(drug_instance = row_number()) %>%
    mutate(date_difference = as.numeric(difftime(lag(Date), Date, units = "days"))) %>%
    mutate(AP1 = case_when(AP1 == "Risperdal" ~ "Risperidone",
                             AP1 == "Paliperidone" ~ "Risperidone",
                             TRUE ~ AP1
                           )) %>%
    mutate(follow_up = case_when(abs(date_difference) >= (Mois-lag(Mois))*30-leeway_time &  abs(date_difference) <= (Mois-lag(Mois))*30+leeway_time ~ "sensible",
                                  is.na(date_difference) ~ "NA",
                                  Mois == 0 ~ "new_regimen",
                                  date_difference < 0 ~ "leeway_exceeds",
                                  TRUE ~ "dupliate")) %>%
    mutate(month_descrepency = case_when(Mois < lag(Mois) | date_difference >= 0 ~ "month_discrepency",
                                TRUE ~ "sensible")) %>%
    mutate(drug_match = check_drug(PatientsRecNum, AP1)) %>%
    mutate(date_difference_first =  as.numeric(difftime(Date, first(Date)), units = "days")) %>%
    ungroup() %>%

    group_by(GEN) %>%
    mutate(AP1_mod = rename_meds(AP1, PatientsRecNum, Date)) %>%
    ungroup() %>%

    group_by(GEN,AP1_mod) %>%
    mutate(AP1_duration= as.numeric(max(Date)-min(Date)))%>% mutate(AP1_duration = na_if(AP1_duration, 0)) %>%
    mutate(Num_followups = n()) %>%
    ungroup() %>% group_by(GEN,AP1) %>%
    filter(grepl('round1$', AP1_mod)) %>% # restrict to only round1
    #filter(Num_followups == max(Num_followups)) %>% # restrict to drug with most number of follow-ups
    mutate(BMI_change = bmi_diff(BMI)) %>%
    mutate(time_between_visits = as.numeric(Date-lag(Date))) %>% replace_na(list(time_between_visits=0)) %>%

    mutate(BMI_change_6mo = bmi_diff_time(BMI, date_difference_first, follow_up_limit)) %>% #mutate(BMI_change2=replace(BMI_change2, date_difference_first > follow_up_limit, NA)) %>%
    mutate(BMI_slope=tryCatch(calc_slope(date_difference_first, BMI, Inf),error=function(e){9999}, warning=function(w){999})) %>%
    mutate(BMI_slope_6mo=tryCatch(calc_slope(date_difference_first, BMI, follow_up_limit),error=function(e){9999}, warning=function(w){999})) %>%
    mutate(BMI_slope_weight=tryCatch(calc_weighted_slope(date_difference_first, BMI, Inf),error=function(e){9999}, warning=function(w){999})) %>%
    mutate(BMI_slope_weight_6mo=tryCatch(calc_weighted_slope(date_difference_first, BMI, follow_up_limit),error=function(e){9999}, warning=function(w){999})) %>%
    mutate(AP1_duration_6mo = calc_duration_short(Date, date_difference_first, follow_up_limit)) %>%

    mutate_at(vars(starts_with("BMI_slope")), ~replace(., . == 9999, NA)) %>% # there should be no errors or changes
    mutate_at(vars(starts_with("BMI_slope")), ~replace(., . == 999, NA)) %>% # warning messages result from "essentially perfect fit: summary may be unreliable" remove these rows

    dplyr::distinct(AP1, .keep_all=T) %>% ungroup() %>%
    rename_at(baseline_vars, function(x) paste0( x, "_start")) %>%
    mutate_at(paste0( baseline_vars, "_start"), destring) %>%
    group_by(GEN) %>%
    mutate(Drug_Number=paste0("Drug_",row_number())) %>%
    pivot_wider(id_cols=c(GEN,sex, ends_with("_ever_drug"), height), names_from=Drug_Number, values_from=c(starts_with("BMI_slope"), starts_with("BMI_change"), "AP1", "Age","Date", paste0(baseline_vars, "_start"), "BMI_slope",
        "AP1_duration", "AP1_duration_6mo", "Num_followups")) %>%
    dplyr::select(matches('GEN|sex|AP1|Age|Date|height|BMI|_start_Drug_1|Num_followups|_ever_drug'))%>%
    ungroup() %>%
    mutate(Age_sq_Drug_1 = Age_Drug_1^2) %>%
    mutate_at(vars(starts_with("AP1_Drug_")) , as.factor) %>%
    left_join(caffeine_munge, by = c("GEN" = "GEN")) %>%
    filter(!(sex != Sexe_caffeine) | any(is.na(sex), is.na(Sexe_caffeine))) %>% #make sure sex from caffeine and pheno data match
    dplyr::select(-Sexe_caffeine) %>% #remove sex from caffeine data
    mutate(Age_caffeine_sq = Age_caffeine^2)
}


# dplyr::select(GEN, follow_up, date_difference, Date, Mois, problems)


read_pcs <- function(pc_dir, study_name, eths){

  pclist = list()
  for (eth in eths)
  {
    if(file.exists(file.path(pc_dir, eth, paste0(study_name,"_",eth,"_projections.txt")))){
      PC_temp<- fread(file.path(pc_dir, eth, paste0(study_name,"_",eth,"_projections.txt")))
      PC_temp <- PC_temp %>%
        separate(FID, c("counter", "GPCR"), sep="_")  %>%
        mutate(eth = eth)  %>%
        mutate(FID = IID)  %>%
        dplyr::select(FID, IID, counter, GPCR, eth, everything())
      pclist[[eth]] <- PC_temp # add it to your list
    }
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
  time_6mo <- as.numeric(unlist(x %>% dplyr::select(starts_with("AP1_duration_6mo_Drug_"))))

  bmi_start <-  as.numeric(unlist(x %>% dplyr::select(starts_with("BMI_start_Drug_"))))

  bmi_change <-  as.numeric(unlist(x %>% dplyr::select(starts_with("BMI_change_Drug_"))))
  bmi_change_6mo <- as.numeric(unlist(x %>% dplyr::select(starts_with("BMI_change_6mo_Drug_"))))

  bmi_slope <-  as.numeric(unlist(x %>% dplyr::select(starts_with("BMI_slope_Drug_"))))
  bmi_slope_weight <- as.numeric(unlist(x %>% dplyr::select(starts_with("BMI_slope_weight_Drug_"))))

  bmi_slope_6mo <-  as.numeric(unlist(x %>% dplyr::select(starts_with("BMI_slope_6mo_Drug_"))))
  bmi_slope_weight_6mo <- as.numeric(unlist(x %>% dplyr::select(starts_with("BMI_slope_weight_6mo_Drug_"))))

  num_followups <- as.numeric(unlist(x %>% dplyr::select(starts_with("Num_followups_Drug_"))))
  age_started <- as.numeric(unlist(x %>% dplyr::select(starts_with("Age_Drug_"))))
  drug_num <- c(1:length(drugs))
  out <- tibble(drug=drugs,time=time,time_6mo=time_6mo,bmi_start=bmi_start, bmi_change=bmi_change, bmi_change_6mo=bmi_change_6mo,
    bmi_slope=bmi_slope, bmi_slope_weight=bmi_slope_weight, bmi_slope_6mo=bmi_slope_6mo, bmi_slope_weight_6mo=bmi_slope_weight_6mo,
    drug_num=drug_num, age_started=age_started, num_followups=num_followups)
  return(out)
}

classify_drugs <- function(x, case_categories, preferential_control_categories, preferential_noncontrol_categories) {
  drug_info <- clean_drugs(x)
  case_match <- unlist(drug_info %>% filter(time >= min_follow_up) %>%
    filter(drug %in% case_categories) %>% dplyr::select(drug_num))[1]
  control_match1 <- unlist(drug_info %>% filter(time >= min_follow_up) %>%
    filter(!drug %in% case_categories) %>% filter(drug %in% preferential_control_categories) %>%
    filter(drug_num==min(drug_num)) %>% dplyr::select(drug_num))[1] ## preferred control

  control_match2 <- unlist(drug_info %>% filter(time >= min_follow_up) %>%
    filter(!drug %in% case_categories) %>% filter(!drug %in% preferential_noncontrol_categories) %>%
    filter(!drug %in% preferential_control_categories) %>%
    filter(drug_num==min(drug_num)) %>% dplyr::select(drug_num))[1]

  control_match3 <- unlist(drug_info %>% filter(time >= min_follow_up) %>%
    filter(!drug %in% case_categories) %>% filter(drug %in% preferential_noncontrol_categories) %>%
    filter(drug_num==min(drug_num)) %>% dplyr::select(drug_num))[1] ## preferred non-control
  control_match <- control_match1
  if(is.na(control_match1)){

    control_match <- ifelse(is.na(control_match2), control_match3,control_match2 )
  }

  case <- unname(ifelse (is.na(case_match), 0, 1))
  if(is.na(control_match) & is.na(case_match)) {case <- NA}
  best_match <- unname(ifelse (is.na(case_match), control_match, case_match))
  bmi_start <- dplyr::pull(drug_info[best_match,"bmi_start"])
  time_temp <- dplyr::pull(drug_info[best_match,"time"])
  time_6mo_temp <- dplyr::pull(drug_info[best_match,"time_6mo"])

  drug_temp <- dplyr::pull(drug_info[best_match,"drug"])
  bmi_change <- dplyr::pull(drug_info[best_match,"bmi_change"])
  bmi_change_6mo <- dplyr::pull(drug_info[best_match,"bmi_change_6mo"])

  bmi_slope <- dplyr::pull(drug_info[best_match,"bmi_slope"])
  bmi_slope_weight <- dplyr::pull(drug_info[best_match,"bmi_slope_weight"])
  bmi_slope_6mo <- dplyr::pull(drug_info[best_match,"bmi_slope_6mo"])
  bmi_slope_weight_6mo <- dplyr::pull(drug_info[best_match,"bmi_slope_weight_6mo"])

  age_started <- dplyr::pull(drug_info[best_match,"age_started"])
  return(list(case=case, drug_num=best_match, drug_name=drug_temp,age_started=age_started, bmi_start=bmi_start, bmi_diff=bmi_change, bmi_diff_6mo=bmi_change_6mo,
    bmi_slope=bmi_slope, bmi_slope_weight=bmi_slope_weight, bmi_slope_6mo=bmi_slope_6mo, bmi_slope_weight_6mo=bmi_slope_weight_6mo,duration=time_temp,duration_6mo=time_6mo_temp))
}


munge_pheno_follow <-  function(pheno_baseline, test_drugs, i) {
  out <- list()
  #for(i in 1:dim(test_drugs)[1])
  #{
    drug_list <- unlist(test_drugs %>% dplyr::select(drugs) %>% dplyr::slice(i))
    drug_class <- unlist(test_drugs %>% dplyr::select(class) %>% dplyr::slice(i))
    col_match <- paste(paste0(drug_list,"_ever_drug"), collapse = "|")
    followup_data <- pheno_baseline %>% rowwise() %>%
      dplyr::do({
            result = as_tibble(.) # result <- pheno_baseline %>% slice(49)
            x =  classify_drugs(result,drug_list,low_inducers, high_inducers)
            result$high_inducer=x$case
            result$high_inducer_drug_num=x$drug_num
            result$high_inducer_drug_name=x$drug_name
            result$bmi_change=x$bmi_diff
            result$bmi_change_6mo=x$bmi_diff_6mo
            result$follow_up_time=x$duration
            result$follow_up_time_6mo=x$duration_6mo
            result$age_started=x$age_started
            result$bmi_start=x$bmi_start
            result$bmi_slope=x$bmi_slope
            result$bmi_slope_6mo=x$bmi_slope_6mo
            result$bmi_slope_weight=x$bmi_slope_weight
            result$bmi_slope_weight_6mo=x$bmi_slope_weight_6mo

            result
        }) %>% ungroup() %>%
      mutate(ever_drug_match = rowSums(dplyr::select(., matches(col_match)) == 1) > 0) %>%
      filter(!(ever_drug_match & high_inducer==0)) %>% ##filter out individuals who have taken this drug but it wasn't followed
      mutate(follow_up_time_sq = follow_up_time^2) %>%
      mutate(follow_up_time_6mo_sq = follow_up_time_6mo^2) %>%
      mutate(age_sq=age_started^2)

    out[[drug_class]] <- followup_data
  #}
  return(out)
}

create_GWAS_pheno <- function(pheno_baseline, pheno_followup, caffeine_vars, test_drug, high_inducers, med_inducers, low_inducers){
  na_to_none <- function(x) ifelse(is.na(x),'NONE',x)
  linear_pheno <- pheno_baseline %>% dplyr::select(matches('^FID$|^IID$|_start_Drug_1'), paste0(caffeine_vars, "_caffeine"))
  linear_covar <- pheno_baseline %>% dplyr::select(matches('^FID$|^IID$|^sex$|^Age_Drug_1$|Age_sq_Drug_1$|^PC[0-9]$|^PC[1][0-9]$|^PC20|^Age_caffeine'))

  list_covar <- list(linear=linear_covar)
  list_pheno <- list(linear=linear_pheno)
  for(sub in names(pheno_followup))
  {
    data_drug <- pheno_followup[[sub]]
    drugs_analyzed <- test_drugs[which(test_drugs$class==sub),] %>% pull(drugs) %>% unlist()
    remove_drugs <- numeric()
    count <- 1
    for(drug_list_name in c("high_inducers", "med_inducers", "low_inducers")){
      drug_list <- get(drug_list_name)
      if(!any(drugs_analyzed %in% drug_list)){
        remove_drug_i <- drug_list
        remove_drugs <- c(remove_drugs, remove_drug_i)
      }
      if(any(drugs_analyzed %in% drug_list)){
        break
        #remove_drug_i <- drug_list[which(!drug_list %in% drugs_analyzed)]
        #remove_drugs <- c(remove_drugs, remove_drug_i)
      }
      count <- count + 1
    }

    data_drug_sensitivity <- data_drug %>% mutate(high_inducer=replace(high_inducer, high_inducer_drug_name %in% remove_drugs, NA)) %>%
        mutate_at(interaction_outcome,  ~replace(high_inducer, high_inducer_drug_name %in% remove_drugs, NA)) %>%
      dplyr::select(FID, IID, high_inducer, interaction_outcome) %>% rename_at(vars(-FID,-IID),function(x) paste0(x,"_sensitivity"))

    all_var <- data_drug %>% dplyr::select(matches('^FID$|^IID$|^bmi_start$|^high_inducer$|^age_|^follow_up_time',ignore.case = F), interaction_outcome) %>%
      left_join(data_drug_sensitivity) %>%
      mutate_at(c("high_inducer", "high_inducer_sensitivity"), as.character) %>%
      mutate_at(vars(high_inducer, high_inducer_sensitivity),  ~(recode(.,"0" = "NoDrug", "1" = "Drug"))) %>%
      rename_at(vars(-FID,-IID),function(x) paste0(x,"_",sub)) %>%
      mutate_if(.predicate = is.character,.funs = na_to_none)

    pheno_var <- all_var %>% dplyr::select(matches('^FID$|^IID$|^high_inducer',ignore.case = F), paste0(interaction_outcome,"_", sub),
      paste0(interaction_outcome,"_sensitivity_", sub))
    covar_var <- all_var %>% dplyr::select(matches('^FID$|^IID$|^age_|^high_inducer|^bmi_start|^follow_up_time',ignore.case = F))
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

write_followup_data <- function(pheno_followup_split, test_drugs, interaction_pheno, interaction_covar){
  for(i in 1:dim(test_drugs)[1])
  {
    drug_list <- unlist(test_drugs %>% dplyr::select(drugs) %>% dplyr::slice(i))
    drug_class <- unlist(test_drugs %>% dplyr::select(class) %>% dplyr::slice(i))
    data_drug <- pheno_followup_split[[drug_class]]
    out_interaction_pheno = write.table(data_drug$pheno, interaction_pheno, row.names=F, quote=F, col.names=T) #file_out(paste0("data/processed/phenotype_data/GWAS_input/",drug_class,"_interaction_pheno_input.txt"))
    out_interaction_covar = write.table(data_drug$covar, interaction_covar, row.names=F, quote=F, col.names=T) #file_out( paste0("data/processed/phenotype_data/GWAS_input/",drug_class,"_interaction_covar_input.txt"))

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
    dir.create(file.path(dir, "META"),showWarnings=F)
    dir.create(file.path(dir, "processed"),showWarnings=F)
  }
  dir.create("analysis/PRS",showWarnings=F)
  return(TRUE)
}

define_baseline_inputs <- function(GWAS_input, baseline_vars, drug_classes, caffeine_vars, interaction_outcome){
  col_match<- paste(c(paste0(baseline_vars,"_start_Drug_1"), paste0(caffeine_vars, "_caffeine")), collapse = "|")
  linear_vars <- colnames(dplyr::select(GWAS_input$full_pheno,matches(col_match)))
  linear_covars <- c(rep(list(c(standard_covars,baseline_covars)),length(baseline_vars)),
    rep(list(c(standard_covars,caffeine_covars)),length(caffeine_vars)))
  for(i in 1:length(drug_classes)){
    for(j in 1:length(interaction_outcome)){
      for(analysis in c("_", "_sensitivity_")){
        outcome <- interaction_outcome[j]
        linear_vars <- c(linear_vars,
          colnames(GWAS_input$full_pheno)[which(colnames(GWAS_input$full_pheno) %in% paste0(outcome,analysis,drug_classes[i]))])

        covars <- c(standard_covars)
        covars <- c(covars,
          colnames((dplyr::select(GWAS_input$full_covar, ends_with(drug_classes[i])) %>%
            dplyr::select(matches('age|bmi_start_|follow_up_time_'))))
          )

        if(grepl("_6mo$", outcome)){
          covars <- covars[-which(covars %in% paste0(c("follow_up_time_", "follow_up_time_sq_"), drug_classes[i]))]
        }
        if(!grepl("_6mo$", outcome)){
          covars <- covars[-which(covars %in% paste0(c("follow_up_time_6mo_", "follow_up_time_6mo_sq_"), drug_classes[i]))]
        }

        linear_covars[[length(linear_covars)+1]]  <- covars

      }
    }
  }


  baseline_gwas_info <- tibble ( pheno = linear_vars,
                                 covars = linear_covars,
                                 output_suffix = linear_vars)
  return(baseline_gwas_info)
}

define_interaction_inputs <- function(GWAS_input, drug_classes, interaction_outcome){
  interaction_vars <- numeric()
  interaction_covars <- list()
  interaction_pams <- numeric()
  output_suffix <- numeric()
  drug <- numeric()
  for(i in 1:length(drug_classes)){
    for(j in 1:length(interaction_outcome)){
      for(analysis in c("_", "_sensitivity_"))
      {
        outcome <- interaction_outcome[j]
        interactoin_pam_temp <- "1-27,50"
        interaction_vars <- c(interaction_vars,
          colnames(GWAS_input$full_pheno)[which(colnames(GWAS_input$full_pheno)==paste0(outcome,analysis,drug_classes[i]))])
        covars <- c(standard_covars)
        covars <- c(covars,
          colnames(dplyr::select(GWAS_input$full_covar, ends_with(drug_classes[i]))))

        high_inducer_sens <- which(grepl('^high_inducer_sensitivity_', covars))
        high_inducer_reg <- which(covars == paste0("high_inducer_", drug_classes[i]))

        if(analysis=="_"){
          covars <- covars[-high_inducer_sens]
        }
        if(analysis=="_sensitivity_"){
          covars[high_inducer_reg] <- covars[high_inducer_sens]
          covars <- covars[-high_inducer_sens]
        }
        if(grepl("_6mo$", outcome)){
          covars <- covars[-which(covars %in% paste0(c("follow_up_time_", "follow_up_time_sq_"), drug_classes[i]))]
        }
        if(!grepl("_6mo$", outcome)){
          covars <- covars[-which(covars %in% paste0(c("follow_up_time_6mo_", "follow_up_time_6mo_sq_"), drug_classes[i]))]
        }

        interaction_covars[[length(interaction_covars)+1]]  <- covars
        interaction_pams <- c(interaction_pams, interactoin_pam_temp)
        output_suffix <- c(output_suffix, paste0(interaction_outcome[j], "_", drug_classes[i] ))
        drug <- c(drug, drug_classes[i])

      }

    }


  }

  interaction_gwas_info <- tibble ( pheno = interaction_vars, ## this is y
                                    covars = interaction_covars, ## these are the x'x
                                    parameters = interaction_pams, ## this is passed to PLINK
                                    output_suffix = output_suffix,
                                    drug = drug) ## this will be the output name
  return(interaction_gwas_info)
}

define_subgroup_inputs <- function(GWAS_input, drug_classes, interaction_outcome){
  subgroup_vars <- numeric()
  subgroup_covars <- list()
  subgroup_pheno <- list()
  output_suffix <- numeric()
  for(i in 1:length(drug_classes)){
    for(j in 1:length(interaction_outcome)){
      for(analysis in c("_", "_sensitivity_"))
      {
        outcome <- interaction_outcome[j]
        pheno <- c(paste0(outcome, "_",drug_classes[i]),paste0("high_inducer_",drug_classes[i]))
        subgroup_vars <- c(subgroup_vars,
          colnames(GWAS_input$full_pheno)[which(colnames(GWAS_input$full_pheno)==paste0("high_inducer",analysis,drug_classes[i]))])
        covars <- c(standard_covars)
        covars <- c(covars,
          colnames(dplyr::select(GWAS_input$full_covar, ends_with(drug_classes[i])) %>%
            dplyr::select(-starts_with("high_inducer"))))
        if(grepl("_6mo$", outcome)){
          covars <- covars[-which(covars %in% paste0(c("follow_up_time_", "follow_up_time_sq_"), drug_classes[i]))]
        }
        if(!grepl("_6mo$", outcome)){
          covars <- covars[-which(covars %in% paste0(c("follow_up_time_6mo_", "follow_up_time_6mo_sq_"), drug_classes[i]))]
        }
        subgroup_covars[[length(subgroup_covars)+1]]  <- covars
        subgroup_pheno[[length(subgroup_pheno)+1]]  <- pheno
        output_suffix <- c(output_suffix, paste0(interaction_outcome[j], "_", drug_classes[i]))
      }

    }

  }
  subgroup_gwas_info <- tibble( pheno = subgroup_pheno,
                                subgroup = subgroup_vars,
                                covars = subgroup_covars,
                                output_suffix = output_suffix)
  return(subgroup_gwas_info)
}

run_gwas <- function(pfile, pheno_file, pheno_name, covar_file, covar_names, eths,
  eth_sample_file, output_dir, output,
  remove_sample_file=NULL,threads=1,type="full",
  subgroup=NULL, subgroup_var="", interaction=FALSE,parameters=NULL, output_suffix="",
  freq_file="analysis/QC/15_final_processing/CEU/PSYMETAB_GWAS.CEU.afreq",
  eth_low_maf_file){

  # eth_sample_file : in the form of "keep_this_ETH.txt"
  # default option: --variance-standardize
  file_name <- output
  file_name <- paste0(file_name,"_",output_suffix)

  sample_file <- read.table(paste0(pfile, ".psam"), header=F)
  nsamples <- dim(sample_file)[1]


  if(type=="subgroup")
  {
    subgroup_commands <- c("--loop-cats", subgroup_var)
  } else subgroup_commands <- NULL

  if(!is.null(remove_sample_file)){
    remove_commands <- c("--remove",remove_sample_file)
    remove_samples <- read.table(remove_sample_file, header=F)
  } else remove_commands <- NULL

  if(type=="interaction")
  {
    analysis_commands <- c("--glm", "interaction", "--parameters", parameters)
    file_name <- paste0(file_name,"_int")
  } else analysis_commands <- c("--glm", "hide-covar")

  general_commands <- unlist(c("--pfile", pfile, "--read-freq", freq_file, "--pheno", pheno_file, "--pheno-name", pheno_name, "--covar", covar_file,
          "--covar-name", unlist(covar_names),
          "--threads", threads, "--variance-standardize"))
  #maf_input <- unlist(c("--pfile", pfile, "--make-pfile", "--threads", threads, "--maf", maf_threshold))


  out <- list()
  for (eth in eths)
  {
    keep_file <- str_replace(eth_sample_file, "ETH", eth)
    low_maf_eth_file <- str_replace(eth_low_maf_file, "ETH", eth)
    eth_samples <- read.table(keep_file, header=F)

    file_name_eth <- paste0(file_name,"_",eth)
    write_dir <- file.path(output_dir, type)
    full_output <- file.path(write_dir,eth,file_name_eth)

    eth_commands <- c("--keep", keep_file, "--out", full_output, "--exclude", low_maf_eth_file)#, "--maf", maf_threshold)

    final_sample_list <- sample_file %>%
      filter(V1 %in% eth_samples$V1) %>%
      filter(!(V1 %in% remove_samples$V1))

    eth_count <- dim(final_sample_list)[1]
    if(eth_count > 100)
    {

      #maf_filter <- processx::run(command="plink2", c(maf_input, eth_commands), error_on_status=F)

      plink_input <- c(general_commands, analysis_commands, remove_commands, subgroup_commands,eth_commands)
      temp_out <- processx::run(command="plink2",plink_input, error_on_status=F)
      out[[eth]] <- temp_out
    }
  }
  return(out)
  #return(paste(plink_input, collapse="|"))
}

meta <- function(output, output_suffix="", eths, pheno, output_dir = "analysis/GWAS", type = "full",
  threads=1, interaction_var = NA){

  write_dir <- file.path(output_dir, type, "META")
  file_name <- output
  file_name <- paste0(file_name,"_",output_suffix)
  file_name_drug <- NA
  file_name_nodrug <- NA
  if(type=="interaction")
  {
    file_name <- paste0(file_name,"_int")
    awk_col7 <- paste0("'NR==1 || $7 == ", '"ADDxhigh_inducer_', interaction_var, '" { print $0 }', "'")

  }

  out <- list()
  gwas_results <- numeric()
  gwas_results_drug <- numeric()
  gwas_results_nodrug <- numeric()
  for(eth in eths){

    file_name_eth <- paste0(file_name,"_",eth)
    eth_dir <- file.path(output_dir, type)

    if(type=="subgroup"){
      file_name_eth_drug <- paste0(file_name_eth, ".Drug")
      file_name_eth_nodrug <- paste0(file_name_eth, ".NoDrug")
      eth_output_drug <- file.path(eth_dir,eth,paste0(file_name_eth_drug, ".", pheno[[1]][1], ".glm.linear"))
      eth_output_nodrug <- file.path(eth_dir,eth,paste0(file_name_eth_nodrug, ".", pheno[[1]][1], ".glm.linear"))

      if(file.exists(eth_output_drug)){
        gwas_results_drug <- c(gwas_results_drug, eth_output_drug)
        file_name_drug <- paste0(file_name, ".Drug")
      }
      if(file.exists(eth_output_nodrug)){
        gwas_results_nodrug <- c(gwas_results_nodrug, eth_output_nodrug)
        file_name_nodrug <- paste0(file_name, ".NoDrug")
      }
    }

    if(type=="interaction"){

      eth_output <- file.path(eth_dir,eth,paste0(file_name_eth, ".", pheno, ".glm.linear"))
      eth_subset_output <- file.path(eth_dir,eth,paste0(file_name_eth, ".", pheno, ".glm.linear.interaction"))

      if(file.exists(eth_output)){
        system(paste0("awk ", awk_col7, " ", eth_output, " > ", eth_subset_output))
        gwas_results <- c(gwas_results, eth_subset_output)
      }
    }

    if(type=="full"){
      eth_output <- file.path(eth_dir,eth,paste0(file_name_eth, ".", pheno, ".glm.linear"))

      if(file.exists(eth_output)){
        gwas_results <- c(gwas_results, eth_output)
      }
    }

  }

  full_output <- file.path(write_dir,file_name)
  full_output_drug <- file.path(write_dir,file_name_drug)
  full_output_nodrug <- file.path(write_dir,file_name_nodrug)

  if(type == "interaction" | type == "full"){
  out[["meta_out"]] <- processx::run(command="plink",c( "--meta-analysis", gwas_results,  "+", "qt", "report-all", "no-map",
    "--meta-analysis-snp-field", "ID", "--out", full_output, "--threads", threads), error_on_status=F)
  }

  if(type=="subgroup"){
    out[["meta_out_drug"]] <- processx::run(command="plink",c( "--meta-analysis", gwas_results_drug,  "+", "qt", "report-all", "no-map",
      "--meta-analysis-snp-field", "ID", "--out", full_output_drug, "--threads", threads), error_on_status=F)
    out[["meta_out_nodrug"]] <- processx::run(command="plink",c( "--meta-analysis", gwas_results_nodrug,  "+", "qt", "report-all", "no-map",
      "--meta-analysis-snp-field", "ID", "--out", full_output_nodrug, "--threads", threads), error_on_status=F)
  }

  return(out)

}

combine_targets <- function(...){
  temp <- substitute(list(...))[-1]
  c(sapply(temp, deparse))
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

define_ukbb_inputs <- function(Neale_summary_dir, ukbb_files, output_folder){

  prs_base_files <- paste0(Neale_summary_dir, "/", ukbb_files$file)
  prs_names <- paste0(output_folder, "/", ukbb_files$name, "_Neale_UKBB")

  prs_info <- tibble( base_file = prs_base_files,
                      out_file = prs_names)
  return(prs_info)

}

run_prsice <- function(base_file, threads=1, memory="7900", out_file, PRSice_dir="/data/sgg2/jenny/bin/PRSice/",
  snp_col="SNP", chr_col="CHR", effect_allele_col="EFFECT_ALLELE",
  other_allele_col="OTHER_ALLELE", beta_or_col="BETA", p_col="PVAL",
  data_type="bgen", bgen_file="", sample_file="", plink_file="",
  pheno_file="", pheno_col="",
  bar_levels="0.00000005,0.0000005,0.000005,0.00005,0.0005,0.005,0.05,0.1,0.2,0.3,0.4,0.5,1", no_regress=TRUE, fastscore=TRUE){
  if(data_type=="plink"){
    data_spec_commands <- c("--target", plink_file)
  }
  if(data_type=="bgen"){
    data_spec_commands <- c("--target", paste0(bgen_file, ",", sample_file), "--type", "bgen")
  }
  general_commands <- c(paste0(PRSice_dir, "PRSice.R"), "--dir", ".", "--prsice", paste0(PRSice_dir, "PRSice_linux"),
    "--base", base_file, "--thread", threads, "--snp", snp_col, "--chr", chr_col, "--A1", effect_allele_col, "--A2",
    other_allele_col, "--stat", beta_or_col, "--pvalue", p_col, "--out", out_file, "--bar-levels", bar_levels, "--fastscore", fastscore)
  if(no_regress==TRUE){
    regress_commands <- c("--no-regress")
  }
  if(no_regress==FALSE){
    regress_commands <- c("--pheno-file", pheno_file, "--pheno-col", pheno_col,
    "--ignore-fid")
  }
  prsice_input <- c(general_commands, data_spec_commands, regress_commands)
  out <- processx::run(command="Rscript", prsice_input, error_on_status=F)
  return(prsice_input)
}

format_prs <- function(all_score_file, out_file){
  prs_data <- fread(all_score_file, data.table=F, header=T)
  prs_format <- prs_data  %>%
    separate(FID, into = c("ID", "GPCR"), sep = "_") %>%
    dplyr::select(-IID)
  write.table(prs_format, out_file, row.names=F, quote=F)
}


full_model_summary <- function(model){

  nr <- length(which(is.na(coef(model))))
  nc <- 4
  rnames <- names(which(summary(model)$aliased))
  cnames <- colnames(summary(model)$coefficients)
  mat_na <- matrix(data = NA,nrow = nr,ncol = nc,
         dimnames = list(rnames,cnames))
  mat_coef <- rbind(summary(model)$coefficients,mat_na)
  return(mat_coef)

}


extract_model_stats <- function(model_summary, coef_list){
  summary <- numeric()
  for(coef in coef_list)
  {
    beta <- model_summary[coef,1]
    se <- model_summary[coef,2]
    pval <- model_summary[coef,4]
    temp <- cbind(beta, se, pval)
    colnames(temp) <- paste0(coef, c("_beta", "_se","_pval"))
    summary <- cbind (summary,  temp)
  }
  return(summary)

}

analyze_prs <- function(prs_file, pheno_file, pc_eth_data, linear_pheno_columns, glm_pheno_columns, covars){

  prs_data <- fread(prs_file, data.table=F)
  pheno_data <- fread(pheno_file, data.table=F)
  colnames(prs_data) <- c("ID", "GEN", "PRS_5e08", "PRS_5e07", "PRS_5e06", "PRS_5e05", "PRS_5e04", "PRS_5e03", "PRS_5e02", "PRS_0.1")

  output <- numeric()
  pc_data <- pc_eth_data[,c("GPCR", "eth", paste0("PC", 1:10))] %>% rename("GEN" = "GPCR")
  eths <- unique(pc_data$eth)
  for(prs in c("PRS_5e08", "PRS_5e07", "PRS_5e06", "PRS_5e05", "PRS_5e04", "PRS_5e03", "PRS_5e02", "PRS_0.1")){
    prs_sub_data <- prs_data[,c("GEN", prs)] %>% rename("prs" = prs)
    for(eth in eths){
      pc_sub_data <- pc_data %>% filter(eth==!!eth) %>% dplyr::select(-eth)
      for(var in linear_pheno_columns){
        model_type <- "linear"
        outcome_data <- pheno_data[,c("GEN", var)] %>% rename("outcome" = var)
        covar_data <- pheno_data[,c("GEN", covars)]

        joint <- reduce(list(outcome_data,prs_sub_data,covar_data,pc_sub_data), inner_join, by = "GEN") %>%
          dplyr::select(-GEN) %>%
          mutate(prs = scale(prs))
        model <- (lm(outcome ~ ., data=joint))
        model_summary <- full_model_summary(model)
        model_stats <- extract_model_stats(model_summary, c("prs"))
        row <- cbind(prs, eth, var, model_type, model_stats)
        output <- rbind(output, row)
      }

      for(var in glm_pheno_columns){
        model_type <- "logistic"
        outcome_data <- pheno_data[,c("GEN", var)] %>% rename("outcome" = var)
        covar_data <- pheno_data[,c("GEN", covars)]

        joint <- reduce(list(outcome_data,prs_sub_data,covar_data,pc_sub_data), inner_join, by = "GEN") %>%
          dplyr::select(-GEN) %>%
          mutate(prs = scale(prs))
        model <- (glm(outcome ~ ., data=joint, family="binomial"))

        model <- (glm(outcome ~ prs + Age_caffeine, data=joint, family="binomial"))

        model_summary <- full_model_summary(model)
        model_stats <- extract_model_stats(model_summary, c("prs"))
        row <- cbind(prs, eth, var, model_type, model_stats)
        output <- rbind(output, row)
      }

    }
  }
  return(output)



}
#### TO BE REVISED

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

    unlink(paste0("analysis/GWAS/subgroup/", eth, "/", output, "_", suffix,"_", eth, ".", outcome_var, ".het"))
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

define_baseline_files <- function(info = baseline_gwas_info, output = "PSYMETAB_GWAS", eths,
  output_dir = "analysis/GWAS", type = "full"){
  file_name <- output
  output_suffix <- info$output_suffix
  file_name <- paste0(file_name,"_",output_suffix)
  pheno <- info$pheno
  pheno_names <- numeric()
  eth_keep <- numeric()
  eth_list <- numeric()
  write_keep <- numeric()
  for(eth in c(eths, "META")){
    eth_dir <- file.path(output_dir, type)
    if(eth!="META"){
      file_name_eth <- paste0(file_name,"_",eth)
      eth_output <- file.path(eth_dir,eth,paste0(file_name_eth, ".", pheno, ".glm.linear"))
      write_file <- gsub(".glm.linear", ".GWAS.txt", eth_output)

    }
    if(eth=="META"){
      file_name_eth <- file_name
      eth_output <- file.path(eth_dir,eth,paste0(file_name_eth, ".meta"))
      write_file <- gsub(".meta", ".GWAS.txt", eth_output)
    }

    eth_list <- c(eth_list, rep(eth, length(which(file.exists(eth_output)))))
    eth_keep <- c(eth_keep, eth_output[which(file.exists(eth_output))])
    pheno_names <- c(pheno_names, pheno[which(file.exists(eth_output))])
    write_keep <- c(write_keep, write_file[which(file.exists(eth_output))])
    write_keep <- gsub(paste0("/", eth, "/"), "/processed/", write_keep)
  }

  out <- tibble(eth = eth_list,
                pheno = pheno_names,
                drug = NA,
                file = eth_keep,
                write_file = write_keep)
  return(out)
}

define_interaction_files <- function(info = interaction_gwas_info, output = "PSYMETAB_GWAS", eths,
  output_dir = "analysis/GWAS", type = "interaction"){

  file_name <- output
  output_suffix <- info$output_suffix
  file_name <- paste0(file_name,"_",output_suffix)
  pheno <- info$pheno
  pheno_names <- numeric()
  eth_keep <- numeric()
  eth_list <- numeric()
  write_keep <- numeric()
  for(eth in c(eths, "META")){
    eth_dir <- file.path(output_dir, type)

    if(eth!="META"){
      file_name_eth <- paste0(file_name,"_int_",eth)
      #eth_output <- file.path(eth_dir,eth,paste0(file_name_eth, ".", pheno, ".glm.linear"))
      eth_output <- file.path(eth_dir,eth,paste0(file_name_eth, ".", pheno, ".glm.linear.interaction"))
      write_file <- gsub(".glm.linear.interaction", ".GWAS.txt", eth_output)
    }
    if(eth=="META"){
      file_name_eth <- paste0(file_name,"_int")
      eth_output <- file.path(eth_dir,eth,paste0(file_name_eth, ".meta"))
      write_file <- gsub(".meta", ".GWAS.txt", eth_output)
    }
    eth_list <- c(eth_list, rep(eth, length(which(file.exists(eth_output)))))
    eth_keep <- c(eth_keep, eth_output[which(file.exists(eth_output))])
    pheno_names <- c(pheno_names, pheno[which(file.exists(eth_output))])
    write_keep <- c(write_keep, write_file[which(file.exists(eth_output))])
    write_keep <- gsub(paste0("/", eth, "/"), "/processed/", write_keep)
  }

  out <- tibble(eth = eth_list,
                pheno = pheno_names,
                drug = NA,
                file = eth_keep,
                write_file = write_keep)

  return(out)

}

define_subgroup_files <- function(info = subgroup_gwas_info, output = "PSYMETAB_GWAS", eths,
  output_dir = "analysis/GWAS", type = "subgroup"){

  file_name <- output
  output_suffix <- info$output_suffix
  file_name <- paste0(file_name,"_",output_suffix)
  pheno <- unlist(map(info$pheno, 1))
  pheno_names <- numeric()
  eth_keep <- numeric()
  eth_list <- numeric()
  drug_keep <- numeric()
  write_keep <- numeric()

  for(eth in c(eths, "META")){
    eth_dir <- file.path(output_dir, type)

    if(eth!="META"){
      file_name_eth <- paste0(file_name,"_",eth)
      file_name_eth_drug <- paste0(file_name_eth, ".Drug")
      file_name_eth_nodrug <- paste0(file_name_eth, ".NoDrug")
      eth_output_drug <- file.path(eth_dir,eth,paste0(file_name_eth_drug, ".", pheno, ".glm.linear"))
      eth_output_nodrug <- file.path(eth_dir,eth,paste0(file_name_eth_nodrug, ".", pheno, ".glm.linear"))
      write_file_drug <- gsub(".glm.linear", ".GWAS.txt", eth_output_drug)
      write_file_nodrug <- gsub(".glm.linear", ".GWAS.txt", eth_output_nodrug)
    }
    if(eth=="META"){
      file_name_eth <- file_name
      file_name_eth_drug <- paste0(file_name_eth, ".Drug")
      file_name_eth_nodrug <- paste0(file_name_eth, ".NoDrug")
      eth_output_drug <- file.path(eth_dir,eth,paste0(file_name_eth_drug, ".meta"))
      eth_output_nodrug <- file.path(eth_dir,eth,paste0(file_name_eth_nodrug, ".meta"))
      write_file_drug <- gsub(".meta", ".GWAS.txt", eth_output_drug)
      write_file_nodrug <- gsub(".meta", ".GWAS.txt", eth_output_nodrug)

    }
    eth_output <- c(eth_output_drug, eth_output_nodrug)
    drug_list <- c(rep("Drug", length(pheno)), rep("NoDrug", length(pheno)))
    drug_keep <- c(drug_keep, drug_list[which(file.exists(eth_output))])
    eth_list <- c(eth_list, rep(eth, length(which(file.exists(eth_output)))))
    eth_keep <- c(eth_keep, eth_output[which(file.exists(eth_output))])
    pheno_names <- c(pheno_names, c(pheno,pheno)[which(file.exists(eth_output))])
    write_keep <- c(write_keep, c(write_file_drug, write_file_nodrug)[which(file.exists(eth_output))])
    write_keep <- gsub(paste0("/", eth, "/"), "/processed/", write_keep)
  }

  out <- tibble(eth = eth_list,
                pheno = pheno_names,
                drug = drug_keep,
                file = eth_keep,
                write_file = write_keep)

  return(out)

}

process_gwas <- function(eth, pheno, drug, file, output = "PSYMETAB_GWAS", output_dir = "analysis/GWAS", type = "full",
  info_file = "analysis/QC/15_final_processing/PSYMETAB_GWAS.info", out_file){

  #outcome_variable,interaction_variable,model,
  info <- fread(info_file)
  write_dir <- file.path(output_dir, type, "processed")

  gwas_result <- fread(file, data.table=F, stringsAsFactors=F)
  if(eth!="META"){
    joint <- full_join(gwas_result, info, by = "ID") %>%
      rename("CHR" = "#CHROM") %>% rename("BP" = "POS.x") %>%
      filter(!is.na(P)) %>%
      dplyr::select(CHR, BP, ID, REF.x, ALT.x, BETA, SE, T_STAT, P) %>%
      rename("REF" = "REF.x") %>% rename("ALT" = "ALT.x") %>%
      rename("SNP" = "ID")

    # gw_sig_result <- joint %>%
    #    mutate_at("P", as.numeric) %>%
    #    filter(P < gw_sig)
  }

  if(eth=="META"){
    joint <- full_join(gwas_result, info, by = c("SNP" = "ID")) %>%
      rename("CHR" = "CHROM") %>% rename("BP" = "POS") %>%
      filter(!is.na(P)) %>%
      dplyr::select(CHR, BP, SNP, REF, ALT, BETA, P)
  }
  fwrite(joint, out_file, sep="\t")
  #return(list(manhattan_file_name=manhattan_file_name, qq_file_name=qq_file_name, title=title))
}

gwas_figures_input <- function(eth, pheno, drug, file, output = "PSYMETAB_GWAS", output_dir = "analysis/GWAS", type = "full",
  info_file = "analysis/QC/15_final_processing/PSYMETAB_GWAS.info", out_file){

  write_dir <- file.path(output_dir, type, "processed")

  title <- ifelse(is.na(drug), paste0(pheno, "_", eth), paste0(pheno, "_", drug, "_", eth))
  manhattan_file_name <- paste0(write_dir, "/", output, "_", title, "_manhattan.png")
  qq_file_name <- paste0(write_dir, "/", output, "_", title, "_qq.png")
  return(list(joint_file=out_file, manhattan_file_name=manhattan_file_name, qq_file_name=qq_file_name, title=title))

}

create_figures <- function(joint_file, manhattan_file_name, qq_file_name, title){

  joint <- fread(joint_file, data.table=F)

  png(manhattan_file_name, width=2000, height=1000, pointsize=18)
  manhattan(joint, main=title)
  dev.off()

  png(qq_file_name, width=2000, height=1000, pointsize=18)
  qq(joint$P, main=title)
  dev.off()

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
