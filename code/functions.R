
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
      previous_drug <- length(drugs_seen)
      previous_drug_i <- max(which(drugs_seen==drug_i))
      if(previous_drug == previous_drug_i){
        if(PatientsRecNum_previous==PatientsRecNum_i){count <- count_previous}
        if(PatientsRecNum_previous!=PatientsRecNum_i){count <- count_previous + 1}
      }
      if(previous_drug > previous_drug_i){
        count <- out_index[previous_drug_i] + 1
      }

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



X_diff <- function(x){
  n_obs <- length(x)
  if(length(x) > 1) {out <- x[n_obs]-x[1]}
  if(length(x) ==1) {out <- NA}
  return(out)
}

X_diff_time <- function(x, date_difference_first, time_limit){
  data <- tibble(X = x, date_difference_first = date_difference_first)
  X_diff(data %>% filter(date_difference_first <= time_limit) %>% pull(X))
}

calc_slope <- function(x, y, time_limit){
  NA_locations <- union(which(is.na(x)), which(is.na(y)))
  if(length(NA_locations) > 1){
    x <- x[-NA_locations]
    y <- y[-NA_locations]
  }

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
  NA_locations <- union(which(is.na(x)), which(is.na(y)))
  if(length(NA_locations) > 1){
    x <- x[-NA_locations]
    y <- y[-NA_locations]
  }

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
  out <- (x > mean(x, na.rm=T)-sd(x, na.rm=T)*5) & (x < mean(x, na.rm=T)+sd(x, na.rm=T)*5)
  return(out)
}

ivt <- function(x){
  out <- qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
  return(out)
}

munge_pheno <- function(pheno_raw, baseline_vars, leeway_time, caffeine_munge, follow_up_6mo, follow_up_3mo, follow_up_1mo, interaction_outcome){
  pheno_raw %>%
    mutate_all(~replace(., . == 999, NA)) %>% filter(!is.na(PatientsTaille) & !is.na(Poids)) %>%
    mutate(Date = as.Date(Date, format = '%d.%m.%y'))  %>%
    filter(!is.na(Date)) %>% arrange(PatientsRecNum, Date)  %>%
    mutate(AP1 = gsub(" ", "_",AP1)) %>% mutate_at("AP1",as.factor) %>% mutate(AP1 = gsub("_.*$","", AP1)) %>% mutate(AP1 = na_if(AP1, "")) %>% ## merge retard/depot with original
    filter(!is.na(AP1)) %>%
    group_by(GEN) %>%  mutate(sex = check_sex(Sexe)) %>%  filter(!is.na(Sexe)) %>% ## if any sex is missing take sex from other entries
    ungroup() %>%
    filter(remove_outliers(PatientsTaille)) %>%  # this removes patients with the following heights (cm): 106 106 106 106 106 106 106 106 116  96  96  90
    #filter(remove_outliers(Poids)) %>%
    group_by(GEN) %>%
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
    filter(!AP1_duration < min_follow_up) %>%
    filter(!is.na(AP1_duration)) %>%
    ungroup() %>% group_by(GEN,AP1) %>%
    mutate(AP1_mod = rename_meds(AP1, PatientsRecNum, Date)) %>% # rename drug rounds after filtering for rounds shorter than minimum (14 days)
    ungroup() %>% group_by(GEN,AP1) %>%
    filter(grepl('round1$', AP1_mod)) %>% # restrict to only round1 of each drug
    #filter(Num_followups == max(Num_followups)) %>% # restrict to drug with most number of follow-ups
      mutate_at(baseline_vars, list(change = ~X_diff(.))) %>%
    mutate(time_between_visits = as.numeric(Date-lag(Date))) %>% replace_na(list(time_between_visits=0)) %>%
    mutate_at(baseline_vars, as.numeric) %>%
    mutate_at(baseline_vars, list(change_6mo = ~X_diff_time(., date_difference_first, follow_up_6mo),
                                  change_3mo = ~X_diff_time(., date_difference_first, follow_up_3mo),
                                  change_1mo = ~X_diff_time(., date_difference_first, follow_up_1mo))) %>%

    mutate_at(baseline_vars, list(slope = ~tryCatch(calc_slope(date_difference_first, ., Inf),error=function(e){9999}, warning=function(w){999}),
                                  slope_6mo= ~tryCatch(calc_slope(date_difference_first, ., follow_up_6mo),error=function(e){9999}, warning=function(w){999}))) %>%

    mutate_at(baseline_vars, list(slope_weight = ~tryCatch(calc_weighted_slope(date_difference_first, ., Inf),error=function(e){9999}, warning=function(w){999}),
                                  slope_weight_6mo = ~tryCatch(calc_weighted_slope(date_difference_first, ., follow_up_6mo),error=function(e){9999}, warning=function(w){999}))) %>%

    mutate(AP1_duration_6mo = calc_duration_short(Date, date_difference_first, follow_up_6mo)) %>%
    mutate(AP1_duration_3mo = calc_duration_short(Date, date_difference_first, follow_up_3mo)) %>%
    mutate(AP1_duration_1mo = calc_duration_short(Date, date_difference_first, follow_up_1mo)) %>%

    mutate_at(vars(matches('_slope$|_slope_6mo$|_slope_weight$|_slope_weight_6mo$')), ~replace(., . == 9999, NA)) %>% # there should be no errors or changes
    mutate_at(vars(matches('_slope$|_slope_6mo$|_slope_weight$|_slope_weight_6mo$')), ~replace(., . == 999, NA)) %>% # warning messages result from "essentially perfect fit: summary may be unreliable" remove these rows

    dplyr::distinct(AP1, .keep_all=T) %>% ungroup() %>%
    rename_at(baseline_vars, function(x) paste0( x, "_start")) %>%
    mutate_at(paste0( baseline_vars, "_start"), destring) %>%
    group_by(GEN) %>%
    mutate(Drug_Number=paste0("Drug_",row_number())) %>%
    pivot_wider(id_cols=c(GEN,sex, ends_with("_ever_drug"), height), names_from=Drug_Number,
                values_from=c(interaction_outcome, "AP1", "Age","Date", paste0(baseline_vars, "_start"),
                              "AP1_duration", "AP1_duration_1mo", "AP1_duration_3mo", "AP1_duration_6mo", "Num_followups")) %>%

    dplyr::select(matches(paste0(paste(interaction_outcome, collapse="|"), '|GEN|sex|AP1|Age|Date|height|_start_Drug_1|Num_followups|_ever_drug'))) %>%
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

define_pheno_drug_combos <- function(drug_classes, baseline_vars){

  df <- (expand.grid(drug_classes, baseline_vars))
  colnames(df) <- c("class", "pheno")
  out <- as_tibble(df)
}

clean_drugs <- function(x, phenotype){
  drugs <- as.character(unlist(x %>% dplyr::select(starts_with("AP1_Drug_"))))
  time <- as.numeric(unlist(x %>% dplyr::select(starts_with("AP1_duration_Drug_"))))
  time_1mo <- as.numeric(unlist(x %>% dplyr::select(starts_with("AP1_duration_1mo_Drug_"))))
  time_3mo <- as.numeric(unlist(x %>% dplyr::select(starts_with("AP1_duration_3mo_Drug_"))))
  time_6mo <- as.numeric(unlist(x %>% dplyr::select(starts_with("AP1_duration_6mo_Drug_"))))

  pheno_start <-  as.numeric(unlist(x %>% dplyr::select(starts_with(paste0(phenotype, "_start_Drug_")))))

  pheno_change <-  as.numeric(unlist(x %>% dplyr::select(starts_with(paste0(phenotype, "_change_Drug_")))))
  pheno_change_1mo <- as.numeric(unlist(x %>% dplyr::select(starts_with(paste0(phenotype, "_change_1mo_Drug_")))))
  pheno_change_3mo <- as.numeric(unlist(x %>% dplyr::select(starts_with(paste0(phenotype, "_change_3mo_Drug_")))))
  pheno_change_6mo <- as.numeric(unlist(x %>% dplyr::select(starts_with(paste0(phenotype, "_change_6mo_Drug_")))))

  pheno_slope <-  as.numeric(unlist(x %>% dplyr::select(starts_with(paste0(phenotype, "_slope_Drug_")))))
  pheno_slope_weight <- as.numeric(unlist(x %>% dplyr::select(starts_with(paste0(phenotype, "_slope_weight_Drug_")))))

  pheno_slope_6mo <-  as.numeric(unlist(x %>% dplyr::select(starts_with(paste0(phenotype, "_slope_6mo_Drug_")))))
  pheno_slope_weight_6mo <- as.numeric(unlist(x %>% dplyr::select(starts_with(paste0(phenotype, "_slope_weight_6mo_Drug_")))))

  num_followups <- as.numeric(unlist(x %>% dplyr::select(starts_with("Num_followups_Drug_"))))
  age_started <- as.numeric(unlist(x %>% dplyr::select(starts_with("Age_Drug_"))))
  drug_num <- c(1:length(drugs))
  out <- tibble(drug=drugs, time=time, time_1mo=time_1mo, time_3mo=time_3mo, time_6mo=time_6mo,
    pheno_start=pheno_start, pheno_change=pheno_change,
    pheno_change_1mo=pheno_change_1mo, pheno_change_3mo=pheno_change_3mo, pheno_change_6mo=pheno_change_6mo,
    pheno_slope=pheno_slope, pheno_slope_weight=pheno_slope_weight, pheno_slope_6mo=pheno_slope_6mo, pheno_slope_weight_6mo=pheno_slope_weight_6mo,
    drug_num=drug_num, age_started=age_started, num_followups=num_followups)
  return(out)
}

classify_drugs <- function(x, phenotype, case_categories, preferential_control_categories, preferential_noncontrol_categories) {
  drug_info <- clean_drugs(x, phenotype)
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
  pheno_start <- dplyr::pull(drug_info[best_match,"pheno_start"])
  time_temp <- dplyr::pull(drug_info[best_match,"time"])
  time_1mo_temp <- dplyr::pull(drug_info[best_match,"time_1mo"])
  time_3mo_temp <- dplyr::pull(drug_info[best_match,"time_3mo"])
  time_6mo_temp <- dplyr::pull(drug_info[best_match,"time_6mo"])

  drug_temp <- dplyr::pull(drug_info[best_match,"drug"])
  pheno_change <- dplyr::pull(drug_info[best_match,"pheno_change"])
  pheno_change_1mo <- dplyr::pull(drug_info[best_match,"pheno_change_1mo"])
  pheno_change_3mo <- dplyr::pull(drug_info[best_match,"pheno_change_3mo"])
  pheno_change_6mo <- dplyr::pull(drug_info[best_match,"pheno_change_6mo"])

  pheno_slope <- dplyr::pull(drug_info[best_match,"pheno_slope"])
  pheno_slope_weight <- dplyr::pull(drug_info[best_match,"pheno_slope_weight"])
  pheno_slope_6mo <- dplyr::pull(drug_info[best_match,"pheno_slope_6mo"])
  pheno_slope_weight_6mo <- dplyr::pull(drug_info[best_match,"pheno_slope_weight_6mo"])

  age_started <- dplyr::pull(drug_info[best_match,"age_started"])
  return(list(case=case, drug_num=best_match, drug_name=drug_temp,age_started=age_started, pheno_start=pheno_start, pheno_diff=pheno_change,
    pheno_change_1mo=pheno_change_1mo, pheno_change_3mo=pheno_change_3mo, pheno_change_6mo=pheno_change_6mo,
    pheno_slope=pheno_slope, pheno_slope_weight=pheno_slope_weight, pheno_slope_6mo=pheno_slope_6mo, pheno_slope_weight_6mo=pheno_slope_weight_6mo,
    duration=time_temp, duration_1mo=time_1mo_temp, duration_3mo=time_3mo_temp, duration_6mo=time_6mo_temp))
}


munge_pheno_follow <-  function(pheno_baseline, test_drugs, class, phenotype, low_inducers, high_inducers) {
  out <- list()
  #for(i in 1:dim(test_drugs)[1])
  #{
    i <- which(test_drugs$class==class)
    drug_list <- unlist(test_drugs %>% dplyr::select(drugs) %>% dplyr::slice(i))
    drug_class <- unlist(test_drugs %>% dplyr::select(class) %>% dplyr::slice(i))
    col_match <- paste(paste0(drug_list,"_ever_drug"), collapse = "|")
    followup_data <- pheno_baseline %>% rowwise() %>%
      dplyr::do({
            result = as_tibble(.) # result <- pheno_baseline %>% slice(49)
            x =  classify_drugs(result, phenotype, drug_list,low_inducers, high_inducers)
            result$high_inducer=x$case
            result$high_inducer_drug_num=x$drug_num
            result$high_inducer_drug_name=x$drug_name
            result$pheno_change=x$pheno_diff
            result$pheno_change_1mo=x$pheno_change_1mo
            result$pheno_change_3mo=x$pheno_change_3mo
            result$pheno_change_6mo=x$pheno_change_6mo
            result$follow_up_time=x$duration
            result$follow_up_time_1mo=x$duration_1mo
            result$follow_up_time_3mo=x$duration_3mo
            result$follow_up_time_6mo=x$duration_6mo
            result$age_started=x$age_started
            result$pheno_start=x$pheno_start
            result$pheno_slope=x$pheno_slope
            result$pheno_slope_6mo=x$pheno_slope_6mo
            result$pheno_slope_weight=x$pheno_slope_weight
            result$pheno_slope_weight_6mo=x$pheno_slope_weight_6mo
            result
        }) %>% ungroup() %>%
      mutate(ever_drug_match = rowSums(dplyr::select(., matches(col_match)) == 1) > 0) %>%
      filter(!(ever_drug_match & high_inducer==0)) %>% ##filter out individuals who have taken this drug but it wasn't followed
      mutate(follow_up_time_sq = follow_up_time^2) %>%
      mutate(follow_up_time_1mo_sq = follow_up_time_1mo^2) %>%
      mutate(follow_up_time_3mo_sq = follow_up_time_3mo^2) %>%
      mutate(follow_up_time_6mo_sq = follow_up_time_6mo^2) %>%
      mutate(age_sq=age_started^2)

    colnames(followup_data)[which(grepl("pheno_", colnames(followup_data)))] <- str_replace(colnames(followup_data)[which(grepl("pheno_", colnames(followup_data)))], "pheno_", paste0(phenotype, "_"))

    out[[paste0(drug_class, ".", phenotype)]] <- followup_data
  #}
  return(out)
}


classify_drug_cases <- function(x, phenotype, case_categories, preferential_control_categories, preferential_noncontrol_categories) {
  drug_info <- clean_drugs(x, phenotype)

  if(case_categories!="Other"){
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
  }

  if(case_categories=="Other"){
    case_match <- unlist(drug_info %>% filter(time >= min_follow_up) %>%
      filter(drug_num==min(drug_num)) %>% dplyr::select(drug_num))[1] ## preferred control
    control_match <- NA
  }

  case <- unname(ifelse (is.na(case_match), 0, 1))
  if(is.na(control_match) & is.na(case_match)) {case <- NA}
  best_match <- unname(ifelse (is.na(case_match), control_match, case_match))
  pheno_start <- dplyr::pull(drug_info[best_match,"pheno_start"])
  time_temp <- dplyr::pull(drug_info[best_match,"time"])
  time_1mo_temp <- dplyr::pull(drug_info[best_match,"time_1mo"])
  time_3mo_temp <- dplyr::pull(drug_info[best_match,"time_3mo"])
  time_6mo_temp <- dplyr::pull(drug_info[best_match,"time_6mo"])

  drug_temp <- dplyr::pull(drug_info[best_match,"drug"])
  pheno_change <- dplyr::pull(drug_info[best_match,"pheno_change"])
  pheno_change_1mo <- dplyr::pull(drug_info[best_match,"pheno_change_1mo"])
  pheno_change_3mo <- dplyr::pull(drug_info[best_match,"pheno_change_3mo"])
  pheno_change_6mo <- dplyr::pull(drug_info[best_match,"pheno_change_6mo"])

  pheno_slope <- dplyr::pull(drug_info[best_match,"pheno_slope"])
  pheno_slope_weight <- dplyr::pull(drug_info[best_match,"pheno_slope_weight"])
  pheno_slope_6mo <- dplyr::pull(drug_info[best_match,"pheno_slope_6mo"])
  pheno_slope_weight_6mo <- dplyr::pull(drug_info[best_match,"pheno_slope_weight_6mo"])

  age_started <- dplyr::pull(drug_info[best_match,"age_started"])
  return(list(case=case, drug_num=best_match, drug_name=drug_temp,age_started=age_started, pheno_start=pheno_start, pheno_diff=pheno_change,
    pheno_change_1mo=pheno_change_1mo, pheno_change_3mo=pheno_change_3mo, pheno_change_6mo=pheno_change_6mo,
    pheno_slope=pheno_slope, pheno_slope_weight=pheno_slope_weight, pheno_slope_6mo=pheno_slope_6mo, pheno_slope_weight_6mo=pheno_slope_weight_6mo,
    duration=time_temp, duration_1mo=time_1mo_temp, duration_3mo=time_3mo_temp, duration_6mo=time_6mo_temp))
}

munge_pheno_case_only <-  function(pheno_baseline, drug_prioritization, phenotype, low_inducers, high_inducers) {
  out <- list()
  for(i in 1:length(drug_prioritization))
  {
    drug <- drug_prioritization[i]
    #drug_list <- unlist(test_drugs %>% dplyr::select(drugs) %>% dplyr::slice(i))
    #drug_class <- unlist(test_drugs %>% dplyr::select(class) %>% dplyr::slice(i))
    col_match <- paste(paste0(drug,"_ever_drug"), collapse = "|")
    followup_data <- pheno_baseline %>% rowwise() %>%
      dplyr::do({
            result = as_tibble(.) # result <- pheno_baseline %>% slice(49)
            x =  classify_drug_cases(result, phenotype, drug,low_inducers, high_inducers)
            result$high_inducer=x$case
            result$high_inducer_drug_num=x$drug_num
            result$high_inducer_drug_name=x$drug_name
            result$pheno_change=x$pheno_diff
            result$pheno_change_1mo=x$pheno_change_1mo
            result$pheno_change_3mo=x$pheno_change_3mo
            result$pheno_change_6mo=x$pheno_change_6mo
            result$follow_up_time=x$duration
            result$follow_up_time_1mo=x$duration_1mo
            result$follow_up_time_3mo=x$duration_3mo
            result$follow_up_time_6mo=x$duration_6mo
            result$age_started=x$age_started
            result$pheno_start=x$pheno_start
            result$pheno_slope=x$pheno_slope
            result$pheno_slope_6mo=x$pheno_slope_6mo
            result$pheno_slope_weight=x$pheno_slope_weight
            result$pheno_slope_weight_6mo=x$pheno_slope_weight_6mo
            result
        }) %>% ungroup() %>%
      mutate(ever_drug_match = rowSums(dplyr::select(., matches(col_match)) == 1) > 0) %>%
      filter(!(ever_drug_match & high_inducer==0)) %>% ##filter out individuals who have taken this drug but it wasn't followed
      filter(high_inducer==1) %>%
      mutate(follow_up_time_sq = follow_up_time^2) %>%
      mutate(follow_up_time_1mo_sq = follow_up_time_1mo^2) %>%
      mutate(follow_up_time_3mo_sq = follow_up_time_3mo^2) %>%
      mutate(follow_up_time_6mo_sq = follow_up_time_6mo^2) %>%
      mutate(age_sq=age_started^2)

    colnames(followup_data)[which(grepl("pheno_", colnames(followup_data)))] <- str_replace(colnames(followup_data)[which(grepl("pheno_", colnames(followup_data)))], "pheno_", paste0(phenotype, "_"))

    out[[paste0(drug, ".", phenotype)]] <- followup_data

    pheno_baseline <- pheno_baseline %>%
      filter(!GPCR %in% followup_data$GPCR )

  }
  return(out)
}

create_GWAS_pheno <- function(pheno_baseline, pheno_followup, caffeine_vars, test_drugs, high_inducers, med_inducers, low_inducers, outcomes){
  na_to_none <- function(x) ifelse(is.na(x),'NONE',x)
  linear_pheno <- pheno_baseline %>% dplyr::select(matches('^FID$|^IID$|_start_Drug_1'), paste0(caffeine_vars, "_caffeine")) %>%
    mutate_at(vars(-FID,-IID), ivt)
  linear_covar <- pheno_baseline %>% dplyr::select(matches('^FID$|^IID$|^sex$|^Age_Drug_1$|Age_sq_Drug_1$|^PC[0-9]$|^PC[1][0-9]$|^PC20|^Age_caffeine'))

  list_covar <- list(linear=linear_covar)
  list_pheno <- list(linear=linear_pheno)
  followup_names <- names(pheno_followup)
  drug_names <- unique(sapply(followup_names, function(x) {unlist(strsplit(x, "[.]"))[1]}))
  phenotypes <- unique(sapply(followup_names, function(x) {unlist(strsplit(x, "[.]"))[2]}))

  #for(sub in names(pheno_followup))
  for(drug in drug_names){
    for(pheno in phenotypes){
      sub <- paste0(drug, ".", pheno)
      data_drug <- pheno_followup[[sub]]
      data_drug <- data_drug %>%
        mutate_at(.vars = paste0(pheno, outcomes), .funs = ivt) # perform IVT on full data set instead of before sensitivity phenotypes are created
                                                                # this means that the data used in the sensitivity models is exactly the same as the full models, just with some participants removed.

      drugs_analyzed <- test_drugs[which(test_drugs$class==drug),] %>% pull(drugs) %>% unlist()
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
        mutate_at(paste0(pheno, outcomes), ~ifelse(is.na(high_inducer), NA, .)) %>%
        dplyr::select(FID, IID, high_inducer, paste0(pheno, outcomes)) %>% rename_at(vars(-FID,-IID),function(x) paste0(x,"_sensitivity"))

      all_var <- data_drug %>% dplyr::select(matches('^FID$|^IID$|^high_inducer$|^age_|^follow_up_time',ignore.case = F), paste0(pheno, "_start"), paste0(pheno, outcomes)) %>%
        left_join(data_drug_sensitivity) %>%
        mutate_at(c("high_inducer", "high_inducer_sensitivity"), as.character) %>%
        mutate_at(vars(high_inducer, high_inducer_sensitivity),  ~(recode(.,"0" = "NoDrug", "1" = "Drug"))) %>%
        rename_at(vars(-FID,-IID),function(x) paste0(x,"_",drug)) %>%
        mutate_if(.predicate = is.character,.funs = na_to_none)

        # mutate_at(.vars = c(paste0(interaction_outcome,"_", sub), paste0(interaction_outcome,"_sensitivity_", sub)), .funs = ivt) ## if you want to perform IVT *after* sensitivity phenotypes are created
        # .funs= ~ifelse(abs(.)>mean(., na.rm=T)+3*sd(., na.rm=T), NA, .)) ## to remove outliers

      pheno_var <- all_var %>% dplyr::select(matches('^FID$|^IID$|^high_inducer',ignore.case = F), paste0(paste0(pheno, outcomes),"_", drug),
        paste0(paste0(pheno, outcomes),"_sensitivity_", drug))

      covar_var <- all_var %>% dplyr::select(matches('^FID$|^IID$|^age_|^high_inducer|^follow_up_time',ignore.case = F), paste0(pheno, "_start_", drug))
      list_pheno[[sub]] <- pheno_var
      list_covar[[sub]] <- covar_var
    }

  }

  full_pheno <- reduce(list_pheno, full_join) %>%
    mutate_if(.predicate = is.character,.funs = na_to_none)
  full_covar <- reduce(list_covar, full_join) %>%
    mutate_if(.predicate = is.character,.funs = na_to_none) %>%
    mutate_at(vars(starts_with("high_inducer_")), function(x) case_when(x == "NONE" ~ NA_real_, x == "Drug" ~ 1, x == "NoDrug" ~ 0))  %>%
    mutate_at("sex", function(x) case_when(x == "NONE" ~ NA_real_, x == "F" ~ 1, x == "M" ~ 0))
  out <- list(full_pheno = full_pheno, full_covar = full_covar)
  return(out)
}


create_GWAS_pheno_case_only <- function(pheno_baseline, pheno_followup, caffeine_vars, test_drugs, high_inducers, med_inducers, low_inducers, outcomes){
  na_to_none <- function(x) ifelse(is.na(x),'NONE',x)
  linear_pheno <- pheno_baseline %>% dplyr::select(matches('^FID$|^IID$|_start_Drug_1'), paste0(caffeine_vars, "_caffeine")) %>%
    mutate_at(vars(-FID,-IID), ivt)
  linear_covar <- pheno_baseline %>% dplyr::select(matches('^FID$|^IID$|^sex$|^Age_Drug_1$|Age_sq_Drug_1$|^PC[0-9]$|^PC[1][0-9]$|^PC20|^Age_caffeine'))

  list_covar <- list(linear=linear_covar)
  list_pheno <- list(linear=linear_pheno)
  followup_names <- names(pheno_followup)
  drug_names <- unique(sapply(followup_names, function(x) {unlist(strsplit(x, "[.]"))[1]}))
  phenotypes <- unique(sapply(followup_names, function(x) {unlist(strsplit(x, "[.]"))[2]}))

  #for(sub in names(pheno_followup))
  for(drug in drug_names){
    for(pheno in phenotypes){
      sub <- paste0(drug, ".", pheno)
      data_drug <- pheno_followup[[sub]]
      data_drug <- data_drug %>%
        mutate_at(.vars = paste0(pheno, outcomes), .funs = ivt) # perform IVT on full data set instead of before sensitivity phenotypes are created
                                                                # this means that the data used in the sensitivity models is exactly the same as the full models, just with some participants removed.

      all_var <- data_drug %>% dplyr::select(matches('^FID$|^IID$|^high_inducer$|^age_|^follow_up_time',ignore.case = F), paste0(pheno, "_start"), paste0(pheno, outcomes)) %>%
        mutate_at(c("high_inducer"), as.character) %>%
        mutate_at(vars(high_inducer),  ~(recode(.,"0" = "NoDrug", "1" = "Drug"))) %>%
        rename_at(vars(-FID,-IID),function(x) paste0(x,"_",drug)) %>%
        mutate_if(.predicate = is.character,.funs = na_to_none)

      pheno_var <- all_var %>% dplyr::select(matches('^FID$|^IID$|^high_inducer',ignore.case = F), paste0(paste0(pheno, outcomes),"_", drug))

      covar_var <- all_var %>% dplyr::select(matches('^FID$|^IID$|^age_|^high_inducer|^follow_up_time',ignore.case = F), paste0(pheno, "_start_", drug))
      list_pheno[[sub]] <- pheno_var
      list_covar[[sub]] <- covar_var
    }

  }

  full_pheno <- reduce(list_pheno, full_join) %>%
    mutate_if(.predicate = is.character,.funs = na_to_none)
  full_covar <- reduce(list_covar, full_join) %>%
    mutate_if(.predicate = is.character,.funs = na_to_none) %>%
    mutate_at(vars(starts_with("high_inducer_")), function(x) case_when(x == "NONE" ~ NA_real_, x == "Drug" ~ 1, x == "NoDrug" ~ 0))  %>%
    mutate_at("sex", function(x) case_when(x == "NONE" ~ NA_real_, x == "F" ~ 1, x == "M" ~ 0))
  out <- list(full_pheno = full_pheno, full_covar = full_covar)
  return(out)
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

residualize_pheno <- function(GWAS_input, pheno, covars, outcome_type, eths, eth_data){
  covars <- unlist(covars)
  pheno_name <- pheno

  if(outcome_type=="subgroup"){
    pheno <- unlist(pheno)
    pheno_name <- pheno[-which(grepl('^high_inducer_', pheno))]
    subgroup_name <- pheno[which(grepl('^high_inducer_', pheno))]

  }

  model_data <- GWAS_input$full_pheno %>% dplyr::select(FID, pheno) %>%
    left_join(GWAS_input$full_covar %>% dplyr::select(FID, IID, covars)) %>%
    rename("outcome" = pheno_name) %>%
    left_join(eth_data %>% dplyr::select(FID, eth)) %>%
    mutate(eth=as.factor(eth))


  if(outcome_type=="baseline"){
    resid_data <- model_data %>% group_by(eth) %>%
    do(  #
      {
        if(all(is.na(.$outcome))){
          resid_eth <- rep(NA, length(dim(.)[1]))
        } else resid_eth <- resid(lm(outcome ~ ., data = subset(., select=c( -FID, -IID, -eth) ), na.action = na.exclude))

        data.frame(., resid_eth)
      }) %>% ungroup() %>%

    dplyr::select(FID, IID, resid_eth) %>%
    rename(!!pheno := resid_eth)
  }

  if(outcome_type=="interaction"){
    high_inducer_name <- covars[which(grepl('^high_inducer_', covars))]
    resid_data <- model_data %>% rename("high_inducer" = high_inducer_name)  %>%
      group_by(eth) %>%
      do(
        { # x <- resid_data %>% filter(eth=="CEU")
          # data <- subset(x, select=c( -FID, -IID, -eth, -high_inducer) )
          # resid_eth <- resid(lm(outcome ~ ., data = data, na.action = na.exclude))
          data <- subset(., select=c( -FID, -IID, -eth, -high_inducer) )
          if(all(is.na(data$outcome))){
            resid_eth <- rep(NA, length(dim(data)[1]))
          } else resid_eth <- resid(lm(outcome ~ ., data = data, na.action = na.exclude))
          data.frame(., resid_eth)
        }) %>% ungroup() %>%

      dplyr::select(FID, IID, resid_eth, high_inducer) %>%
      rename(!!pheno := resid_eth) %>% rename(!!high_inducer_name := high_inducer)

  }


  if(outcome_type=="subgroup"){

    resid_data <- model_data %>%
      rename("subgroup" = subgroup_name) %>% mutate(subgroup=as.factor(subgroup)) %>%
      #filter(subgroup=="Drug" & eth=="YRI")
      #resid_group <- resid(lm(outcome ~ ., data = subset(resid_data, select=c( -FID, -IID, -eth, -subgroup) ), na.action = na.exclude))

      group_by(eth, subgroup) %>%
      do(
        {
          data <- subset(., select=c( -FID, -IID, -eth, -subgroup) )
          #print(as.character(.$subgroup[1]))
          #print(as.character(.$eth[1]))
          if(all(.$subgroup=="NONE") | all(!complete.cases(data)))
          {
            #print("YES")
            resid_group <- rep(NA, dim(.)[1])

          } else {
            data <- subset(., select=c( -FID, -IID, -eth, -subgroup) )
            if(all(is.na(data$outcome))){
              resid_group <- rep(NA, length(dim(data)[1]))
            } else resid_group <- resid(lm(outcome ~ ., data = data, na.action = na.exclude))

            #print("NO")
          }
          data.frame(., resid_group)
        }
      ) %>% ungroup() %>%

          # resid_eth <- tryCatch(resid(lm(outcome ~ ., data = subset(., select=c( -FID, -IID, -eth, -subgroup) ), na.action = na.exclude)),
          #   error=function(e){NA})

          ## some of these models give errors:
          ## check: table(resid_data$subgroup, resid_data$eth)
      dplyr::select(FID, IID, resid_group, subgroup) %>%
      rename(!!pheno_name := resid_group) %>% rename(!!subgroup_name := subgroup)

  }

  out <- list()
  out[[pheno_name]] <- resid_data
  return(out)


}

create_analysis_dirs <- function(top_level, eths){
  dir.create(top_level,showWarnings=F)

  to_create <- c(file.path(top_level, "interaction"), file.path(top_level, "full"), file.path(top_level, "subgroup"),  file.path(top_level, "case_only"))
  for(dir in to_create){
    dir.create(dir,showWarnings=F)
    for (eth in eths){
      dir.create(file.path(dir, eth),showWarnings=F)
    }
    dir.create(file.path(dir, "META"),showWarnings=F)
    dir.create(file.path(dir, "processed"),showWarnings=F)
  }
  dir.create("analysis/PRS",showWarnings=F)
  dir.create("analysis/GWAS/UKBB",showWarnings=F)
  dir.create("analysis/GWAS/case_only/META_DRUGS")

  return(TRUE)
}

define_baseline_models <- function(GWAS_input, baseline_vars, drug_classes, caffeine_vars, interaction_outcome){
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

        for(pheno_temp in baseline_vars){
          if(grepl(pheno_temp, outcome)){
            phenotype <- pheno_temp
          }
        }

        # covars <- c(covars,
        #   colnames((dplyr::select(GWAS_input$full_covar, ends_with(drug_classes[i])) %>%
        #     dplyr::select(matches('age|bmi_start_|follow_up_time_'))))
        #   )

        outcome_ending <- substrRight(outcome, 2)
        if(outcome_ending=="mo"){
          month_str <- substrRight(outcome, 3)
          covars <- c(covars,
            colnames((dplyr::select(GWAS_input$full_covar, ends_with(drug_classes[i]))) %>%
              dplyr::select(matches('age'), paste0(phenotype, "_start_", drug_classes[i])))
            )
          covars <- c(covars,
            colnames((dplyr::select(GWAS_input$full_covar, ends_with(drug_classes[i]))) %>%
              dplyr::select(matches('follow_up_time_')) %>%
              dplyr::select(contains(month_str)))
            )

        }
        if(!outcome_ending=="mo"){
          covars <- c(covars,
            colnames((dplyr::select(GWAS_input$full_covar, ends_with(drug_classes[i])) %>%
              dplyr::select(matches('age|follow_up_time_'), paste0(phenotype, "_start_", drug_classes[i]))))
            )

          covars <- covars[-which(grepl("_[1-9]mo_", covars))]

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

define_interaction_models <- function(GWAS_input, drug_classes, interaction_outcome, baseline_vars){
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
        interactoin_pam_temp <- "1-3"
        interaction_vars <- c(interaction_vars,
          colnames(GWAS_input$full_pheno)[which(colnames(GWAS_input$full_pheno)==paste0(outcome,analysis,drug_classes[i]))])
        covars <- c(standard_covars)

        for(pheno_temp in baseline_vars){
          if(grepl(pheno_temp, outcome)){
            phenotype <- pheno_temp
          }
        }

        outcome_ending <- substrRight(outcome, 2)
        if(outcome_ending=="mo"){
          month_str <- substrRight(outcome, 3)
          covars <- c(covars,
            colnames((dplyr::select(GWAS_input$full_covar, ends_with(drug_classes[i]))) %>%
              dplyr::select(matches('age|high_inducer_'), paste0(phenotype, "_start_", drug_classes[i])))
            )
          covars <- c(covars,
            colnames((dplyr::select(GWAS_input$full_covar, ends_with(drug_classes[i]))) %>%
              dplyr::select(matches('follow_up_time_')) %>%
              dplyr::select(contains(month_str)))
            )

        }
        if(!outcome_ending=="mo"){
          covars <- c(covars,
            colnames((dplyr::select(GWAS_input$full_covar, ends_with(drug_classes[i])) %>%
              dplyr::select(matches('age|high_inducer_|follow_up_time_'), paste0(phenotype, "_start_", drug_classes[i]))))
            )

          covars <- covars[-which(grepl("_[1-9]mo_", covars))]

        }


        high_inducer_sens <- which(grepl('^high_inducer_sensitivity_', covars))
        high_inducer_reg <- which(covars == paste0("high_inducer_", drug_classes[i]))


        if(analysis=="_"){
          covars <- covars[-high_inducer_sens]
        }
        if(analysis=="_sensitivity_"){
          interactoin_pam_temp <- "1-3"
          covars[high_inducer_reg] <- covars[high_inducer_sens]
          covars <- covars[-high_inducer_sens]
        }


        interaction_covars[[length(interaction_covars)+1]]  <- covars
        interaction_pams <- c(interaction_pams, interactoin_pam_temp)
        output_suffix <- c(output_suffix, paste0(interaction_outcome[j], analysis, drug_classes[i] ))
        drug <- c(drug, substring(paste0(analysis, drug_classes[i]), 2))

      }

    }


  }

  interaction_gwas_info <- tibble ( pheno = interaction_vars, ## this is y
                                    covars = interaction_covars, ## these are the x's
                                    parameters = interaction_pams, ## this is passed to PLINK
                                    output_suffix = output_suffix, ## this will be the output name
                                    drug = drug) ## this is the name used to grab interaction column
  return(interaction_gwas_info)
}

define_subgroup_models <- function(GWAS_input, drug_classes, interaction_outcome, baseline_vars){
  subgroup_vars <- numeric()
  subgroup_covars <- list()
  subgroup_pheno <- list()
  output_suffix <- numeric()
  for(i in 1:length(drug_classes)){
    for(j in 1:length(interaction_outcome)){
      for(analysis in c("_", "_sensitivity_"))
      {
        outcome <- interaction_outcome[j]
        pheno <- c(paste0(outcome, analysis,drug_classes[i]),paste0("high_inducer",analysis,drug_classes[i]))
        subgroup_vars <- c(subgroup_vars,
          colnames(GWAS_input$full_pheno)[which(colnames(GWAS_input$full_pheno)==paste0("high_inducer",analysis,drug_classes[i]))])
        covars <- c(standard_covars)

        for(pheno_temp in baseline_vars){
          if(grepl(pheno_temp, outcome)){
            phenotype <- pheno_temp
          }
        }

        outcome_ending <- substrRight(outcome, 2)
        if(outcome_ending=="mo"){
          month_str <- substrRight(outcome, 3)
          covars <- c(covars,
            colnames((dplyr::select(GWAS_input$full_covar, ends_with(drug_classes[i]))) %>%
              dplyr::select(matches('age'), paste0(phenotype, "_start_", drug_classes[i])))
            )
          covars <- c(covars,
            colnames((dplyr::select(GWAS_input$full_covar, ends_with(drug_classes[i]))) %>%
              dplyr::select(matches('follow_up_time_')) %>%
              dplyr::select(contains(month_str)))
            )

        }
        if(!outcome_ending=="mo"){
          covars <- c(covars,
            colnames((dplyr::select(GWAS_input$full_covar, ends_with(drug_classes[i])) %>%
              dplyr::select(matches('age|follow_up_time_'), paste0(phenotype, "_start_", drug_classes[i]))))
            )

          covars <- covars[-which(grepl("_[1-9]mo_", covars))]

        }


        subgroup_covars[[length(subgroup_covars)+1]]  <- covars
        subgroup_pheno[[length(subgroup_pheno)+1]]  <- pheno
        output_suffix <- c(output_suffix, paste0(interaction_outcome[j], analysis, drug_classes[i]))
      }

    }

  }
  subgroup_gwas_info <- tibble( pheno = subgroup_pheno,
                                subgroup = subgroup_vars,
                                covars = subgroup_covars,
                                output_suffix = output_suffix)
  return(subgroup_gwas_info)
}

define_case_only_models <- function(GWAS_input, drug_prioritization, interaction_outcome, baseline_vars){
  case_only_covars <- list()
  case_only_pheno <- numeric()
  output_suffix <- numeric()
  for(i in 1:length(drug_prioritization)){
    for(j in 1:length(interaction_outcome)){
      #for(analysis in c("_", "_sensitivity_"))
      #{
        analysis <- "_"
        outcome <- interaction_outcome[j]
        pheno <- paste0(outcome, analysis,drug_prioritization[i])
        covars <- c(standard_covars)

        for(pheno_temp in baseline_vars){
          if(grepl(pheno_temp, outcome)){
            phenotype <- pheno_temp
          }
        }

        outcome_ending <- substrRight(outcome, 2)
        if(outcome_ending=="mo"){
          month_str <- substrRight(outcome, 3)
          covars <- c(covars,
            colnames((dplyr::select(GWAS_input$full_covar, ends_with(drug_prioritization[i]))) %>%
              dplyr::select(matches('age'), paste0(phenotype, "_start_", drug_prioritization[i])))
            )
          covars <- c(covars,
            colnames((dplyr::select(GWAS_input$full_covar, ends_with(drug_prioritization[i]))) %>%
              dplyr::select(matches('follow_up_time_')) %>%
              dplyr::select(contains(month_str)))
            )

        }
        if(!outcome_ending=="mo"){
          covars <- c(covars,
            colnames((dplyr::select(GWAS_input$full_covar, ends_with(drug_prioritization[i])) %>%
              dplyr::select(matches('age|follow_up_time_'), paste0(phenotype, "_start_", drug_prioritization[i]))))
            )

          covars <- covars[-which(grepl("_[1-9]mo_", covars))]

        }


        case_only_covars[[length(case_only_covars)+1]]  <- covars
        case_only_pheno[[length(case_only_pheno)+1]]  <- pheno
        output_suffix <- c(output_suffix, paste0(interaction_outcome[j], analysis, drug_prioritization[i]))
      #}

    }

  }
  case_only_gwas_info <- tibble( pheno = case_only_pheno,
                                covars = case_only_covars,
                                output_suffix = output_suffix)
  return(case_only_gwas_info)
}

extract_colnames <- function(data, pattern, additional_columns = c("FID", "IID")){

  columns <- colnames(data)
  columns_extract <- columns[which(grepl(pattern, columns))]
  data_out <- data[,c(additional_columns, columns_extract)]
  return(data_out)
}


residualize_pheno_case_only <- function(GWAS_input, pheno, covars, outcome_type, eths, eth_data){
  covars <- unlist(covars)
  pheno_name <- pheno

  model_data <- GWAS_input$full_pheno %>% dplyr::select(FID, pheno) %>%
    left_join(GWAS_input$full_covar %>% dplyr::select(FID, IID, covars)) %>%
    rename("outcome" = pheno_name) %>%
    left_join(eth_data %>% dplyr::select(FID, eth)) %>%
    mutate(eth=as.factor(eth))


  if(outcome_type=="case_only"){
    resid_data <- model_data %>% group_by(eth) %>%
    do(  #
      {
        if(all(is.na(.$outcome))){
          resid_eth <- rep(NA, length(dim(.)[1]))
        } else resid_eth <- resid(lm(outcome ~ ., data = subset(., select=c( -FID, -IID, -eth) ), na.action = na.exclude))

        data.frame(., resid_eth)
      }) %>% ungroup() %>%

    dplyr::select(FID, IID, resid_eth) %>%
    rename(!!pheno := resid_eth)
  }

  out <- list()
  out[[pheno_name]] <- resid_data
  return(out)

}


write_interaction_resid <- function(interaction_resid_data, drug_name){

  pattern <- paste0("_", drug_name,  "$")
  pattern_sensitivity <- paste0("_sensitivity_", drug_name,  "$")
  file_name_list <- numeric()

  for(pattern_test in c(pattern_sensitivity, pattern)){

    data_out <- extract_colnames(interaction_resid_data, pattern_test)
    if(!grepl("sensitivity", pattern_test, fixed = TRUE)){
      data_out <- data_out[,-which(grepl("sensitivity", colnames(data_out)))]
    }
    drug_covar_out <- data_out[,c(1,2,which(grepl("high_inducer_", colnames(data_out))))]
    data_out <- data_out[,-which(grepl("high_inducer_", colnames(data_out)))]
    drug_file <- ifelse(grepl("sensitivity", pattern_test, fixed = TRUE), paste0("sensitivity_", drug_name), drug_name)
    file_name <- paste0("data/processed/phenotype_data/GWAS_input/interaction_input_", drug_file, "_resid.txt")
    covar_file_name <- paste0("data/processed/phenotype_data/GWAS_input/interaction_covar_input_", drug_file, "_resid.txt")

    write.table(data_out, file_name, row.names = F, quote = F, col.names = T)
    write.table(drug_covar_out, covar_file_name, row.names = F, quote = F, col.names = T)

    file_name_list <- c(file_name_list, file_name, covar_file_name)

  }

  return(file_name_list)

}

write_subgroup_resid <- function(subgroup_resid_data, drug_name){

  pattern <- paste0("_", drug_name,  "$")
  pattern_sensitivity <- paste0("_sensitivity_", drug_name,  "$")

  file_name_list <- numeric()
  for(pattern_test in c(pattern_sensitivity, pattern)){
    data_out <- extract_colnames(subgroup_resid_data, pattern_test)
    if(!grepl("sensitivity", pattern_test, fixed = TRUE)){
      data_out <- data_out[,-which(grepl("sensitivity", colnames(data_out)))]
    }
    drug_file <- ifelse(grepl("sensitivity", pattern_test, fixed = TRUE), paste0("sensitivity_", drug_name), drug_name)
    file_name <- paste0("data/processed/phenotype_data/GWAS_input/subgroup_input_", drug_file, "_resid.txt")
    write.table(data_out, file_name, row.names = F, quote = F, col.names = T)
    file_name_list <- c(file_name_list, file_name)

  }
  return(file_name_list)

}


define_interaction_inputs <- function(drug_classes){

  interaction_pheno_file <- numeric()
  interaction_covar_file <- numeric()
  interaction_covars <- list()
  interaction_pams <- numeric()
  output_suffix <- numeric()
  drug <- numeric()
  for(i in 1:length(drug_classes)){
    drug_name <- drug_classes[i]
    for(analysis in c("_", "_sensitivity_"))
    {
      drug_file <- ifelse(analysis=="_sensitivity_", paste0("sensitivity_", drug_name), drug_name)
      file_name <- paste0("data/processed/phenotype_data/GWAS_input/interaction_input_", drug_file, "_resid.txt")
      covar_file_name <- paste0("data/processed/phenotype_data/GWAS_input/interaction_covar_input_", drug_file, "_resid.txt")
      covars <- paste0("high_inducer", analysis, drug_name)
      interactoin_pam_temp <- "1-3"

      interaction_pheno_file <- c(interaction_pheno_file, file_name)
      interaction_covar_file <- c(interaction_covar_file, covar_file_name)

      interaction_covars[[length(interaction_covars)+1]]  <- covars
      interaction_pams <- c(interaction_pams, interactoin_pam_temp)
      output_suffix <- c(output_suffix, gsub("^.", "", paste0(analysis, drug_name)))
    }

  }


  interaction_gwas_inputs <- tibble( pheno_file = interaction_pheno_file, ## this is pheno file
                                     covar_file = interaction_covar_file, ## this is the covar file
                                     covars = interaction_covars, ## these are the x'x
                                     parameters = interaction_pams, ## this is passed to PLINK
                                     output_suffix = output_suffix, ## this will be the output name
                                    ) ## this is the name used to grab interaction column
  return(interaction_gwas_inputs)

}

define_subgroup_inputs <- function(drug_classes){

  subgroup_pheno_file <- numeric()
  subgroup_vars <- numeric()
  subgroup_pheno <- list()
  output_suffix <- numeric()
  for(i in 1:length(drug_classes)){
    drug_name <- drug_classes[i]
    for(analysis in c("_", "_sensitivity_"))
    {
      drug_file <- ifelse(analysis=="_sensitivity_", paste0("sensitivity_", drug_name), drug_name)
      file_name <- paste0("data/processed/phenotype_data/GWAS_input/subgroup_input_", drug_file, "_resid.txt")

      subgroup_pheno_file <- c(subgroup_pheno_file, file_name)
      subgroup_vars <- c(subgroup_vars, paste0("high_inducer", analysis, drug_name))
      output_suffix <- c(output_suffix, gsub("^.", "", paste0(analysis, drug_name)))
    }

  }

  subgroup_gwas_inputs <- tibble( pheno_file = subgroup_pheno_file,
                                subgroup = subgroup_vars,
                                output_suffix = output_suffix)
  return(subgroup_gwas_inputs)
}

run_gwas <- function(pfile, pheno_file, pheno_name=NULL, covar_file=NULL, covar_names=NULL, eths,
  eth_sample_file, output_dir, output,
  remove_sample_file=NULL,threads=1,type="full",
  subgroup=NULL, subgroup_var="", interaction=FALSE,parameters=NULL, output_suffix="",
  freq_file="analysis/QC/15_final_processing/CEU/PSYMETAB_GWAS.CEU.afreq",
  eth_low_maf_file){

  # eth_sample_file : in the form of "keep_this_ETH.txt"
  # default option: --variance-standardize
  file_name <- output
  if(output_suffix!=""){
    file_name <- paste0(file_name,"_",output_suffix)
  }

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
    analysis_commands <- c("--glm", "interaction", "skip", "--parameters", parameters)
    file_name <- paste0(file_name,"_int")
  } else analysis_commands <- c("--glm", "hide-covar", "skip")

  if(!is.null(covar_names)){
    covar_commands <- unlist(c("--covar", covar_file, "--covar-name", unlist(covar_names)))
  } else covar_commands <- NULL

  if(!is.null(pheno_name)){
    pheno_commands <- unlist(c("--pheno-name", pheno_name))
  } else pheno_commands <- NULL

  general_commands <- unlist(c("--pfile", pfile, "--read-freq", freq_file,
          "--threads", threads, "--variance-standardize", "--pheno", pheno_file))
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

      plink_input <- c(general_commands, analysis_commands, remove_commands, subgroup_commands,eth_commands,covar_commands,pheno_commands)
      temp_out <- processx::run(command="plink2",plink_input, error_on_status=F)
      out[[eth]] <- temp_out
    }
  }
  return(out)
  #return(paste(plink_input, collapse="|"))
}

# run_gwas_freq_counts is the same function as run_gwas but with frequency observations and counts recorded using the 'cols' argument,
# see here: https://www.cog-genomics.org/plink/2.0/general_usage#colset and https://www.cog-genomics.org/plink/2.0/formats#glm_linear

run_gwas_freq_counts <- function(pfile, pheno_file, pheno_name=NULL, covar_file=NULL, covar_names=NULL, eths,
  eth_sample_file, output_dir, output,
  remove_sample_file=NULL,threads=1,type="full",
  subgroup=NULL, subgroup_var="", interaction=FALSE,parameters=NULL, output_suffix="",
  freq_file="analysis/QC/15_final_processing/CEU/PSYMETAB_GWAS.CEU.afreq",
  eth_low_maf_file){

  # eth_sample_file : in the form of "keep_this_ETH.txt"
  # default option: --variance-standardize
  file_name <- output
  if(output_suffix!=""){
    file_name <- paste0(file_name,"_",output_suffix)
  }

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
    analysis_commands <- c("--glm", "cols=+a1count,+a1freq", "interaction", "skip", "--parameters", parameters)
    file_name <- paste0(file_name,"_int")
  } else analysis_commands <- c("--glm", "cols=+a1count,+a1freq", "hide-covar", "skip")

  if(!is.null(covar_names)){
    covar_commands <- unlist(c("--covar", covar_file, "--covar-name", unlist(covar_names)))
  } else covar_commands <- NULL

  if(!is.null(pheno_name)){
    pheno_commands <- unlist(c("--pheno-name", pheno_name))
  } else pheno_commands <- NULL

  general_commands <- unlist(c("--pfile", pfile, "--read-freq", freq_file,
          "--threads", threads, "--variance-standardize", "--pheno", pheno_file))
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

      plink_input <- c(general_commands, analysis_commands, remove_commands, subgroup_commands,eth_commands,covar_commands,pheno_commands)
      temp_out <- processx::run(command="plink2",plink_input, error_on_status=F)
      out[[eth]] <- temp_out
    }
  }
  return(out)
  #return(paste(plink_input, collapse="|"))
}

meta <- function(output, output_suffix="", eths, pheno_list, output_dir = "analysis/GWAS", type = "full",
  threads=1, interaction_var = NA){

  write_dir <- file.path(output_dir, type, "META")
  file_name <- output
  if(output_suffix!=""){
    file_name <- paste0(file_name,"_",output_suffix)
  }

  file_name_drug <- NA
  file_name_nodrug <- NA
  if(type=="interaction")
  {
    file_name <- paste0(file_name,"_int")
    awk_col7 <- paste0("'NR==1 || $7 == ", '"ADDxhigh_inducer_', interaction_var, '" { print $0 }', "'")

  }

  out <- list()

  for(pheno in pheno_list){
    gwas_results <- numeric()
    gwas_results_drug <- numeric()
    gwas_results_nodrug <- numeric()
    for(eth in eths){

      file_name_eth <- paste0(file_name,"_",eth)
      eth_dir <- file.path(output_dir, type)

      if(type=="subgroup"){
        file_name_eth_drug <- paste0(file_name_eth, ".Drug")
        file_name_eth_nodrug <- paste0(file_name_eth, ".NoDrug")
        eth_output_drug <- file.path(eth_dir,eth,paste0(file_name_eth_drug, ".", pheno, "_", output_suffix, ".glm.linear"))
        eth_output_nodrug <- file.path(eth_dir,eth,paste0(file_name_eth_nodrug, ".", pheno, "_", output_suffix, ".glm.linear"))

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

        eth_output <- file.path(eth_dir,eth,paste0(file_name_eth, ".", pheno, "_", output_suffix, ".glm.linear"))
        eth_subset_output <- file.path(eth_dir,eth,paste0(file_name_eth, ".", pheno, "_", output_suffix, ".glm.linear.interaction"))

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

    full_output <- file.path(write_dir,paste0(file_name, ".", pheno))
    full_output_drug <- file.path(write_dir,paste0(file_name_drug, ".", pheno))
    full_output_nodrug <- file.path(write_dir,paste0(file_name_nodrug, ".", pheno))

    if(type == "interaction" | type == "full"){
      out[[paste0("meta_out_", pheno)]] <- processx::run(command="plink",c( "--meta-analysis", gwas_results,  "+", "qt", "report-all", "no-map",
        "--meta-analysis-snp-field", "ID", "--out", full_output, "--threads", threads), error_on_status=F)
      if(length(gwas_results)!=0){
        out[[paste0("meta_out_", pheno, "_compress")]] <- processx::run(command="gzip", gwas_results)
      }

    }

    if(type=="subgroup"){
      out[[paste0("meta_out_drug_", pheno)]] <- processx::run(command="plink",c( "--meta-analysis", gwas_results_drug,  "+", "qt", "report-all", "no-map",
        "--meta-analysis-snp-field", "ID", "--out", full_output_drug, "--threads", threads), error_on_status=F)
      out[[paste0("meta_out_nodrug_", pheno)]] <- processx::run(command="plink",c( "--meta-analysis", gwas_results_nodrug,  "+", "qt", "report-all", "no-map",
        "--meta-analysis-snp-field", "ID", "--out", full_output_nodrug, "--threads", threads), error_on_status=F)
      if(length(gwas_results_drug)!=0){
        out[[paste0("meta_out_drug_", pheno, "_compress")]] <- processx::run(command="gzip", gwas_results_drug)
      }
      if(length(gwas_results_nodrug)!=0){
        out[[paste0("meta_out_nodrug_", pheno, "_compress")]] <- processx::run(command="gzip", gwas_results_nodrug)
      }
    }

  }

  return(out)

}

meta_case_only_eths <- function(output, output_suffix="", eths, pheno, output_dir = "analysis/GWAS", type = "case_only",
  threads=1){

    write_dir <- file.path(output_dir, type, "META")
    file_name <- output
    if(output_suffix!=""){
      file_name <- paste0(file_name,"_",output_suffix)
    }

    out <- list()

    gwas_results <- numeric()

    for(eth in eths){

      file_name_eth <- paste0(file_name,"_",eth)
      eth_dir <- file.path(output_dir, type)

      eth_output <- file.path(eth_dir,eth,paste0(file_name_eth, ".", pheno, ".glm.linear"))

      if(file.exists(eth_output)){
        gwas_results <- c(gwas_results, eth_output)
      }


    }

    full_output <- file.path(write_dir,paste0(file_name, ".", pheno))


    out[[paste0("meta_out_", pheno)]] <- processx::run(command="plink",c( "--meta-analysis", gwas_results,  "+", "qt", "report-all", "no-map",
      "--meta-analysis-snp-field", "ID", "--out", full_output, "--threads", threads), error_on_status=F)
    # if(length(gwas_results)!=0){
    #   out[[paste0("meta_out_", pheno, "_compress")]] <- processx::run(command="gzip", gwas_results)
    # }

    return(out)


}

meta_case_only_drugs <- function(output, output_suffix="", eths, drug_groups, outcome, output_dir = "analysis/GWAS", type = "case_only",
  threads=1, interaction_var = NA){

  write_dir <- file.path(output_dir, type, "META_DRUGS")
  file_name <- output
  if(output_suffix!=""){
    file_name <- paste0(file_name,"_",output_suffix)
  }

  out <- list()

  for(eth in eths){
    gwas_results <- numeric()

    for(drug in drug_groups){

      file_name_eth <- paste0(file_name,"_",eth)
      eth_dir <- file.path(output_dir, type)

      outcome_output <- file.path(eth_dir,eth,paste0(file_name_eth, ".", outcome, "_", drug, ".glm.linear"))

      if(file.exists(outcome_output)){
        gwas_results <- c(gwas_results, outcome_output)
      }

    }


    full_output <- file.path(write_dir,paste0(file_name, "_", eth, ".", outcome))


    out[[paste0("meta_out_", eth, "_", outcome)]] <- processx::run(command="plink",c( "--meta-analysis", gwas_results,  "+", "qt", "report-all", "no-map",
      "--meta-analysis-snp-field", "ID", "--out", full_output, "--threads", threads), error_on_status=F)
    if(length(gwas_results)!=0){
      out[[paste0("meta_out_", outcome, "_", eth, "_compress")]] <- processx::run(command="gzip", gwas_results)

    }

  }

  return(out)


}

munge_sqc <- function(sqc,fam){
  sqc_sub <- sqc[,c(26:65)]
  colnames(sqc_sub) <- c(sapply(1:40, function(x) {paste0("PC_",x)}))
  sqc_sub$ID <- fam[,1]
  return(sqc_sub)
}

get_british_ids <- function(qc_data, fam){
  qc_data$ID <- fam[,1]
  IDs <- qc_data %>% filter(in_white_british_ancestry_subset==1) %>% pull(ID)
  return(IDs)

}

bmi_followup_exists <- function(...){
  x <- c(...)
  out <- FALSE
  if(length(which(!is.na(x)))==1){
    out <- FALSE
  }
  if(length(which(!is.na(x)))>1){
    out <- TRUE
  }
  return(out)

}

ukbb_slope <- function(date1, date2, date3, val1, val2, val3){

  dates <- c(date1, date2, date3)
  vals <- c(val1, val2, val3)

  complete_indices <- union(which(!is.na(dates)), which(!is.na(vals)))
  dates <- dates[complete_indices]
  vals <- vals[complete_indices]

  date_difference <- numeric()
  for(i in 1:length(dates)){
    date <- dates[i]
    temp <- as.numeric(difftime(date, dates[1], units="days"))
    date_difference <- c(date_difference, temp)
  }


  #beta <- NA
  beta <- summary(lm(vals~date_difference))$coefficients[2,1]
  return(beta)
}

ukbb_weighted_slope <- function(date1, date2, date3, val1, val2, val3){

  dates <- c(date1, date2, date3)
  vals <- c(val1, val2, val3)

  complete_indices <- union(which(!is.na(dates)), which(!is.na(vals)))
  dates <- dates[complete_indices]
  vals <- vals[complete_indices]

  date_difference <- numeric()
  for(i in 1:length(dates)){
    date <- dates[i]
    temp <- as.numeric(difftime(date, dates[1], units="days"))
    date_difference <- c(date_difference, temp)
  }

  #beta <- NA
  beta <- summary(lm(vals~date_difference))$coefficients[2,1]
  se <- summary(lm(vals~date_difference))$coefficients[2,2]

  std_beta <- beta/se
  return(std_beta)

}



munge_ukbb_drug_users <- function(ukbb_org, drug_codes_interest, medication_columns){

  data_munge <- ukbb_org %>%
    cbind(mapply(classify_ukbb_drugs, drug_codes_interest$codes, list(.[,medication_columns])) %>% set_colnames(paste0(drug_codes_interest$molecule, "_users"))) %>%
    dplyr::select(eid, paste0(drug_codes_interest$molecule, "_users"))

}

count_ukbb_drug_users <- function(ukbb_munge_drug_users, drug_codes_interest){
  data <- ukbb_munge_drug_users %>% dplyr::select(paste0(drug_codes_interest$molecule, "_users")) %>% map_df(~length(which(.==1))) %>% t()

  return(data)
}

create_UKBB_v2_snp_list <- function(UKBB_processed){
  chr_num = tibble(chr = 1:22)
  output <- chr_num %>% mutate(v2_snp_list_files = paste0(UKBB_processed, "v2_snp_list/", "snp_list_chr", chr, ".txt"))
  return(output)
}

munge_ukbb_bmi_slope <- function(ukbb_org, ukb_key, date_followup, bmi_var){

  cols_date <- dplyr::pull(ukb_key[which(ukb_key[,1] == date_followup),"col.name"])
  cols_pheno <-  dplyr::pull(ukb_key[which(ukb_key[,1] == bmi_var),"col.name"])
  #cols_field =  dplyr::pull(ukb_key[which(ukb_key[,1] == !!bmi_var),"field.tab"]),
  #cols_sex <- dplyr::pull(ukb_key[which(ukb_key[,1] == sex_var),"col.name"])
  #cols_age <- dplyr::pull(ukb_key[which(ukb_key[,1] == age_var),"col.name"])
  ## GET AGE COLUMN
  cols_interest <- c(cols_date, cols_pheno)
  ukbb_sub <- ukbb_org %>% dplyr::select(eid, cols_interest) %>%
    rowwise %>% mutate(bmi_followup_value = bmi_followup_exists(c(!!!syms(cols_pheno)))) %>%
    filter(bmi_followup_value==TRUE)

  ukbb_sub_slope <- ukbb_sub %>%
    mutate(bmi_slope = ukbb_slope(date_of_attending_assessment_centre_f53_0_0, date_of_attending_assessment_centre_f53_1_0, date_of_attending_assessment_centre_f53_2_0,
        body_mass_index_bmi_f21001_0_0, body_mass_index_bmi_f21001_1_0, body_mass_index_bmi_f21001_2_0)) %>%
    mutate(bmi_wt_slope = ukbb_weighted_slope(date_of_attending_assessment_centre_f53_0_0, date_of_attending_assessment_centre_f53_1_0, date_of_attending_assessment_centre_f53_2_0,
        body_mass_index_bmi_f21001_0_0, body_mass_index_bmi_f21001_1_0, body_mass_index_bmi_f21001_2_0))


  output <- list()

  for(col in c("bmi_slope", "bmi_wt_slope")){
    out_i <- ukbb_sub_slope %>% dplyr::select(eid, !!sym(col)) %>% filter(!is.na(!!sym(col)))
    output[[paste(col)]] <- out_i
  }

  return(output)

}

munge_ukbb_drug_users_bmi <- function(ukbb_org, ukb_key, ukbb_munge_drug_users, bmi_var){
  cols_pheno <-  dplyr::pull(ukb_key[which(ukb_key[,1] == bmi_var),"col.name"])
  #cols_sex <- dplyr::pull(ukb_key[which(ukb_key[,1] == sex_var),"col.name"])
  #cols_age <- dplyr::pull(ukb_key[which(ukb_key[,1] == age_var),"col.name"])

  ukbb_join <- full_join(ukbb_org %>% dplyr::select(eid, cols_pheno[1]), ukbb_munge_drug_users, by="eid")
  data_out <- ukbb_join %>%

    mutate_at(vars(ends_with("_users")), list(BMI = ~ case_when(
                                                . == 1 ~ !!!syms(cols_pheno[1]),
                                                . == 0 ~ NA_real_ )))

  BMI_cols <- data_out %>% dplyr::select(ends_with("_BMI")) %>% colnames()
  #cols_interest <- c(BMI_cols, cols_sex, cols_age)

  data_BMI <- data_out %>% dplyr::select(eid, BMI_cols)

  output <- list()

  for(col in BMI_cols){
    out_i <- data_BMI %>% dplyr::select(eid, !!sym(col)) %>% filter(!is.na(!!sym(col))) %>% as_tibble
    output[[paste(col)]] <- out_i
  }

  return(output)
}

filter_ukbb <- function(ukb_with_data_df, relatives, exclusion_list, british_subset){

  # ukb_with_data_df should be a dataframe with only 2 columns: (1) eid and (2) variable of interest
  data <- ukb_with_data_df[[1]]
  related_IDs_remove <- ukb_gen_samples_to_remove(relatives, ukb_with_data = as.integer(data$eid))

  ukbb_filter <- data %>% filter(!eid %in% !!related_IDs_remove) %>% filter(!eid %in% !!exclusion_list) %>%
    filter(eid %in% !!british_subset)

  out <- list()
  out[[colnames(data)[2]]] <- ukbb_filter
  return(out)
}

resid_ukbb <- function(ukbb_filter, ukbb_org, ukb_key, sqc_munge, sex_var, age_var){

  data <- ukbb_filter[[1]]
  outcome <- colnames(data)[2]

  ## join with genetic data
  ukbb_merge <- merge(data, ukbb_org, by="eid")

  ukbb_plus_geno_covars <- merge(ukbb_merge, sqc_munge, by.x="eid", by.y="ID")
  #cols_pheno <-  dplyr::pull(ukb_key[which(ukb_key[,1] == bmi_var),"col.name"])
  cols_sex <- dplyr::pull(ukb_key[which(ukb_key[,1] == sex_var),"col.name"])
  cols_age <- dplyr::pull(ukb_key[which(ukb_key[,1] == age_var),"col.name"])

  #cols_interest <- c("eid", outcome, cols_pheno[1], cols_sex, cols_age, colnames(ukbb_merge)[grepl("PC", colnames(ukbb_merge))])
  cols_interest <- c("eid", outcome, cols_sex, cols_age, colnames(ukbb_merge)[grepl("PC", colnames(ukbb_merge))])

  data_sub <- ukbb_plus_geno_covars[, cols_interest] %>% rename(outcome = outcome)
  data_sub <- data_sub %>% mutate(outcome_ivt = ivt(outcome))
  resid_data <- resid(lm(outcome_ivt ~ ., data = subset(data_sub, select=c( -eid, -outcome) ), na.action = na.exclude))
  out <- as.data.frame(cbind(data_sub[["eid"]], resid_data))
  colnames(out) <- c("eid", paste0(outcome, "_resid"))

  out_list <- list()
  out_list[[paste0(outcome, "_resid")]] <- out
  return(out_list)

}


# define function to classify drugs into based on list of IDs
classify_ukbb_drugs <- function(risk_codes, ...){

  t <- cbind(...)
  #return(t)
  return(apply(t, 1, function(x) ifelse(any(x %in% risk_codes), 1, 0)))
}


order_bgen <- function(bgen_file, data){

  ord = match(bgen_file$ID_1[-1], data$eid)

  bgen_merge <- left_join(bgen_file %>% slice(-1), data, by=c("ID_1" = "eid"))

  add_na <- bgen_merge %>% mutate_at(vars(-ID_1, -ID_2, -missing), ~replace_na(., -999))

  output <- add_na %>% dplyr::select(-ID_1, -ID_2, -missing)
  return(output)
}

make_ukbb_chunks <- function(v2_snp_list_file, chunk_size=1e6){

  v2_snp_list <- fread(v2_snp_list_file, data.table=F)
  min_positon <- min(v2_snp_list$position)
  max_position <- max(v2_snp_list$position)
  chr <- as.character(v2_snp_list$chr[1])
  if(nchar(as.character(v2_snp_list$chr[1]))==1){
    chr_char <- paste0("0", as.character(v2_snp_list$chr[1])) } else chr_char <- as.character(v2_snp_list$chr[1])

  out <- numeric()

  chunk_starts <- seq(min_positon, max_position, chunk_size)
  for(i in 1:length(chunk_starts)){
    chunk_num <- i
    start <- chunk_starts[i]
    if(i!=length(chunk_starts)){
      end <- chunk_starts[i+1]-1
    }
    if(i==length(chunk_starts)){
      end <- max_position
    }
      out_i <- cbind(chunk_num, chr, chr_char, start, end)
      out <- rbind(out, out_i)
  }
  out <- as_tibble(out) %>% mutate_all(as.character)

}


launch_bgenie <- function(chr, phenofile, UKBB_dir, chr_char, start_pos, end_pos, chunk_num){

  cat(paste0("Running chr: ", chr, ".\n"))
  #cancel_if(file.exists(paste0("analysis/GWAS/UKBB/chr", chr, "_chunk", chunk_num, ".out")))
  system(paste0("/data/sgg3/jonathan/bgenie_v1.3/bgenie_v1.3_static1 ",
                "--bgen ", UKBB_dir, "imp/_001_ukb_imp_chr", chr, "_v2.bgen ",
                "--pheno ", phenofile, " ",
                "--range ", chr_char, " ", start_pos, " ", end_pos, " ",
                "--pvals --out analysis/GWAS/UKBB/chr", chr, "_chunk", chunk_num, ".out"))
  if(!file.exists(paste0("analysis/GWAS/UKBB/chr", chr, "_chunk", chunk_num, ".out.gz"))){
    file.create(paste0("analysis/GWAS/UKBB/chr", chr, "_chunk", chunk_num, ".out.gz"))
  }
  # Eleonora's command line
  # /data/sgg3/jonathan/bgenie_v1.3/bgenie_v1.3_static1 --bgen
  # /data/sgg3/eleonora/projects/UKBB_GWAS/UK10K_SNPrs/CHR11/chr11.bgen ## this script would not include SNPs specific to HRC_list
  # --pheno ../phenofile --pvals --out chr11.out
}

unzip_bgenie <- function(chr, chunk_num){
  # system(paste0("gzip -dk analysis/GWAS/UKBB/chr", chr, ".out.gz")) ## keep flag doesn't exist on this system

  file <- paste0("analysis/GWAS/UKBB/chr", chr, "_chunk", chunk_num, ".out")
  #cancel_if(file.exists(file))
  if(file.info(paste0(file, ".gz"))$size!=0){
    system(paste0("gunzip < ", paste0(file, ".gz"), " > ", file))} else file.create(file)

}

process_bgenie <- function(directory, extension=".out", HRC_panel){

  current_dir <- getwd()
  setwd(directory)
  system(paste0("tail -q -n +2 *", extension, " > temp"))

  header_file <- fread(list.files(pattern=paste0("\\", extension, "$"))[1], nrows=2, data.table=F)
  bgenie_columns <- colnames(header_file)
  full_data <- fread("temp", data.table=F, header=F)
  setwd(current_dir)
  colnames(full_data) <- bgenie_columns

  HRC_list <- fread(HRC_panel, data.table=F)
  HRC_filtered_data <- full_data %>% filter(rsid %in%  HRC_list$ID)
  return(HRC_filtered_data)

}

read_bgenie_output <- function(bgen_file){

  out <- fread(bgen_file, data.table=F)
  return(out)

}

corr_bmi_slope <- function(bmi_data, ukbb_org, ukb_key, bmi_var){

  data <- bmi_data[[1]]
  cols_sex <- dplyr::pull(ukb_key[which(ukb_key[,1] == sex_var),"col.name"])
  cols_age <- dplyr::pull(ukb_key[which(ukb_key[,1] == age_var),"col.name"])

  cols_pheno <-  dplyr::pull(ukb_key[which(ukb_key[,1] == bmi_var),"col.name"])
  cols_pheno_first <- cols_pheno[1]
  ukbb_sub <- ukbb_org %>% dplyr::select(eid, !!cols_pheno_first, !!cols_sex, !!cols_age) %>%  mutate(outcome = ivt(body_mass_index_bmi_f21001_0_0))
  full_data <- inner_join(data, ukbb_sub)
}

extract_file_info <- function(file, eth){
  drug_class <-  gsub(paste0(".*PSYMETAB_GWAS_(.+)_", eth ,".*"), "\\1", file)
  pheno_temp <- gsub(".*Drug[.](.+)[.]GWAS.txt*", "\\1", file)
  pheno <- gsub(paste0("_", drug_class), "", pheno_temp)
  variable <- unlist(regmatches(pheno, regexpr("[_]", pheno), invert = TRUE))[1]
  outcome <- unlist(regmatches(pheno, regexpr("[_]", pheno), invert = TRUE))[2]

  return(c(drug_class=drug_class, variable=variable, outcome=outcome))
  #return(c(outcome))
}

extract_sig_results <- function(file, threshold, eth){

  result <- numeric()

  data <- fread(file)
  af_file_temp <- str_replace(file, "PSYMETAB_GWAS", "PSYMETAB_GWAS_FREQ")
  af_file_temp2 <- str_replace(af_file_temp, ".GWAS.txt", ".glm.linear")
  af_file <- str_replace(af_file_temp2, "/processed/", paste0("/", eth, "/"))

  file_info <- extract_file_info(file, eth)

  af_data <- fread(af_file)
  data_gwsig <- filter(data, P< threshold)
  if(dim(data_gwsig)[1]!=0){

    data_gwsig$file <- file
    data_gwsig$drug_class <- file_info["drug_class"]
    data_gwsig$variable <- file_info["variable"]
    data_gwsig$outcome <- file_info["outcome"]

    result <- rbind(result, data_gwsig)
  }

  AF_merge <- left_join(result, af_data %>% dplyr::select(ID, A1, A1_CT, A1_FREQ, OBS_CT), by = c("SNP" = "ID"))
  return(AF_merge)

}

merge_psy_UKBB <- function(psy_GWAS, UKBB_GWAS_file){

  ukbb_gwas <- fread(UKBB_GWAS_file)
  ukbb_gwas_sub <- ukbb_gwas %>% dplyr::select(chr, rsid, pos, a_0, a_1, af, info) %>% rename_all(paste0, "_UKBB")
  output <- left_join(psy_GWAS, ukbb_gwas_sub, by = c("SNP" = "rsid_UKBB"))
  return(output)

}

add_AF <- function(data, AF_file){

  AF_data <- fread(AF_file, data.table=F)
  AF_merge <- left_join(data, AF_data %>% dplyr::select(ID, ALT_FREQS), by = c("SNP" = "ID"))


}

## PLINK beta in terms of A1
## Bgenie beta in terms of a_1

check_palindromic <- function(A1, A2){
	(A1 == "T" & A2 == "A") |
	(A1 == "A" & A2 == "T") |
	(A1 == "G" & A2 == "C") |
	(A1 == "C" & A2 == "G")
}

flip_alleles <- function(x){
	x <- toupper(x)
	x <- gsub("C", "g", x)
	x <- gsub("G", "c", x)
	x <- gsub("A", "t", x)
	x <- gsub("T", "a", x)
	return(toupper(x))
}


match_alleles <- function(A1, A2, B1, B2){

  match_description <- rep("simple_match", length(A1))
  status1 <- (A1 == B1) & (A2 == B2)
  to_swap <- (A1 == B2) & (A2 == B1)

  # If B's alleles are the wrong way round then swap
	Btemp <- B1[to_swap]
	B1[to_swap] <- B2[to_swap]
	B2[to_swap] <- Btemp
  match_description[to_swap] <- "swap_alleles"

	# Check again
	status1 <- (A1 == B1) & (A2 == B2)
	palindromic <- check_palindromic(A1, A2)

	# If NOT palindromic and alleles DON'T match then try flipping
	i <- !palindromic & !status1
	B1[i] <- flip_alleles(B1[i])
	B2[i] <- flip_alleles(B2[i])
	status1 <- (A1 == B1) & (A2 == B2)
	#jlog[['flipped_alleles_basic']] <- sum(i)

	# If still NOT palindromic and alleles DON'T match then try swapping
	i <- !palindromic & !status1
	to_swap <- (A1 == B2) & (A2 == B1)
  match_description[to_swap] <- "swap_alleles"
	Btemp <- B1[to_swap]
	B1[to_swap] <- B2[to_swap]
	B2[to_swap] <- Btemp

	# Any SNPs left with unmatching alleles need to be removed
	status1 <- (A1 == B1) & (A2 == B2)
	remove <- !status1
  match_description[remove] <- "non-match"
  return(list(match_description=match_description, palindromic=palindromic))

  ## since all datasets are imputed and on positive strand, we aren't really worried about palindromic SNPs
}

get_plink_a2 <- function(plink_glm){

  a2 <- plink_glm$ALT
  ref_a1_mis_match <- !plink_glm$REF == plink_glm$A1
  a2[ref_a1_mis_match] <- plink_glm$REF[ref_a1_mis_match]
  return(a2)

}

harmonize_ukbb_plink_data <- function(data, SNP_col="SNP", REF1_col="A1", ALT1_col="A2", REF2_col="a_1_UKBB", ALT2_col="a_0_UKBB"){

  data$A2 <- get_plink_a2(data)
  data <- data %>% filter(!is.na(bmi_slope_resid_beta_UKBB))

  ## PLINK beta in terms of A1
  ## Bgenie beta in terms of a_1

  temp <- match_alleles(data$A1, data$A2, data$a_1_UKBB, data$a_0_UKBB)
  data$palindromic <- temp$palindromic
  data$match_description <- temp$match_description

  return(data)
}

flip_beta <- function(betas_to_flip, match_description){
  new_beta <- betas_to_flip
  i <- match_description=="swap_alleles"
  new_beta[i] <- new_beta[i]*(-1)
  return(new_beta)
}

calc_het <- function(data, snp_col, beta1_col, beta2_col, se1_col, se2_col, match_description_col){

  data$beta2_harmonized <- flip_beta(data[[beta2_col]], data[[match_description_col]])

  data$het_pval <- NA
  for(snp in 1:dim(data)[1] ){
    rsid <- as.character(data[snp,snp_col])
    beta1 <-  as.numeric(as.character(data[snp,beta1_col]))
    beta2 <-   as.numeric(as.character(data[snp,"beta2_harmonized"]))
    SE1 <-  as.numeric(as.character(data[snp,se1_col]))
    SE2 <-  as.numeric(as.character(data[snp,se2_col]))
    se <- sqrt( (SE1^2) + (SE2^2) )
    t <- (beta1-beta2)/se
    p_het <- 2*pnorm(-abs(t))
    #temp <- cbind(rsid, p_het)
    data$het_pval[snp] <- p_het
    #het_out <- rbind(het_out,temp)

  }

  return(data)

}

prune_psy <-function(data){

  # unique_files <- unique(data$file)
  #
  result <- numeric()
  # for(i in 1:length(unique_files)){
  #
  # file_i <- unique_files[i]
  # t <- subset(data, data$file==file_i)
  t <- data[order(data$P),] %>% filter(!is.na(chr_UKBB))
  j <- 1

  while(j <= dim(t)[1]){
    chr <- t[["CHR"]][j]
    bp <- t[["BP"]][j]
    remove_indices <- which(t[["CHR"]]==chr & t[["BP"]] >= (bp - 500000) & t[["BP"]] <= (bp + 500000))
    number_signals <- length(remove_indices)
    remove_indices <- remove_indices[-which(remove_indices==j)]
    if(length(remove_indices)!=0){
      t <- t[-remove_indices, ]
    }
    t$number_signals <- number_signals
    j <- j+1

  }

  #result <- rbind(result, t)

  # }

  return(t)

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
    z_stat <- qnorm(.025,lower.tail=FALSE)
    L95 <- beta-z_stat*se
    U95 <- beta+z_stat*se
    temp <- cbind(beta, se, pval, L95, U95)
    colnames(temp) <- paste0(coef, c("_beta", "_se","_pval", "_L95", "_U95"))
    summary <- cbind (summary,  temp)
  }
  return(summary)

}


analyze_prs <- function(prs_file, pheno_file, pc_file, linear_pheno_columns, glm_pheno_columns, covars){

  out <- list()
  prs_data <- fread(prs_file, data.table=F)
  pheno_data <- fread(pheno_file, data.table=F)
  pc_eth_data <- fread(pc_file, data.table=F)
  colnames(prs_data) <- c("ID", "GEN", "PRS_5e08", "PRS_5e07", "PRS_5e06", "PRS_5e05", "PRS_5e04", "PRS_5e03", "PRS_5e02", "PRS_0.1")

  output <- numeric()

  pc_data <- pc_eth_data[,c("GPCR", "eth", paste0("PC", 1:20))] %>% rename("GEN" = "GPCR")
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
          mutate(prs = scale(prs)) %>%
          mutate(outcome = ivt(outcome))

        model <- (lm(outcome ~ ., data=joint))
        model_summary <- full_model_summary(model)
        model_stats <- extract_model_stats(model_summary, c("prs"))
        n_model <- nobs(model)
        row <- cbind(prs, eth, var, model_type, model_stats, n_model)
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


        model_summary <- full_model_summary(model)
        model_stats <- extract_model_stats(model_summary, c("prs"))
        n_model <- nobs(model)

        row <- cbind(prs, eth, var, model_type, model_stats, n_model)
        output <- rbind(output, row)
      }

    }
  }

  out[[prs_file]] <- as_tibble(output)

  return(out)

}



write_ukbb_prs_analysis <- function(data, out_file_name){

  filter_data <- data %>% filter(eth=="CEU") %>% filter(prs=="PRS_0.1")
  write.csv(filter_data, out_file_name, row.names=F)
  return(out_file_name)
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

define_baseline_files <- function(info = baseline_gwas_process, output = "PSYMETAB_GWAS", eths,
  output_dir = "analysis/GWAS", type = "full"){
  file_name <- output
  output_suffix <- info$output_suffix
  pheno <- info$pheno
  pheno_names <- numeric()
  eth_keep <- numeric()
  eth_list <- numeric()
  write_keep <- numeric()
  log_keep <- numeric()
  for(eth in c(eths, "META")){
    eth_dir <- file.path(output_dir, type)
    if(eth!="META"){
      file_name_eth <- paste0(file_name,"_",eth)
      eth_output <- file.path(eth_dir,eth,paste0(file_name_eth, ".", pheno, ".glm.linear.gz"))
      log_file <- file.path(eth_dir,eth,paste0(file_name_eth, ".", "log"))
      processed_file <- gsub(".glm.linear.gz", ".GWAS.txt", eth_output)

    }
    if(eth=="META"){
      file_name_eth <- file_name
      eth_output <- file.path(eth_dir,eth,paste0(file_name_eth, ".", pheno, ".meta"))
      log_file <- file.path(eth_dir,eth,paste0(file_name_eth, ".", pheno, ".log"))
      processed_file <- gsub(".meta", ".GWAS.txt", eth_output)
    }

    eth_list <- c(eth_list, rep(eth, length(which(file.exists(eth_output)))))
    eth_keep <- c(eth_keep, eth_output[which(file.exists(eth_output))])
    pheno_names <- c(pheno_names, pheno[which(file.exists(eth_output))])
    if(eth!="META"){
      log_keep <- c(log_keep, rep(log_file, length(which(file.exists(eth_output)))))
    }
    if(eth=="META"){
      log_keep <- c(log_keep, log_file[which(file.exists(eth_output))])
    }

    write_keep <- c(write_keep, processed_file[which(file.exists(eth_output))])
    write_keep <- gsub(paste0("/", eth, "/"), "/processed/", write_keep)
  }

  out <- tibble(eth = eth_list,
                pheno = pheno_names,
                drug = NA,
                log_file = log_keep,
                plink_file = eth_keep,
                processed_file = write_keep)
  return(out)
}

define_interaction_files <- function(info = interaction_gwas_process, output = "PSYMETAB_GWAS", eths,
  output_dir = "analysis/GWAS", type = "interaction", pheno_list){

  file_name <- output
  output_suffix <- info$output_suffix
  #file_name <- paste0(file_name,"_",output_suffix)
  #pheno <- pheno_list
  pheno_names <- numeric()
  eth_keep <- numeric()
  eth_list <- numeric()
  write_keep <- numeric()
  drug_class <- numeric()
  log_keep <- numeric()
  for(eth in c(eths, "META")){
    eth_dir <- file.path(output_dir, type)
    for(pheno in pheno_list){

      if(eth!="META"){
        file_name_eth <- paste0(file_name,"_", output_suffix, "_int_",eth)
        #eth_output <- file.path(eth_dir,eth,paste0(file_name_eth, ".", pheno, ".glm.linear"))
        eth_combos <- apply(expand.grid(pheno, output_suffix), 1, paste, collapse="_")
        eth_output <- file.path(eth_dir,eth,paste0(file_name_eth, ".", eth_combos,  ".glm.linear.interaction.gz"))
        log_file <- file.path(eth_dir,eth,paste0(file_name_eth, ".", "log"))
        processed_file <- gsub(".glm.linear.interaction.gz", ".GWAS.txt", eth_output)
      }
      if(eth=="META"){
        file_name_eth <- paste0(file_name,"_", output_suffix, "_int")
        eth_output <- file.path(eth_dir,eth,paste0(file_name_eth, ".", pheno, ".meta"))
        processed_file <- gsub(".meta", ".GWAS.txt", eth_output)
        log_file <- file.path(eth_dir,eth,paste0(file_name_eth, ".", pheno, ".log"))
      }
      eth_list <- c(eth_list, rep(eth, length(which(file.exists(eth_output)))))
      eth_keep <- c(eth_keep, eth_output[which(file.exists(eth_output))])
      write_keep <- c(write_keep, processed_file[which(file.exists(eth_output))])
      write_keep <- gsub(paste0("/", eth, "/"), "/processed/", write_keep)
      pheno_names <- c(pheno_names, rep(pheno, length(which(file.exists(eth_output)))))
      drug_class <- c(drug_class, output_suffix[which(file.exists(eth_output))])
      log_keep <- c(log_keep, log_file[which(file.exists(eth_output))])

    }


  }

  out <- tibble(eth = eth_list,
                pheno = pheno_names,
                drug_class = drug_class,
                drug = NA,
                log_file = log_keep,
                plink_file = eth_keep,
                processed_file = write_keep)

  return(out)

}

define_subgroup_files <- function(info = subgroup_gwas_process, output = "PSYMETAB_GWAS", eths,
  output_dir = "analysis/GWAS", type = "subgroup", pheno_list){

    file_name <- output
    output_suffix <- info$output_suffix
    pheno_names <- numeric()
    eth_keep <- numeric()
    eth_list <- numeric()
    drug_keep <- numeric()
    write_keep <- numeric()
    drug_class <- numeric()
    log_keep <- numeric()

    for(eth in c(eths, "META")){
      for(pheno in pheno_list){
        eth_dir <- file.path(output_dir, type)

        if(eth!="META"){
          file_name_eth <- paste0(file_name,"_", output_suffix, "_",eth)
          eth_combos <- apply(expand.grid(pheno, output_suffix), 1, paste, collapse="_")

          file_name_eth_drug <- paste0(file_name_eth, ".Drug")
          file_name_eth_nodrug <- paste0(file_name_eth, ".NoDrug")
          eth_output_drug <- file.path(eth_dir,eth,paste0(file_name_eth_drug, ".", eth_combos, ".glm.linear.gz"))
          eth_output_nodrug <- file.path(eth_dir,eth,paste0(file_name_eth_nodrug, ".", eth_combos, ".glm.linear.gz"))
          processed_file_drug <- gsub(".glm.linear.gz", ".GWAS.txt", eth_output_drug)
          processed_file_nodrug <- gsub(".glm.linear.gz", ".GWAS.txt", eth_output_nodrug)

          log_file <- file.path(eth_dir,eth,paste0(file_name_eth, ".", "log"))
        }
        if(eth=="META"){
          file_name_eth <- paste0(file_name,"_", output_suffix)
          file_name_eth_drug <- paste0(file_name_eth, ".Drug")
          file_name_eth_nodrug <- paste0(file_name_eth, ".NoDrug")
          eth_output_drug <- file.path(eth_dir,eth,paste0(file_name_eth_drug, ".", pheno, ".meta"))
          eth_output_nodrug <- file.path(eth_dir,eth,paste0(file_name_eth_nodrug, ".", pheno, ".meta"))
          processed_file_drug <- gsub(".meta", ".GWAS.txt", eth_output_drug)
          processed_file_nodrug <- gsub(".meta", ".GWAS.txt", eth_output_nodrug)

          log_file_drug <- file.path(eth_dir,eth,paste0(file_name_eth_drug, ".", pheno, ".log"))
          log_file_nodrug <- file.path(eth_dir,eth,paste0(file_name_eth_nodrug, ".", pheno, ".log"))
          log_file <- c(log_file_drug, log_file_nodrug)
        }

        eth_output <- c(eth_output_drug, eth_output_nodrug)
        drug_list <- c(rep("Drug", length(output_suffix)), rep("NoDrug", length(output_suffix)))

        drug_keep <- c(drug_keep, drug_list[which(file.exists(eth_output))])
        eth_list <- c(eth_list, rep(eth, length(which(file.exists(eth_output)))))
        eth_keep <- c(eth_keep, eth_output[which(file.exists(eth_output))])
        pheno_names <- c(pheno_names, rep(pheno, length(which(file.exists(eth_output)))))
        drug_class_list <- c(output_suffix, output_suffix)
        drug_class <- c(drug_class, drug_class_list[which(file.exists(eth_output))])

        if(eth!="META"){
          log_keep <- c(log_keep, c(log_file, log_file)[which(file.exists(eth_output))])
        }
        if(eth=="META"){
          log_keep <- c(log_keep, log_file[which(file.exists(eth_output))])
        }


        write_keep <- c(write_keep, c(processed_file_drug, processed_file_nodrug)[which(file.exists(eth_output))])
        write_keep <- gsub(paste0("/", eth, "/"), "/processed/", write_keep)
      }

  }

  out <- tibble(eth = eth_list,
                pheno = pheno_names,
                drug_class = drug_class,
                drug = drug_keep,
                log_file = log_keep,
                plink_file = eth_keep,
                processed_file = write_keep)

  return(out)

}


define_case_only_files <- function(info = case_only_gwas_process, output = "PSYMETAB_GWAS", eths,
  output_dir = "analysis/GWAS", type = "case_only", pheno_list = interaction_outcome){
  file_name <- output
  output_suffix <- info$output_suffix
  pheno <- info$pheno
  pheno_names <- numeric()
  eth_keep <- numeric()
  eth_list <- numeric()
  write_keep <- numeric()
  log_keep <- numeric()
  for(eth in c(eths,  "META_DRUGS")){
    eth_dir <- file.path(output_dir, type)
    if(eth!="META" | eth!="META_DRUGS"){
      file_name_eth <- paste0(file_name,"_",eth)
      eth_output <- file.path(eth_dir,eth,paste0(file_name_eth, ".", pheno, ".glm.linear.gz"))
      log_file <- file.path(eth_dir,eth,paste0(file_name_eth, ".", "log"))
      processed_file <- gsub(".glm.linear.gz", ".GWAS.txt", eth_output)

    }
    if(eth=="META"){
      file_name_eth <- file_name
      eth_output <- file.path(eth_dir,eth,paste0(file_name_eth, ".", pheno, ".meta"))
      log_file <- file.path(eth_dir,eth,paste0(file_name_eth, ".", pheno, ".log"))
      processed_file <- gsub(".meta", ".GWAS.txt", eth_output)
    }

    if(eth=="META_DRUGS"){
      # pheno_without_drug <- unique(sub("_[^_]+$", "", pheno))
      pheno_eths <- apply(expand.grid(eths, pheno_list), 1, paste, collapse=".")
      file_name_eth <- file_name
      eth_output <- file.path(eth_dir,eth,paste0(file_name_eth, "_", pheno_eths, ".meta"))
        # apply(expand.grid(eths, pheno_without_drug), 1, paste, collapse=".")
      log_file <- file.path(eth_dir,eth,paste0(file_name_eth, "_", pheno_eths, ".log"))
      processed_file <- gsub(".meta", ".GWAS.txt", eth_output)
    }

    eth_list <- c(eth_list, rep(eth, length(which(file.exists(eth_output)))))
    eth_keep <- c(eth_keep, eth_output[which(file.exists(eth_output))])

    if(eth=="META_DRUGS"){
      pheno_names <- c(pheno_names, pheno_eths[which(file.exists(eth_output))])
    } else pheno_names <- c(pheno_names, pheno[which(file.exists(eth_output))])

    if(eth!="META" & eth!="META_DRUGS"){
      log_keep <- c(log_keep, rep(log_file, length(which(file.exists(eth_output)))))
    }
    if(eth=="META" | eth=="META_DRUGS"){
      log_keep <- c(log_keep, log_file[which(file.exists(eth_output))])
    }

    write_keep <- c(write_keep, processed_file[which(file.exists(eth_output))])
    write_keep <- gsub(paste0("/", eth, "/"), "/processed/", write_keep)
  }

  out <- tibble(eth = eth_list,
                pheno = pheno_names,
                drug = NA,
                log_file = log_keep,
                plink_file = eth_keep,
                processed_file = write_keep)
  return(out)
}

process_gwas <- function(eth, pheno, drug, file, output = "PSYMETAB_GWAS", output_dir = "analysis/GWAS", type = "full",
  info_file = "analysis/QC/15_final_processing/PSYMETAB_GWAS.info", out_file){

  #outcome_variable,interaction_variable,model,
  info <- fread(info_file)
  write_dir <- file.path(output_dir, type, "processed")

  gwas_result <- fread(file, data.table=F, stringsAsFactors=F)
  if(dim(gwas_result)[1]!=0){
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
  }

  if(dim(gwas_result)[1]==0){

    joint <- as.data.frame(rbind(c(1, 1, "rs12345", "A", "G", 1, 1, 1, 0.001)))
    joint <- joint[rep(seq_len(nrow(joint)), each = 10), ]
    colnames(joint) <- c("CHR", "BP", "SNP", "REF", "ALT", "BETA", "SE", "T_STAT", "P")
    fwrite(joint, out_file, sep="\t")
  }

  #return(list(manhattan_file_name=manhattan_file_name, qq_file_name=qq_file_name, title=title))
}

gwas_figures_input <- function(eth, pheno, drug, file, output = "PSYMETAB_GWAS", output_dir = "analysis/GWAS", type = "full",
  info_file = "analysis/QC/15_final_processing/PSYMETAB_GWAS.info", out_file, drug_class = NULL){

  write_dir <- file.path(output_dir, type, "processed")

  if(type == "interaction" | type == "subgroup"){
    pheno <- paste0(pheno, "_", drug_class)
  }
  title <- ifelse(is.na(drug), paste0(pheno, "_", eth), paste0(pheno, "_", drug, "_", eth))

  manhattan_file_name <- paste0(write_dir, "/", output, "_", title, "_manhattan.png")
  qq_file_name <- paste0(write_dir, "/", output, "_", title, "_qq.png")
  return(list(joint_file=out_file, manhattan_file_name=manhattan_file_name, qq_file_name=qq_file_name, title=title))

}

create_figures <- function(joint_file, manhattan_file_name, qq_file_name, title, gw_sig){

  joint <- fread(joint_file, data.table=F)

  gw_sig_filter <- joint[which(joint$P < gw_sig),]
  #if(dim(gw_sig_filter)[1]!=0){
    joint_sig_filter <- joint[which(joint$P < 0.01),]
    png(manhattan_file_name, width=2000, height=1000, pointsize=18)
    manhattan(joint_sig_filter, main=title)
    dev.off()
  #}

  rm(joint_sig_filter)

  joint_rand_filter <- joint[sample(nrow(joint), ceiling(0.1*dim(joint)[1])), ]

  rm(joint)
  png(qq_file_name, width=2000, height=1000, pointsize=18)
  qq(joint_rand_filter$P, main=title)
  dev.off()

}

gw_sig_extract <- function(gwas_file, gw_sig_nominal, eth, pheno, drug_class, drug){

  output <- list()
  gwas <- fread(gwas_file, data.table=F)
  gwas$eth <- eth
  gwas$pheno <- pheno
  gwas$drug_class <- drug_class
  gwas$drug <- drug

  gw_sig_filter <- gwas[which(gwas$P < gw_sig_nominal),]


  output_name <- paste(eth, pheno, drug_class, drug, sep=".")
  if(is.na(drug)){
    output_name <- paste(eth, pheno, drug_class, sep=".")
  }
  output[[output_name]] <- gw_sig_filter
  return(output)

}

clean_gwas_summary <- function(gwas_data){

  output_name <- names(gwas_data)
  gwas_data <- gwas_data[[1]]

  if(dim(gwas_data)[1]==0){
    gwas_data$variable <- character()
    gwas_data$outcome <- character()
  }
  if(dim(gwas_data)[1]!=0){

    pheno_name <- gwas_data$pheno[1]

    variable <- str_extract(string=pheno_name, pattern="[^_]*")
    outcome <- str_replace(pheno_name, paste0(variable, "_"), "")

    gwas_data$variable <- variable
    gwas_data$outcome <- outcome

  }

  output <- list()
  output[[output_name]] <- gwas_data
  return(output)
}

clean_case_only_summary <- function(case_only_data){

  output_name <- names(case_only_data)
  case_only_data <- case_only_data[[1]]

  if(dim(case_only_data)[1]==0){

    case_only_data$group <- character()
    case_only_data$variable <- character()
    case_only_data$outcome <- character()
  }
  if(dim(case_only_data)[1]!=0){
    pheno_name <- case_only_data$pheno[1]
    group <- case_only_data$eth[1]

    if(group == "META_DRUGS"){
      eth <- str_extract(string=pheno_name, pattern="[^\\.]*")
      drug <- NA
      pheno <-  str_extract(string=pheno_name, pattern="[^\\.]+$")

    } else {
      eth <- group
      drug <-  drug <- sub('.*\\_', '', pheno_name)
      pheno <- sub("_[^_]+$", "", pheno_name)
    }

    outcome <- str_extract(string=pheno, pattern="[^_]+$")


    variable <- str_extract(string=pheno, pattern="[^_]*")
    outcome <- str_replace(pheno, paste0(variable, "_"), "")

    case_only_data$drug <- drug
    case_only_data$eth <- eth
    case_only_data$group <- group
    case_only_data$variable <- variable
    case_only_data$outcome <- outcome
  }

  output <- list()
  output[[output_name]] <- case_only_data
  return(output)


}

summarize_case_only_meta <- function(snp, eth, variable, outcome, study_name, drug_prioritization, directory){

  pheno <- paste0(eth, ".", variable, "_", outcome)
  meta_file <- paste0(study_name, "_", pheno, ".meta")
  meta_gwas <- fread(file.path(directory, "META_DRUGS", meta_file), data.table=F)

  meta_row <- which(meta_gwas$SNP==snp)
  P_META <- as.numeric(meta_gwas[meta_row, "P"])
  BETA_META <- as.numeric(meta_gwas[meta_row, "BETA"])
  SE_META <-abs(BETA_META/qnorm(P_META/2))

  drug_result <- numeric()
  for(drug in drug_prioritization){
    drug_file <- paste0(study_name, "_", pheno, "_", drug, ".GWAS.txt")
    if(file.exists(file.path(directory, "processed", drug_file))){
      drug_gwas <- fread(file.path(directory, "processed", drug_file), data.table=F)
      row_num <-  which(drug_gwas$SNP==snp)
      if(length(row_num)!=0){
        drug_p <- as.numeric(drug_gwas[row_num, "P"])
        drug_beta <- as.numeric(drug_gwas[row_num, "BETA"])
        drug_se <- as.numeric(drug_gwas[row_num, "SE"])
      } else {
        drug_p <- NA
        drug_beta <- NA
        drug_se <- NA
      }
    }
    if(!file.exists(file.path(directory, "processed", drug_file))){
      drug_p <- NA
      drug_beta <- NA
      drug_se <- NA
    }
    drug_result <- cbind(drug_result, drug_p, drug_beta, drug_se)

  }

  colnames(drug_result) <- apply(expand.grid( c("P_", "BETA_", "SE_"), drug_prioritization),
                               1, paste, collapse="")


  result <- cbind(snp, eth, variable, outcome, drug_result, P_META, BETA_META, SE_META)

}

count_GWAS_n <- function(psam_file, pheno_file, output_suffix, subgroup=NA, covar_file=NA, covars=NA, eths,
  eth_sample_file, related_ids_file){

  sample_file <- read.table(psam_file, header=F)
  nsamples <- dim(sample_file)[1]

  pheno_data <- fread(pheno_file, data.table=F)
  remove_samples <- read.table(related_ids_file, header=F)

  if(!is.na(covar_file)){
    covar_list <- unlist(covars)
    covar_data <- fread(covar_file, data.table=F)
    pheno_data <- left_join(pheno_data, covar_data %>% dplyr::select(FID, IID, !!covar_list))
  }

  output <- numeric()

  for (eth in eths){
    keep_file <- str_replace(eth_sample_file, "ETH", eth)
    eth_samples <- read.table(keep_file, header=F)

    final_sample_list <- sample_file %>%
      filter(V1 %in% eth_samples$V1) %>%
      filter(!(V1 %in% remove_samples$V1))

    eth_count <- dim(final_sample_list)[1]

    if(eth_count > 100){
      pheno_data_eth <- pheno_data %>% filter(FID %in% final_sample_list$V1)
      colnames_pheno_data_eth <- colnames(pheno_data_eth)

      # if subgroup
      if(!is.na(subgroup)){
        vars_analyzed <- colnames_pheno_data_eth[-which(colnames_pheno_data_eth %in% c("FID", "IID", subgroup))]
        for(var in vars_analyzed){
          for(drug_class in c("Drug", "NoDrug")){

            pheno_i <- pheno_data_eth[which(pheno_data_eth[[subgroup]]==drug_class),c("FID", var)]
            pheno_i_complete <- pheno_i[complete.cases(pheno_i),]
            n_samples <- dim(pheno_i_complete)[1]

            row <- cbind(pheno_file, eth, var, output_suffix, drug_class, n_samples)
            output <- rbind(output, row)
          }
        }
      }

      # if interaction
      if(!is.na(covar_file)){
        vars_analyzed <- colnames_pheno_data_eth[-which(colnames_pheno_data_eth %in% c("FID", "IID"))]
        for(var in vars_analyzed){
          drug_class <- NA
          pheno_i <- pheno_data_eth[,c("FID", var, covar_list)]
          pheno_i_complete <- pheno_i[complete.cases(pheno_i),]
          n_samples <- dim(pheno_i_complete)[1]
          row <- cbind(pheno_file, eth, var, output_suffix, drug_class, n_samples)
          output <- rbind(output, row)

        }
      }

      # if baseline
      if(output_suffix=="baseline"){
        vars_analyzed <- colnames_pheno_data_eth[-which(colnames_pheno_data_eth %in% c("FID", "IID"))]
        for(var in vars_analyzed){
          drug_class <- NA
          pheno_i <- pheno_data_eth[,c("FID", var)]
          pheno_i_complete <- pheno_i[complete.cases(pheno_i),]
          n_samples <- dim(pheno_i_complete)[1]
          row <- cbind(pheno_file, eth, var, output_suffix, drug_class, n_samples)
          output <- rbind(output, row)
        }

      }

    }
  }
  colnames(output) <- c("pheno_file", "eth", "variable_name", "suffix", "drug_class", "n_sample")
  output <- as_tibble(output)
  return(output)
}


extract_bgen <- function(snps, file){

  num_snps <- length(snps)
  bgen_data <- rbgen::bgen.load(file=file, rsids = snps)
  snp_map <- bgen_data$variants
  snp_map$AF <- NA
  geno <- numeric()

  for(i in 1:dim(snp_map)[1])
  {
      snp <- as.character(snp_map[["rsid"]][i])
      snp_temp <- bgen_data$data[snp,,]
      snp_temp <- as.data.frame(snp_temp)
      snp_temp[,4] <- snp_temp[,2]+2*snp_temp[,3]
      snp_map$AF[i] <- mean(snp_temp[,4], rm.na=T)/2
      snp_out <- snp_temp[,4]
      geno <- cbind(geno, snp_out)
  }
  colnames(geno) <- bgen_data$variants$rsid
  return(list(geno_data=geno,snp_map=snp_map))
}

load_geno <- function(bgen_sample_file,snp_data, UKBB_dir){

  ####### Extract SNP froms bgen files
  snp_map <- numeric()
  geno_data <- numeric()
  for(chr in 1:22)
  {
    snps <- pull(snp_data[which(snp_data[,"chr"]==chr),], rsid)
    if(length(snps)==0) next
    bgen_file=paste(UKBB_dir, "/imp/", "_001_ukb_imp_chr", chr, "_v2.bgen", sep="")
    bgen_chr_extract <- extract_bgen(unique(snps),bgen_file)
    geno_data <- cbind(geno_data, bgen_chr_extract$geno_data)
    snp_map <- rbind(snp_map,bgen_chr_extract$snp_map)
    cat(paste0("Finished loading chr: ", chr, ".\n"))

  }
  row.names(geno_data) <- bgen_sample_file[-1,1]
  return(list(geno_data, snp_map))
}

UKBB_GWAS <- function(ukbb_geno, phenotype_file, phenotype_col, phenotype_description,IV_file,sample_file,pheno_cov,grouping_var,household_time,output){

  geno_data <- ukbb_geno[[1]]
  snp_map <- ukbb_geno[[2]]

  ####### load IV list
  IV_data <- fread(IV_file, header=T, data.table=F)

  ####### load phenotypes: GRS and raw
  phenotype <- fread( phenotype_file,header=T, data.table=F, select=c("IID",phenotype_col))
  # grs_pheno <- fread( paste0(pheno_dir,"/GRS/",outcome_sex, "/",trait_ID,"_",outcome_sex,".best"),header=T, data.table=F)

  ## extract best threshold found by PRSice
  # grs_summary <- read.table(paste0(pheno_dir,"/GRS/",outcome_sex, "/",trait_ID,"_",outcome_sex,".summary"), header=T)
  # grs_best <- grs_summary[["Threshold"]][1]
  # joint_pheno <- merge(raw_pheno, grs_pheno, by="IID")
  # cat(paste0("The best threshold chosen by PRSice was: ", grs_best, "\n"))

  outcome_GWAS <- numeric()

  cat(paste0("Running regression for each GWAS in all participants and bins of household pairs
  (divided by length of time in household) using specified phenotype as outcomes...\n"))

  for(k in 1:dim(snp_map)[1])
  {
    snp <- IV_data[["rsid"]][k]
    geno_sub <- as.data.frame(IV_geno[,which(colnames(IV_geno)==snp)])
    colnames(geno_sub)[1] <- "geno"
    geno_sub$IID <- as.numeric(row.names(geno_sub))
    temp1 <- merge(pheno_cov,geno_sub, by.x=index, by.y="IID")
    temp2 <- merge(temp1, phenotype, by.x=opp_index, by.y="IID")
    temp3 <- merge(temp2, household_time[,c("HOUSEHOLD_MEMBER1",grouping_var)], by="HOUSEHOLD_MEMBER1")
    # final_data <- merge(temp3, grs_pheno[,c("IID","PRS")], by.x=opp_index, by.y="IID")
    final_data <- temp3
    # colnames(final_data) <- c(colnames(pheno_cov), "geno_index", "pheno_household", "time_together_bin", "PRS_household"  )###,"GRS_0.01_household","GRS_0.001_household")
    colnames(final_data) <- c(colnames(pheno_cov), "geno_index", "phenotype", grouping_var)###,"GRS_0.01_household","GRS_0.001_household")

    outcome <- "phenotype"
    pheno_run <- final_data[,c(grep("_age", names(final_data)),grep("_PC_", names(final_data)), which(names(final_data)=="geno_index" |names(final_data)==outcome))]
    colnames(pheno_run)[which(colnames(pheno_run)==outcome)] <- "outcome"
    pheno_run$outcome <- scale(pheno_run$outcome)
    mod <- glm(outcome ~ ., data=pheno_run, family="gaussian") # check that "grouping_var" is removed

    model_summary <- full_model_summary(mod)

    mod_extract <- extract_model_stats(model_summary, c("geno_index"))
    n <- length(fitted(mod))

    k_GWAS <- numeric()
    k_GWAS <- rbind(k_GWAS,cbind(phenotype_file, outcome, as.character(snp), "all", mod_extract,n))

    household_intervals <- levels(household_time[[grouping_var]])
    for(bin in household_intervals)
    {
      bin_sub <- final_data[which(final_data[[grouping_var]]==bin),]
      bin_pheno_run <- bin_sub[,c(grep("_age", names(bin_sub)),grep("_PC_", names(bin_sub)), which(names(bin_sub)=="geno_index" |names(bin_sub)==outcome))]
      colnames(bin_pheno_run)[which(colnames(bin_pheno_run)==outcome)] <- "outcome"
      bin_pheno_run$outcome <- scale(bin_pheno_run$outcome)
      bin_pheno_run <- bin_pheno_run[complete.cases(bin_pheno_run),] #if all values are the same then NA's will be produced
      if(dim(bin_pheno_run)[1]!=0){
        bin_mod <- glm(outcome ~ ., data=bin_pheno_run, family="gaussian")

        bin_model_summary <- full_model_summary(bin_mod)

        bin_mod_extract <- extract_model_stats(bin_model_summary, c("geno_index"))
        bin_n <- length(fitted(bin_mod))
        bin_row_summary <- cbind(phenotype_file, outcome, as.character(snp), bin, bin_mod_extract,bin_n)

      } else bin_row_summary <- cbind(phenotype_file, outcome, as.character(snp), bin, NA,NA,NA,0)

      k_GWAS <- rbind(k_GWAS, bin_row_summary)

      #assign(paste0("outcome_GWAS_bin",bin), rbind(get(paste0("outcome_GWAS_bin",bin)), bin_row_summary))

    }

    outcome_GWAS <- rbind(outcome_GWAS, k_GWAS)
    cat(paste0("Finished running regression for ", k, " of ", dim(IV_data)[1], " SNPs.\n"))

  }

  outcome_GWAS <- as.data.frame(outcome_GWAS)
  colnames(outcome_GWAS)[2] <- "phenotype_file"
  colnames(outcome_GWAS)[2] <- "outcome"
  colnames(outcome_GWAS)[3] <- "SNP"
  outcome_GWAS$outcome <- phenotype_col
  colnames(outcome_GWAS)[4] <- grouping_var
  outcome_gwas_out <- merge(outcome_GWAS, snp_map,by.x="SNP", by.y="rsid")

  out <- list(list(outcome_gwas_out = outcome_gwas_out))
  return(out)
  #write.csv(outcome_gwas_out, paste0(pheno_dir, "/household_GWAS/outcome_",outcome_sex,"/", phenotype_description,"_" ,exposure_sex, "-",outcome_sex, "_GWAS.csv"), row.names=F)
  #cat(paste0("Finished processing and writing GWAS files for ", outcome_sex, "s.\n"))

}


sort_ukbb_comparison <- function(ukbb_comparison, subgroup_GWAS_count, interaction_outcome, drug_classes){

  ukbb_comparison$n_sample <- NA
  for(i in 1:dim(ukbb_comparison)[1]){

    file <- basename(as.character(ukbb_comparison$file[i]))

    interaction_temp <- interaction_outcome[which(sapply(interaction_outcome, function(x) {grepl(paste0(x, "_"), file)}))]
    drug <- drug_classes[which(sapply(drug_classes, function(x) {grepl(paste0("_", x, "_"), file)}))]

    interaction_drug_temp <- paste0(interaction_temp, "_", drug)
    interaction_exact <- which(sapply(interaction_drug_temp, function(x) {grepl(paste0(".", x, "."), file)}))

    row <- which(subgroup_GWAS_count$eth=="CEU" & subgroup_GWAS_count$drug_class=="Drug" & subgroup_GWAS_count$variable_name==paste0(interaction_temp[interaction_exact], "_", drug) & subgroup_GWAS_count$suffix==drug)

    count <- subgroup_GWAS_count$n_sample[row]
    ukbb_comparison$n_sample[i] <- count
  }

  return(ukbb_comparison)

}
