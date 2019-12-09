
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

munge_qc_pheno <- function(raw_pheno){

}

read_pcs <- function(pc_dir){

  pclist = list()
  for (eth in eths)
  {
    PC_temp<- fread(pc_dir, eth, paste0(study_name,"_",eth,"_projections.txt")))
    print(eth)
    print(dim(PC_temp))
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
