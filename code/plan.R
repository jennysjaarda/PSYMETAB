
# Workflow Plans
#
# ## Data sources
#
# The raw data are in `.xlsx` files. Some of the data are pre-processed
# (mean and sd of count and length of gemmae, 30 min averages of PPFD).
#
# The files include some additional analsyes and other data that won't
# be used here.

plan <- drake_plan (

  # Load and clean data ----

  # ### Microclimate data
  #
  # Microclimate variables include PPFD (photon flux density, in μmol of
  # light per sq m per sec), rel. humidity (%), and temperature (°C)
  # measured once every 30 min.
  #
  # PPFD is 30 min-averages of values taken every 4 minutes
  # (raw 4 minute values not included).
  #
  # There are two sites, Okutama (site of the independent gametophyte
  # colony) and Uratakao (site of the sporophyte population). Additional
  # sites were also measured, but not included in this analysis
  # because of too much missing data due to mechanical failures.
  #
  # For Okutama, the raw data for temperature, rel. hum,
  # and light (PPFD) are in different columns
  # in the `xlsx` file, so read in each variable separately.

  qc_info = read_excel(
    file_in("data/raw/phenotype_data/QC_sex_eth.xlsx"),
    sheet = 1),
  id_code = read.csv("data/raw/ID_key.csv", header=T),
  bed_conversion = basename(original_plink_data) # convert plink to bed file
  fam = read.table(".fam") # read the .fam file

  run(command="plink",c( "--file", "PSYMETAB_GWAS", "--freq", "--out", "test"), error_on_status=F)



  )
