library(processx)
library(pacman)

library(dplyr)
library(drake)
library(ggplot2)
library(knitr)
library(purrr)
library(rlang)
library(tibble)
library(processx)

library(pacman)
library(clustermq)
Packages <- c("dplyr", "data.table", "R.utils", "rbgen","devtools",
"TwoSampleMR","MendelianRandomization","MRInstruments", "rmeta","mr.raps","ukbtools","rslurm",
"stringr","bit64","optparse", "reader","foreach","doParallel","parallel","MASS","hyperSpec","tidyverse",
"ggthemes","tryCatchLog", "futile.logger","stats4","bbmle","readxl", "qqman", "taRifx", "DataExplorer")
p_load(Packages, character.only=T)
