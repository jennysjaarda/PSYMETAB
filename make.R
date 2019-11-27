# make.R

# Master script for running PSYMETAB analyses.
# This project uses drake to manage workflows.
# For more information about drake, see
# https://ropensci.github.io/drake/

# Setup ----

# Set working directory
setwd(here::here())

# Load packages
source("code/packages.R")

# Load functions and plans  ----
source("code/functions.R")
source("code/plan.R")

# Run analyses ----
make(plan)
