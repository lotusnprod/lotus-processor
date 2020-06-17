# title: "editing bio"

#library(Hmisc)
library(data.table)
library(dplyr)
library(jsonlite)
library(parallel)
library(pbmcapply)
library(readr)
library(splitstackshape)
library(stringi)
library(stringr)
library(zoo)
library(tidyr)

# loading paths
source("paths.R")

source("functions/log.R")
source("functions/helpers.R")
source("functions/parallel.R")

log_debug("Ready to start")

log_debug("Creating biocleaning")
#######################################################
#######################################################
source("2_curating/2_editing/bio/functions/biocleaning.R")
source("2_curating/2_editing/bio/functions/gnfinder_cleaning.R")
source("2_curating/2_editing/bio/functions/manipulating_taxo.R")
#######################################################
#######################################################
log_debug("Finished loading functions")


# cleaning original organism
###COMMENT### temporary, will write path accordingly if validated that way
source("2_curating/2_editing/bio/subscripts/1_cleaningOriginalOrganism.R")

# translating organism
source("2_curating/2_editing/bio/subscripts/2_translatingOrganism.R")

# cleaning translated organism 
source("2_curating/2_editing/bio/subscripts/3_cleaningTranslatedOrganism.R")

# cleaning taxonomy
source("2_curating/2_editing/bio/subscripts/4_cleaningTaxonomy.R")
