# title: "editing bio"

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
source("2_curating/2_editing/organism/functions/biocleaning.R")
source("2_curating/2_editing/organism/functions/gnfinder_cleaning.R")
source("2_curating/2_editing/organism/functions/manipulating_taxo.R")
#######################################################
#######################################################
log_debug("Finished loading functions")


# cleaning original organism
###COMMENT### temporary, will write path accordingly if validated that way
source("2_curating/2_editing/organism/subscripts/1_cleaningOriginal.R")

# translating organism
source("2_curating/2_editing/organism/subscripts/2_translating.R")

# cleaning translated organism 
source("2_curating/2_editing/organism/subscripts/3_cleaningTranslated.R")

# cleaning taxonomy
source("2_curating/2_editing/organism/subscripts/4_cleaningTaxonomy.R")
