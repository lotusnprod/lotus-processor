# title: "editing bio"

# loading paths
source("paths.R")

# loading functions
source("functions/bio.R")

# cleaning original organism
###COMMENT### temporary, will write path accordingly if validated that way
source("2_curating/2_editing/bio/subscripts/1_cleaningOriginalOrganism.R")

# translating organism
source("2_curating/2_editing/bio/subscripts/2_translatingOrganism.R")

# cleaning translated organism 
source("2_curating/2_editing/bio/subscripts/3_cleaningTranslatedOrganism.R")

# cleaning taxonomy
source("2_curating/2_editing/bio/subscripts/4_cleaningTaxonomy.R")
