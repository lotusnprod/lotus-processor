# title: "editing references"

# loading paths
source("paths.R")

# loading functions
source("functions/reference.R")

# splitting references
###COMMENT### temporary, will write path accordingly if validated that way
source("2_curating/2_editing/reference/subscripts/1_splitting.R")

# translating references
## DOI
source("2_curating/2_editing/reference/subscripts/2_translating/doi.R")

## manual
source("2_curating/2_editing/reference/subscripts/2_translating/manual.R")

## other sources
source("2_curating/2_editing/reference/subscripts/2_translating/other.R")

## pubmed
source("2_curating/2_editing/reference/subscripts/2_translating/pubmed.R")

## text
source("2_curating/2_editing/reference/subscripts/2_translating/text.R")

# integrating references
source("2_curating/2_editing/reference/subscripts/3_integrating.R")

# cleaning references
source("2_curating/2_editing/reference/subscripts/4_cleaning.R")
