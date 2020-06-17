# title: "editing chemo"

# loading paths
source("paths.R")

# loading functions
source("functions/chemo.R")

# translating names
###COMMENT### temporary, will write path accordingly if validated that way
source("2_curating/2_editing/chemo/subscripts/1_translating/names.R")

# translating smiles
###will be probably fixed with the  environment, temporarily done so
pathToPython <- "/Users/rutza/opt/anaconda3/bin/python"

system(
  command = paste(
    pathToPython,
    "2_curating/2_editing/chemo/subscripts/1_translating/smiles.py"
  )
)

# curating and enriching structures
# source("2_curating/2_editing/chemo/subscripts/2_curatingAndEnriching/...")

# more enriching ...
