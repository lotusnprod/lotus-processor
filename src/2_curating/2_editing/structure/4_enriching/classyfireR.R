cat("This script adds chemical taxonomy to structures dictionary \n")

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")

cat("... libraries \n")
library(classyfireR)
library(tidyverse)
source("r/vroom_safe.R")

cat("loading files ... \n")
cat("...  counted structures \n")
structureCounted <-
  vroom_read_safe(path = pathDataInterimTablesCleanedStructureStereoCounted)

structuresForClassification <- structureCounted %>%
  distinct(smilesSanitized, inchikeySanitized) %>%
  sample_n(1000)

inchikeys <- structuresForClassification$inchikeySanitized

clasification_list_inchikey <-
  purrr::map(inchikeys, get_classification)

cat("exporting partial results to Rdata for the moment \n")

save(list = ls(.GlobalEnv), file = "../data/interim/temp.Rdata")

smiles <- structuresForClassification$smilesSanitized

for (i in seq_len(length(clasification_list_inchikey))) {
  smiles[[i]] <-
    ifelse(
      test = is.null(clasification_list_inchikey[[i]]),
      yes = smiles[[i]],
      no = NA
    )
}

names(smiles) <- structuresForClassification$inchikeySanitized

smiles <- smiles[!is.na(smiles)]

Sys.sleep(10)

cat("following error is expected because of strange behaviour of the function \n")
try({
  classification_list_smiles <- submit_query(
    label = "query_test",
    input = smiles,
    type = "STRUCTURE"
  )
})

Sys.sleep(10)

cat("this one too \n ")
if (exists("classification_list_smiles") != TRUE) {
  try({
    classification_list_smiles <- submit_query(
      label = "query_test",
      input = smiles,
      type = "STRUCTURE"
    )
  })
}

Sys.sleep(10)

cat("once again for safety \n ")
if (exists("classification_list_smiles") != TRUE) {
  try({
    classification_list_smiles <- submit_query(
      label = "query_test",
      input = smiles,
      type = "STRUCTURE"
    )
  })
}

Sys.sleep(10)

cat("should work now \n")
if (exists("classification_list_smiles") != TRUE) {
  classification_list_smiles <- submit_query(
    label = "query_test",
    input = smiles,
    type = "STRUCTURE"
  )
}
cat("worked! \n")

## see accessor methods later on

cat("exporting Rdata for the moment before deciding what to do \n")

save(list = ls(.GlobalEnv), file = "../data/interim/temp.Rdata")

end <- Sys.time()

cat("Script finished in", format(end - start), "\n")
