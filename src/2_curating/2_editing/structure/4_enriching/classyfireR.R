# load("../data/interim/temp.Rdata")

cat("This script adds chemical taxonomy to structures dictionary \n")

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")

cat("... libraries \n")
library(classyfireR)
library(pbmcapply)
library(tidyverse)
source("r/vroom_safe.R")

cat("loading files ... \n")
cat("...  counted structures \n")
structureCounted <-
  vroom_read_safe(path = pathDataInterimTablesCleanedStructureStereoCounted)

structuresForClassification <- structureCounted %>%
  distinct(inchikeySanitized)

inchikeys <- structuresForClassification$inchikeySanitized

clasification_list_inchikey <-
  purrr::map(inchikeys, get_classification)

# cat("exporting partial results to Rdata for the moment \n")

# save(list = ls(.GlobalEnv), file = "../data/interim/temp.Rdata")
#
# smiles <- structuresForClassification$smilesSanitized
#
# for (i in seq_len(length(clasification_list_inchikey))) {
#   smiles[[i]] <-
#     ifelse(
#       test = is.null(clasification_list_inchikey[[i]]),
#       yes = smiles[[i]],
#       no = NA
#     )
# }
#
# names(smiles) <- structuresForClassification$inchikeySanitized
#
# smiles <- smiles[!is.na(smiles)]
#
# Sys.sleep(10)
#
# cat("following error is expected because of strange behaviour of the function \n")
# try({
#   classification_list_smiles <- submit_query(label = "query_test",
#                                              input = smiles,
#                                              type = "STRUCTURE")
# })
#
# Sys.sleep(10)
#
# cat("this one too \n ")
# if (exists("classification_list_smiles") != TRUE) {
#   try({
#     classification_list_smiles <- submit_query(label = "query_test",
#                                                input = smiles,
#                                                type = "STRUCTURE")
#   })
# }
#
# Sys.sleep(10)
#
# cat("once again for safety \n ")
# if (exists("classification_list_smiles") != TRUE) {
#   try({
#     classification_list_smiles <- submit_query(label = "query_test",
#                                                input = smiles,
#                                                type = "STRUCTURE")
#   })
# }
#
# Sys.sleep(10)
#
# cat("should work now \n")
# if (exists("classification_list_smiles") != TRUE) {
#   classification_list_smiles <- submit_query(label = "query_test",
#                                              input = smiles,
#                                              type = "STRUCTURE")
# }
cat("worked! \n")

## see accessor methods later on
# classyfireR::alternative_parents(object = clasification_list_inchikey[[1]])
# classyfireR::chebi(object = clasification_list_inchikey[[1]])
# classyfireR::classification(object = clasification_list_inchikey[[1]])
# classyfireR::description(object = clasification_list_inchikey[[1]])
# classyfireR::direct_parent(object = clasification_list_inchikey[[1]])
# classyfireR::meta(object = clasification_list_inchikey[[1]])
# classyfireR::show(object = clasification_list_inchikey[[1]])

X <- seq_along(clasification_list_inchikey)

get_alternative_parents <- function(X) {
  tryCatch(
    {
      alternative_parents <-
        bind_cols(
          "inchikey" = clasification_list_inchikey[[X]]@meta[["inchikey"]],
          "chemontId" = clasification_list_inchikey[[X]]@alternative_parents[["chemont_id"]]
        )
      return(alternative_parents)
    },
    error = function(e) {
      "Error"
    }
  )
}

get_chebi <- function(X) {
  tryCatch(
    {
      chebi <-
        bind_cols(
          "inchikey" = clasification_list_inchikey[[X]]@meta[["inchikey"]],
          "chebi" = clasification_list_inchikey[[X]]@predicted_chebi
        )
      return(chebi)
    },
    error = function(e) {
      "Error"
    }
  )
}

get_direct_parent <- function(X) {
  tryCatch(
    {
      direct_parent <-
        bind_cols(
          "inchikey" = clasification_list_inchikey[[X]]@meta[["inchikey"]],
          "directParent" = clasification_list_inchikey[[X]]@direct_parent[["chemont_id"]]
        )
      return(direct_parent)
    },
    error = function(e) {
      "Error"
    }
  )
}

alternative_parents <- invisible(
  pbmclapply(
    FUN = get_alternative_parents,
    X = X,
    mc.preschedule = TRUE,
    mc.set.seed = TRUE,
    mc.silent = TRUE,
    mc.cores = (parallel::detectCores() - 2),
    mc.cleanup = TRUE,
    mc.allow.recursive = TRUE,
    ignore.interactive = TRUE
  )
)

alternative_parents <-
  bind_rows(alternative_parents[alternative_parents != "Error"]) %>%
  mutate(inchikey = gsub(
    pattern = "InChIKey=",
    replacement = "",
    x = inchikey,
    fixed = TRUE
  ))

chebi <- invisible(
  pbmclapply(
    FUN = get_chebi,
    X = X,
    mc.preschedule = TRUE,
    mc.set.seed = TRUE,
    mc.silent = TRUE,
    mc.cores = (parallel::detectCores() - 2),
    mc.cleanup = TRUE,
    mc.allow.recursive = TRUE,
    ignore.interactive = TRUE
  )
)

chebi <-
  bind_rows(chebi[chebi != "Error"]) %>%
  mutate(inchikey = gsub(
    pattern = "InChIKey=",
    replacement = "",
    x = inchikey,
    fixed = TRUE
  ))

direct_parent <- invisible(
  pbmclapply(
    FUN = get_direct_parent,
    X = X,
    mc.preschedule = TRUE,
    mc.set.seed = TRUE,
    mc.silent = TRUE,
    mc.cores = (parallel::detectCores() - 2),
    mc.cleanup = TRUE,
    mc.allow.recursive = TRUE,
    ignore.interactive = TRUE
  )
)

direct_parent <-
  bind_rows(direct_parent[direct_parent != "Error"]) %>%
  mutate(inchikey = gsub(
    pattern = "InChIKey=",
    replacement = "",
    x = inchikey,
    fixed = TRUE
  ))

cat("exporting Rdata for the moment before deciding what to do \n")

vroom_write_safe_append(
  x = alternative_parents,
  path = "../data/interim/dictionaries/structure/classyfire/alternative_parents.tsv.gz"
)

vroom_write_safe_append(
  x = direct_parent,
  path = "../data/interim/dictionaries/structure/classyfire/direct_parent.tsv.gz"
)

vroom_write_safe_append(
  x = chebi,
  path = "../data/interim/dictionaries/structure/chebi/chebi.tsv.gz"
)

# save(list = ls(.GlobalEnv), file = "../data/interim/temp.Rdata")

end <- Sys.time()

cat("Script finished in", format(end - start), "\n")
