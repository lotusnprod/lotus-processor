# load("../data/interim/temp.Rdata")
source("r/log_debug.R")
source("r/parallel.R")
log_debug("This script adds chemical taxonomy to structures dictionary")

start <- Sys.time()

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")

log_debug("... libraries")
library(classyfireR)
library(dplyr)
library(pbmcapply)
library(purrr)
library(readr)

log_debug("loading files ...")
log_debug("...  counted structures")
structureCounted <-
  read_delim(
    file = pathDataInterimTablesProcessedStructureStereoCounted,
    delim = "\t"
  )

old <-
  read_delim(
    file = file.path(
      pathDataInterimDictionariesStructureDictionaryClassyfire,
      "direct_parent.tsv.gz"
    ),
    delim = "\t"
  ) %>%
  distinct(inchikey)

structuresForClassification <-
  anti_join(
    structureCounted,
    old %>% select(inchikeySanitized = inchikey)
  ) %>%
  filter(!is.na(inchikeySanitized)) %>%
  distinct(inchikeySanitized)

inchikeys <- structuresForClassification$inchikeySanitized

clasification_list_inchikey <-
  purrr::map(inchikeys, get_classification)

log_debug("worked!")

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
    mc.cores = numCores,
    mc.cleanup = TRUE,
    mc.allow.recursive = TRUE,
    ignore.interactive = TRUE,
    mc.style = "txt",
    mc.substyle = 1
  )
)

alternative_parents <-
  bind_rows(alternative_parents[alternative_parents != "Error"])

if (nrow(alternative_parents != 0)) {
  alternative_parents <- alternative_parents %>%
    mutate(inchikey = gsub(
      pattern = "InChIKey=",
      replacement = "",
      x = inchikey,
      fixed = TRUE
    ))
}

chebi <- invisible(
  pbmclapply(
    FUN = get_chebi,
    X = X,
    mc.preschedule = TRUE,
    mc.set.seed = TRUE,
    mc.cores = numCores,
    mc.cleanup = TRUE,
    mc.allow.recursive = TRUE,
    ignore.interactive = TRUE,
    mc.style = "txt",
    mc.substyle = 1
  )
)

chebi <- bind_rows(chebi[chebi != "Error"])

if (nrow(chebi != 0)) {
  chebi <- chebi %>%
    mutate(inchikey = gsub(
      pattern = "InChIKey=",
      replacement = "",
      x = inchikey,
      fixed = TRUE
    ))
}

direct_parent <- invisible(
  pbmclapply(
    FUN = get_direct_parent,
    X = X,
    mc.preschedule = TRUE,
    mc.set.seed = TRUE,
    mc.cores = numCores,
    mc.cleanup = TRUE,
    mc.allow.recursive = TRUE,
    ignore.interactive = TRUE,
    mc.style = "txt",
    mc.substyle = 1
  )
)

direct_parent <-
  bind_rows(direct_parent[direct_parent != "Error"])

if (nrow(direct_parent != 0)) {
  direct_parent <- direct_parent %>%
    mutate(inchikey = gsub(
      pattern = "InChIKey=",
      replacement = "",
      x = inchikey,
      fixed = TRUE
    ))
}

log_debug("exporting")

if (file.exists("../data/interim/dictionaries_full/structure/classyfire/alternative_parents.tsv.gz")) {
  write_delim(
    x = alternative_parents,
    delim = "\t",
    file = "../data/interim/dictionaries_full/structure/classyfire/alternative_parents.tsv.gz",
    na = "",
    append = TRUE
  )
} else {
  write_delim(
    x = alternative_parents,
    delim = "\t",
    file = "../data/interim/dictionaries_full/structure/classyfire/alternative_parents.tsv.gz",
    na = ""
  )
}

if (file.exists("../data/interim/dictionaries_full/structure/classyfire/direct_parent.tsv.gz")) {
  write_delim(
    x = direct_parent,
    delim = "\t",
    file = "../data/interim/dictionaries_full/structure/classyfire/direct_parent.tsv.gz",
    na = "",
    append = TRUE
  )
} else {
  write_delim(
    x = direct_parent,
    delim = "\t",
    file = "../data/interim/dictionaries_full/structure/classyfire/direct_parent.tsv.gz",
    na = ""
  )
}

if (file.exists("../data/interim/dictionaries_full/structure/chebi/chebi.tsv.gz")) {
  write_delim(
    x = chebi,
    delim = "\t",
    file = "../data/interim/dictionaries_full/structure/chebi/chebi.tsv.gz",
    na = "",
    append = TRUE
  )
} else {
  write_delim(
    x = chebi,
    delim = "\t",
    file = "../data/interim/dictionaries_full/structure/chebi/chebi.tsv.gz",
    na = ""
  )
}

# save(list = ls(.GlobalEnv), file = "../data/interim/temp.Rdata")

end <- Sys.time()

log_debug("Script finished in", format(end - start))
