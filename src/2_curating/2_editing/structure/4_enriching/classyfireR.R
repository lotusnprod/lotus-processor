# load("../data/interim/temp.Rdata")
source("r/log_debug.R")
log_debug("This script adds chemical taxonomy to structures dictionary")

start <- Sys.time()

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")

log_debug("... libraries")
library(classyfireR)
library(dplyr)
library(future)
library(future.apply)
library(progressr)
library(purrr)
library(readr)

source("r/progressr.R")

log_debug("opening cache...")
create_dir(export = pathDataInterimDictionariesStructureDictionaryClassyfireDB)
ClassyFireCache <-
  classyfireR::open_cache(dbname = pathDataInterimDictionariesStructureDictionaryClassyfireDB)

log_debug("loading files ...")
log_debug("...  counted structures")
structureCounted <-
  readr::read_delim(
    file = pathDataInterimTablesProcessedStructureStereoCounted,
    delim = "\t",
    col_types = cols(.default = "c")
  )

#' if you want to query again what was not classified
#' in the dictionary
# structureCounted <-
#   readr::read_delim(
#     file = pathDataInterimDictionariesStructureMetadata,
#     delim = "\t",
#     col_types = cols(.default = "c"),
#     locale = locales
#   ) |>
#   dplyr::mutate(inchikeySanitized = structureCleanedInchikey)

if (file.exists(
  pathDataInterimDictionariesStructureDictionaryClassyfireDirectParent
)) {
  old <-
    readr::read_delim(
      file = pathDataInterimDictionariesStructureDictionaryClassyfireDirectParent,
      delim = "\t",
      col_types = cols(.default = "c")
    ) |>
    dplyr::distinct(inchikey)


  structuresForClassification <- structureCounted |>
    dplyr::anti_join(old |>
      dplyr::select(inchikeySanitized = inchikey)) |>
    dplyr::filter(!is.na(inchikeySanitized)) |>
    dplyr::distinct(inchikeySanitized)
} else {
  structuresForClassification <- structureCounted |>
    dplyr::filter(!is.na(inchikeySanitized)) |>
    dplyr::distinct(inchikeySanitized)
}

## if you want to query again what was was classified
## in the dictionary to put it in the local cache
## (in case you lost it)
# structuresForClassification <- old |>
#   dplyr::select(inchikeySanitized = inchikey) |>
#   dplyr::filter(!is.na(inchikeySanitized)) |>
#   dplyr::distinct(inchikeySanitized)

# RSQLite::dbListObjects(ClassyFireCache)

## If you want to put things you have in your cache to your dictionary
# incache <- RSQLite::dbGetQuery(conn = ClassyFireCache,
#                                statement = "SELECT *
#   FROM classyfire;")
# structuresForClassification <- dplyr::inner_join(
#   structuresForClassification,
#   incache,
#   by = c("inchikeySanitized" = "InChikey"))

## if you want to query again what was was classified
## in the dictionary to put it in the local cache
## (in case you lost it)
# incache <- RSQLite::dbGetQuery(conn = ClassyFireCache,
#                                statement = "SELECT *
#   FROM classyfire;")
# structuresForClassification <- dplyr::anti_join(
#   structuresForClassification,
#   incache,
#   by = c("inchikeySanitized" = "InChikey"))

inchikeys <- structuresForClassification$inchikeySanitized

log_debug("classifying...")
clasification_list_inchikey <-
  purrr::map(inchikeys, get_classification, conn = ClassyFireCache)

log_debug("worked!")

## see accessor methods later on
# classyfireR::alternative_parents(object = clasification_list_inchikey[[1]])
# classyfireR::chebi(object = clasification_list_inchikey[[1]])
# classyfireR::classification(object = clasification_list_inchikey[[1]])
# classyfireR::description(object = clasification_list_inchikey[[1]])
# classyfireR::direct_parent(object = clasification_list_inchikey[[1]])
# classyfireR::meta(object = clasification_list_inchikey[[1]])
# classyfireR::show(object = clasification_list_inchikey[[1]])

xs <- seq_along(clasification_list_inchikey)

get_alternative_parents <- function(xs) {
  p <- progressr::progressor(along = xs)
  future.apply::future_lapply(
    future.seed = TRUE,
    X = xs,
    FUN = function(x) {
      tryCatch(
        {
          p(sprintf("x=%g", as.numeric(x))) ## little hack
          alternative_parents <-
            dplyr::bind_cols(
              "inchikey" = clasification_list_inchikey[[x]]@meta[["inchikey"]],
              "chemontId" = clasification_list_inchikey[[x]]@alternative_parents[["chemont_id"]]
            )
          return(alternative_parents)
        },
        error = function(e) {
          "Error"
        }
      )
    }
  )
}

get_chebi <- function(xs) {
  p <- progressr::progressor(along = xs)
  future.apply::future_lapply(
    future.seed = TRUE,
    X = xs,
    FUN = function(x) {
      tryCatch(
        {
          p(sprintf("x=%g", as.numeric(x))) ## little hack
          chebi <-
            dplyr::bind_cols(
              "inchikey" = clasification_list_inchikey[[x]]@meta[["inchikey"]],
              "chebi" = clasification_list_inchikey[[x]]@predicted_chebi
            )
          return(chebi)
        },
        error = function(e) {
          "Error"
        }
      )
    }
  )
}

get_direct_parent <- function(xs) {
  p <- progressr::progressor(along = xs)
  future.apply::future_lapply(
    future.seed = TRUE,
    X = xs,
    FUN = function(x) {
      tryCatch(
        {
          p(sprintf("x=%g", as.numeric(x))) ## little hack
          direct_parent <-
            dplyr::bind_cols(
              "inchikey" = clasification_list_inchikey[[x]]@meta[["inchikey"]],
              "directParent" = clasification_list_inchikey[[x]]@direct_parent[["chemont_id"]]
            )
          return(direct_parent)
        },
        error = function(e) {
          "Error"
        }
      )
    }
  )
}

alternative_parents <- get_alternative_parents(xs = xs) |>
  progressr::with_progress()

if (!is_empty(alternative_parents)) {
  alternative_parents <-
    dplyr::bind_rows(alternative_parents[alternative_parents != "Error"])
} else {
  alternative_parents <- data.frame()
}

if (nrow(alternative_parents != 0)) {
  alternative_parents <- alternative_parents |>
    dplyr::mutate(inchikey = gsub(
      pattern = "InChIKey=",
      replacement = "",
      x = inchikey,
      fixed = TRUE
    ))
}

chebi <- get_chebi(xs = xs) |>
  progressr::with_progress()

if (!is_empty(chebi)) {
  chebi <- dplyr::bind_rows(chebi[chebi != "Error"])
} else {
  chebi <- data.frame()
}

if (nrow(chebi != 0)) {
  chebi <- chebi |>
    dplyr::mutate(inchikey = gsub(
      pattern = "InChIKey=",
      replacement = "",
      x = inchikey,
      fixed = TRUE
    ))
}

direct_parent <- get_direct_parent(xs = xs) |>
  progressr::with_progress()

if (!is_empty(direct_parent)) {
  direct_parent <-
    dplyr::bind_rows(direct_parent[direct_parent != "Error"])
} else {
  direct_parent <- data.frame()
}

if (nrow(direct_parent != 0)) {
  direct_parent <- direct_parent |>
    dplyr::mutate(inchikey = gsub(
      pattern = "InChIKey=",
      replacement = "",
      x = inchikey,
      fixed = TRUE
    ))
}

log_debug("exporting")
create_dir(export = pathDataInterimDictionariesStructureDictionaryChebiFile)

if (file.exists(
  pathDataInterimDictionariesStructureDictionaryClassyfireAlternativeParent
)) {
  readr::write_delim(
    x = alternative_parents,
    delim = "\t",
    file = pathDataInterimDictionariesStructureDictionaryClassyfireAlternativeParent,
    na = "",
    append = TRUE
  )
} else {
  readr::write_delim(
    x = alternative_parents,
    delim = "\t",
    file = pathDataInterimDictionariesStructureDictionaryClassyfireAlternativeParent,
    na = ""
  )
}

if (file.exists(pathDataInterimDictionariesStructureDictionaryClassyfireDirectParent)) {
  readr::write_delim(
    x = direct_parent,
    delim = "\t",
    file = pathDataInterimDictionariesStructureDictionaryClassyfireDirectParent,
    na = "",
    append = TRUE
  )
} else {
  readr::write_delim(
    x = direct_parent,
    delim = "\t",
    file = pathDataInterimDictionariesStructureDictionaryClassyfireDirectParent,
    na = ""
  )
}

if (file.exists(pathDataInterimDictionariesStructureDictionaryChebiFile)) {
  readr::write_delim(
    x = chebi,
    delim = "\t",
    file = pathDataInterimDictionariesStructureDictionaryChebiFile,
    na = "",
    append = TRUE
  )
} else {
  readr::write_delim(
    x = chebi,
    delim = "\t",
    file = pathDataInterimDictionariesStructureDictionaryChebiFile,
    na = ""
  )
}

# save(list = ls(.GlobalEnv), file = "../data/interim/temp.Rdata")

end <- Sys.time()

log_debug("Script finished in", format(end - start))
