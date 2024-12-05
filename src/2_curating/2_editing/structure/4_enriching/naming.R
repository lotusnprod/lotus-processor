source("r/log_debug.R")
log_debug("This script adds chemical names (and CIDs) to structures dictionary")

start <- Sys.time()

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")

log_debug("... libraries")
library(dplyr)
library(httr)
library(readr)

log_debug("loading files ...")
log_debug("...  counted structures")
structureCounted <-
  readr::read_delim(
    file = pathDataInterimTablesProcessedStructureStereoCounted,
    delim = "\t"
  )

log_debug("keeping inchikey only ...")
inchikeyDictionary <- structureCounted |>
  dplyr::distinct(inchikeySanitized) |>
  dplyr::select(inchikey = inchikeySanitized)
## if some names are missing in structures metadata
# inchikeyDictionary <- readr::read_delim(
#     file = pathDataInterimDictionariesStructureMetadata,
#     delim = "\t",
#     col_types = cols(.default = "c"),
#     locale = locales
#   ) |>
#   dplyr::filter(is.na(structureCleaned_nameTraditional)) |>
#   dplyr::distinct(structureCleanedInchikey) |>
#   dplyr::select(inchikey = structureCleanedInchikey)

chunk <- 5000
n <- nrow(inchikeyDictionary)
r <- rep(1:ceiling(n / chunk), each = chunk)[1:n]
lists <- split(inchikeyDictionary, r)

lists_ready <-
  lapply(
    X = lists,
    FUN = function(x) {
      paste0(
        "inchikey=",
        paste(
          x |>
            dplyr::pull() |>
            as.character(),
          collapse = ","
        )
      )
    }
  )

url <-
  "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/property/Title,IUPACName,inchikey/JSON"

post_to_pubchem <- function(url, vals) {
  content <- httr::RETRY(
    verb = "POST",
    url = url,
    body = vals
  ) |>
    httr::content()

  df <- content[["PropertyTable"]][["Properties"]] |>
    dplyr::bind_rows()

  return(df)
}

list_results <-
  lapply(
    X = lists_ready,
    FUN = post_to_pubchem,
    url = url
  )


empty_df <- tibble(
  inchikeySanitized = character(),
  structureCleaned_nameTraditional = character(),
  structureCleaned_nameIupac = character(),
  structureCleaned_cid = character()
)

results <- tryCatch(
  {
    list_results |>
      dplyr::bind_rows() |>
      dplyr::distinct() |>
      dplyr::select(
        inchikeySanitized = InChIKey,
        structureCleaned_nameTraditional = Title,
        structureCleaned_nameIupac = IUPACName,
        structureCleaned_cid = CID
      )
  },
  error = function(e) {
    message("An error occurred, returning an empty dataframe: ", e$message)
    empty_df
  }
)

## if some names are missing in structures metadata
# part_1 <- structureMetadata |>
#   dplyr::filter(is.na(structureCleaned_nameTraditional)) |>
#   dplyr::select(-structureCleaned_nameTraditional,
#          -structureCleaned_nameIupac) |>
#   dplyr::left_join(results, by = c("structureCleanedInchikey" = "inchikeySanitized"))
# part_2 <- structureMetadata |>
#   dplyr::filter(!is.na(structureCleaned_nameTraditional))
# structureMetadata <- rbind(part_1, part_2)

structureNamed_cleaned <-
  dplyr::left_join(
    structureCounted,
    results
  )

## old legacy version with molconvert which generated a lot of (+/-) without stereo
# structureNamed_defined <- structureNamed |>
#   dplyr::filter(count_unspecified_atomic_stereocenters == 0) |>
#   dplyr::mutate_all(as.character)
#
# structureNamed_undefined <- structureNamed |>
#   dplyr::filter(count_unspecified_atomic_stereocenters != 0) |>
#   dplyr::mutate(
#     structureCleaned_nameTraditional = gsub(
#       pattern = "(-)-",
#       replacement = "",
#       x = structureCleaned_nameTraditional,
#       fixed = TRUE
#     )
#   ) |>
#   dplyr::mutate(
#     structureCleaned_nameTraditional = gsub(
#       pattern = "(+)-",
#       replacement = "",
#       x = structureCleaned_nameTraditional,
#       fixed = TRUE
#     )
#   ) |>
#   dplyr::mutate(
#     structureCleaned_nameTraditional = gsub(
#       pattern = "(-)",
#       replacement = "",
#       x = structureCleaned_nameTraditional,
#       fixed = TRUE
#     )
#   ) |>
#   dplyr::mutate(
#     structureCleaned_nameTraditional = gsub(
#       pattern = "(+)",
#       replacement = "",
#       x = structureCleaned_nameTraditional,
#       fixed = TRUE
#     )
#   ) |>
#   dplyr::mutate_all(as.character)
#
# structureNamed_cleaned <-
#   dplyr::bind_rows(structureNamed_defined, structureNamed_undefined)

if (mode == "custom") {
  library(tidyr)
  structureNamed <-
    readr::read_delim(
      file = pathDataInterimTablesOriginalStructureFull,
      delim = "\t"
    ) %>%
    ## using https://stackoverflow.com/a/58837832 to allow multiple smiles/inchi/nominal
    dplyr::group_by(structureType) %>%
    dplyr::mutate(row = dplyr::row_number()) %>%
    tidyr::pivot_wider(
      names_from = structureType,
      values_from = structureValue
    ) %>%
    tidyr::unnest(cols = c(inchi, smiles, nominal)) %>%
    tidyr::pivot_longer(cols = colnames(.)[grepl(pattern = "inchi|smiles", x = colnames(.))]) %>%
    dplyr::distinct(
      structureType = name,
      structureValue = value,
      structureCleaned_nameTraditional = nominal
    ) %>%
    dplyr::filter(!is.na(structureValue))

  structureTranslated <-
    readr::read_delim(file = pathDataInterimTablesTranslatedStructureFinal)

  structureNamed_cleaned <- structureNamed_cleaned |>
    dplyr::select(-structureCleaned_nameTraditional) |>
    dplyr::left_join(structureTranslated) |>
    dplyr::left_join(structureNamed) |>
    dplyr::select(-structureType, -structureValue)
}

log_debug("ensuring directories exist")
log_debug("exporting ...")
log_debug(pathDataInterimTablesProcessedStructureNamed)

readr::write_delim(
  x = structureNamed_cleaned,
  delim = "\t",
  file = pathDataInterimTablesProcessedStructureNamed,
  na = ""
)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
