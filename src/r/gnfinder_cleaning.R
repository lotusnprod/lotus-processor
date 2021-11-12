source("paths.R")
library(jsonlite)
library(readr)
library(stringr)
source("r/biocleaning.R")

#' Title
#'
#' @param num
#' @param organismCol
#'
#' @return
#' @export
#'
#' @examples
gnfinder_cleaning <- function(num, organismCol) {
  if (organismCol == "organismValue") {
    inpath_organism_f <-
      paste0(
        pathDataInterimTablesOriginalOrganism,
        "/",
        str_pad(
          string = num,
          width = 6,
          pad = 0
        ),
        ".tsv"
      )

    inpath_gnfinder_f <-
      paste0(
        pathDataInterimTablesProcessedOrganismOriginal,
        "/",
        str_pad(
          string = num,
          width = 6,
          pad = 0
        ),
        ".json"
      )
  }

  if (organismCol == "organismInterim") {
    inpath_organism_f <-
      paste0(
        pathDataInterimTablesTranslatedOrganism,
        "/",
        str_pad(
          string = num,
          width = 6,
          pad = 0
        ),
        ".tsv"
      )

    inpath_gnfinder_f <-
      paste0(
        pathDataInterimTablesProcessedOrganismTranslated,
        "/",
        str_pad(
          string = num,
          width = 6,
          pad = 0
        ),
        ".json"
      )
  }

  ## because gnfinder reads it that way for counting chars
  data_bio <- read_delim(
    file = inpath_organism_f,
    delim = "\t",
    quote = "",
    escape_double = FALSE,
    trim_ws = FALSE,
    escape_backslash = FALSE,
    na = "",
    col_types = cols(.default = "c")
  ) %>%
    mutate_all(as.character)

  if (fromJSON(
    txt = inpath_gnfinder_f,
    simplifyDataFrame = TRUE
  )$metadata$totalNames != 0) {
    gnfound <- data.frame(fromJSON(
      txt = inpath_gnfinder_f,
      simplifyDataFrame = TRUE
    ))
  } else {
    gnfound <- data.frame(
      names.start = 0,
      names.verification = NA
    )
    b <- data.frame(preferredResults = NA)
    c <- data.frame(
      list(
        datasourceId = NA,
        dataSourceTitleShort = NA,
        curation = NA,
        recordId = NA,
        outlink = NA,
        entryDate = NA,
        matchedName = NA,
        matchedCardinality = NA,
        matchCanonicalSimple = NA,
        matchedCanonicalFull = NA,
        currentRecordId = NA,
        currentName = NA,
        currentCardinality = NA,
        currentCanonicalSimple = NA,
        currentCanonicalFull = NA,
        isSynonym = NA,
        classificationPath = NA,
        classificationRanks = NA,
        classificationIds = NA,
        editDistance = NA,
        stemEditDistance = NA,
        matchType = NA,
        localId = NA
      )
    )
    b$preferredResults <- c
    gnfound$names.verification <- b
  }

  if (nrow(gnfound) != 0) {
    data_bio_clean <- biocleaning(
      gnfound = gnfound,
      names = data_bio,
      organismCol = organismCol
    )
  } else {
    data_bio_clean <- data_bio %>%
      mutate(
        nchar = NA,
        sum = NA,
        value_min = NA,
        value_max = NA,
        canonicalname = NA,
        canonicalnameCurrent = NA,
        taxonId = NA,
        dbTaxo = NA,
        taxonomy = NA,
        rank = NA,
        ids = NA,
        dbQuality = NA
      )
  }
  return(data_bio_clean)
}
