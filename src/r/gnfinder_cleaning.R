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
        pathDataInterimTablesCleanedOrganismOriginal,
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
        pathDataInterimTablesCleanedOrganismTranslated,
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
  )$
    metadata$
    totalNames != 0) {
    gnfound <- data.frame(fromJSON(
      txt = inpath_gnfinder_f,
      simplifyDataFrame = TRUE
    ))
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
