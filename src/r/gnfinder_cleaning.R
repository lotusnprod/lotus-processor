source("paths.R")
library(jsonlite)
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
  if (organismCol == "organismOriginal") {
    inpath_organism_f <- paste(
      pathDataInterimTablesOriginalOrganism,
      "/",
      str_pad(
        string = num,
        width = 6,
        pad = 0
      ),
      ".tsv",
      sep = ""
    )

    inpath_gnfinder_f <-
      paste(
        pathDataInterimTablesCleanedOrganismOriginal,
        "/",
        str_pad(
          string = num,
          width = 6,
          pad = 0
        ),
        ".json",
        sep = ""
      )
  }

  if (organismCol == "organismInterim") {
    inpath_organism_f <- paste(
      pathDataInterimTablesTranslatedOrganism,
      "/",
      str_pad(
        string = num,
        width = 6,
        pad = 0
      ),
      ".tsv",
      sep = ""
    )

    inpath_gnfinder_f <-
      paste(
        pathDataInterimTablesCleanedOrganismTranslated,
        "/",
        str_pad(
          string = num,
          width = 6,
          pad = 0
        ),
        ".json",
        sep = ""
      )
  }

  data_bio <- vroom(
    file = inpath_organism_f,
    delim = "\t",
    escape_double = TRUE,
    escape_backslash = TRUE,
    col_types = cols(.default = "c"),
    num_threads = 1
  ) %>%
    mutate_all(as.character)

  ## because gnfinder reads it that way for counting chars
  data_bio_2 <- vroom(
    file = inpath_organism_f,
    delim = "\t",
    quote = "",
    escape_double = FALSE,
    escape_backslash = FALSE,
    col_types = cols(.default = "c"),
    num_threads = 1
  ) %>%
    mutate_all(as.character)

  data_bio <- data_bio[!is.na(data_bio[, organismCol]), ]

  data_bio_2 <- data_bio_2[!is.na(data_bio_2[, paste0("\"", organismCol, "\"")]), ]

  if (fromJSON(
    txt = inpath_gnfinder_f,
    simplifyDataFrame = TRUE
  )$metadata$totalNames != 0) {
    gnfound <- data.frame(fromJSON(
      txt = inpath_gnfinder_f,
      simplifyDataFrame = TRUE
    ))

    data_bio_clean <- biocleaning(
      gnfound = gnfound,
      names = data_bio,
      names_quotes = data_bio_2,
      organismCol = organismCol
    )
  }
  else {
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
