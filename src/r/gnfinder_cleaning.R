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
        str_pad(
          string = num,
          width = 6,
          pad = 0
        ),
        ".json",
        sep = ""
      )
  }

  gnfound <- data.frame(fromJSON(
    txt = inpath_gnfinder_f,
    simplifyDataFrame = TRUE
  ))

  data_bio <- vroom_read_safe(path = inpath_organism_f) %>%
    mutate_all(as.character)

  data_bio <- data_bio[!is.na(data_bio[, organismCol]), ]

  data_bio_clean <- biocleaning(
    gnfound = gnfound,
    names = data_bio,
    organismCol = organismCol
  )

  return(data_bio_clean)
}
