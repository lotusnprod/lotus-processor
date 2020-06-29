source("2_curating/2_editing/bio/functions/biocleaning.R")

gnfinder_cleaning <- function(num, organismCol) {
  if (organismCol == "organismOriginal") {
    inpath_organism_f <- paste(
      pathDataInterimTablesOriginalGnfinder,
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
        pathDataInterimTablesCleanedGnfinderOriginal,
        str_pad(
          string = num,
          width = 6,
          pad = 0
        ),
        ".json",
        sep = "")
  }

  if (organismCol == "organismInterim") {
    inpath_organism_f <- paste(
      pathDataInterimTablesTranslatedGnfinder,
      str_pad(
        string = num,
        width = 6,
        pad = 0
      ),
      ".tsv",
      sep = "")

    inpath_gnfinder_f <-
      paste(
        pathDataInterimTablesCleanedGnfinderTranslated,
        str_pad(
          string = num,
          width = 6,
          pad = 0
        ),
        ".json",
        sep = "")
   }

  gnfound <- data.frame(fromJSON(txt = inpath_gnfinder_f,
                                 simplifyDataFrame = TRUE))

  data_bio <- read_delim(
    file = inpath_organism_f,
    delim = "\t",
    escape_double = FALSE,
    trim_ws = FALSE
  ) %>%
    mutate_all(as.character)


  data_bio <- data_bio[!is.na(data_bio[, organismCol]),]

  data_bio_clean <- biocleaning(gnfound = gnfound,
                                names = data_bio,
                                organismCol = organismCol)

  return(data_bio_clean)
}
