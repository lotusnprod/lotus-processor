# title: "editing bio"

# loading paths
source("paths.R")

# loading functions
source("functions/bio.R")

gnfinder_cleaning <- function(num, organismCol) {
  if (organismCol == "organismOriginal") {
    inpath_organism_f <- paste(
      pathOriginalOrganismDistinct,
      "originalOrganismGnfinderUntil_",
      num,
      ".tsv",
      sep = "")

    inpath_gnfinder_f <-
      paste(
        pathCleanedOrganismOriginalDirJson,
        "originalOrganismGnfinderUntil_",
        num,
        ".json",
        sep = "")
  }

  if (organismCol == "organismInterim") {
    inpath_organism_f <- paste(
      pathTranslatedOrganismDistinct,
      "translatedOrganismGnfinderUntil_",
      num,
      ".tsv",
      sep = "")

    inpath_gnfinder_f <-
      paste(
        pathCleanedOrganismTranslatedDirJson,
        "translatedOrganismGnfinderUntil_",
        num,
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

  data_bio_clean <- biocleaning(x = gnfound,
                                y = data_bio,
                                organismCol = organismCol)

  return(data_bio_clean)
}

# cleaning original organism
###COMMENT### temporary, will write path accordingly if validated that way
source("2_curating/2_editing/bio/subscripts/1_cleaningOriginalOrganism.R")

# translating organism
source("2_curating/2_editing/bio/subscripts/2_translatingOrganism.R")

# cleaning translated organism 
source("2_curating/2_editing/bio/subscripts/3_cleaningTranslatedOrganism.R")

# cleaning taxonomy
source("2_curating/2_editing/bio/subscripts/4_cleaningTaxonomy.R")
