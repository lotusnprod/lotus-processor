#######################################################
####################   Functions   ####################
#######################################################

library(dplyr)
# library(ChemmineR)
library(parallel)
library(pbmcapply)
library(readr)
library(rvest)
# library(webchem)

#######################################################
#######################################################

pubchem2inchi <- function(i)
{
  tryCatch({
    cpd <-
      data_translated_pubchem[i, "structure_original_numerical_pubchem"]
    url <-
      paste(
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
        cpd,
        "/property/InChI/txt",
        sep = ""
      )
    url <- gsub(pattern = "\\s",
                replacement = "%20",
                x = url)
    read_html(url) %>%
      html_text()
  }
  , error = function(e) {
    NA
  })
}

#######################################################
#######################################################

name2inchi <- function(i)
  # {
  #   tryCatch({
  #     x <- cts_convert(
  #       query = dataTranslatedNominal[i, "nameCleaned"],
  #       from = "Chemical Name",
  #       to = "InChI Code",
  #       verbose = FALSE,
  #       choices = 1
  #     )
  #     return(x)
  #   }
#   , error = function(e) {
#     NA
#   })
# }
{
  tryCatch({
    cpd <- dataTranslatedNominal[i, "nameCleaned"]
    url <-
      paste("https://cactus.nci.nih.gov/chemical/structure/",
            cpd,
            "/stdinchi",
            sep = "")
    url <- gsub(pattern = "\\s",
                replacement = "%20",
                x = url)
    read_html(url) %>%
      html_text()
  }
  , error = function(e) {
    NA
  })
}

#######################################################
#######################################################

preparing_name <- function(x) {
  x$nameCleaned <- x$structureOriginalNominal
  x$nameCleaned <- gsub("α", "alpha", x$nameCleaned)
  x$nameCleaned <- gsub("Α", "alpha", x$nameCleaned)
  x$nameCleaned <- gsub("β", "beta", x$nameCleaned)
  x$nameCleaned <- gsub("Β", "beta", x$nameCleaned)
  x$nameCleaned <- gsub("γ", "gamma", x$nameCleaned)
  x$nameCleaned <- gsub("Γ", "gamma", x$nameCleaned)
  x$nameCleaned <- gsub("δ", "delta", x$nameCleaned)
  x$nameCleaned <- gsub("Δ", "delta", x$nameCleaned)
  x$nameCleaned <- gsub("ε", "epsilon", x$nameCleaned)
  x$nameCleaned <- gsub("Ε", "epsilon", x$nameCleaned)
  x$nameCleaned <- gsub("- ", "-", x$nameCleaned)
  x$nameCleaned <- gsub("–", "-", x$nameCleaned)
  x$nameCleaned <- gsub("\\) ", "\\)", x$nameCleaned)
  x$nameCleaned <- trimws(x$nameCleaned)
  
  return(x)
}

#######################################################
#######################################################
