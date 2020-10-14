#######################################################
####################   Functions   ####################
#######################################################

library(parallel)
library(pbmcapply)
library(rvest)
library(tidyverse)
library(webchem)

#######################################################
#######################################################

# pubchem2inchi <- function(i)
# {
#   tryCatch({
#     cpd <-
#       data_translated_pubchem[i, "structure_original_numerical_pubchem"]
#     url <-
#       paste(
#         "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
#         cpd,
#         "/property/InChI/txt",
#         sep = ""
#       )
#     url <- gsub(pattern = "\\s",
#                 replacement = "%20",
#                 x = url)
#     read_html(url) %>%
#       html_text()
#   }
#   , error = function(e) {
#     NA
#   })
# }

#######################################################
#######################################################

name2inchi_cactus <- function(i)
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
  x$nameCleaned <- x$structureOriginal_nominal
  x$nameCleaned <- gsub("\\u03b1", "alpha", x$nameCleaned,
                        fixed = TRUE)
  x$nameCleaned <- gsub("\\u03b2", "beta", x$nameCleaned,
                        fixed = TRUE)
  x$nameCleaned <- gsub("\\u03b3", "gamma", x$nameCleaned,
                        fixed = TRUE)
  x$nameCleaned <- gsub("\\u03b4", "delta", x$nameCleaned,
                        fixed = TRUE)
  x$nameCleaned <- gsub("\\u03b5", "epsilon", x$nameCleaned,
                        fixed = TRUE)
  x$nameCleaned <- gsub("\\u03c9", "omega", x$nameCleaned,
                        fixed = TRUE)
  x$nameCleaned <- gsub("(\\u00b1)-", "", x$nameCleaned,
                        fixed = TRUE)
  x$nameCleaned <- gsub("\\u00f6", "ö", x$nameCleaned,
                        fixed = TRUE)
  x$nameCleaned <- gsub("\\u2192", "->", x$nameCleaned,
                        fixed = TRUE)
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

name2inchi_cts <- function(i)
{
  tryCatch({
    x <- cts_convert(
      query = dataTranslatedNominal_cts[i, "nameCleaned"],
      from = "Chemical Name",
      to = "InChI Code",
      verbose = FALSE,
      match = "first"
    )[[1]]
    return(x)
  }
  , error = function(e) {
    NA
  })
}

#######################################################
#######################################################
