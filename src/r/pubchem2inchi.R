source("paths.R")
library(rvest)

#' Title
#'
#' @param i
#'
#' @return
#' @export
#'
#' @examples
pubchem2inchi <- function(i) {
  tryCatch(
    {
      cpd <-
        data_translated_pubchem[i, "structure_original_numerical_pubchem"]
      url <-
        paste(
          "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
          cpd,
          "/property/InChI/txt",
          sep = ""
        )
      url <- gsub(
        pattern = "\\s",
        replacement = "%20",
        x = url
      )
      read_html(url) %>%
        html_text()
    },
    error = function(e) {
      NA
    }
  )
}
