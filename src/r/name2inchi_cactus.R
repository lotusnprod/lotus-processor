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
name2inchi_cactus <- function(i) {
  tryCatch(
    {
      cpd <- dataForCactus[i, "nameCleaned"]
      url <-
        paste0("https://cactus.nci.nih.gov/chemical/structure/", cpd, "/stdinchi")
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
