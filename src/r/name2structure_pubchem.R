library(future)
library(future.apply)
library(progressr)
library(rvest)

source("r/progressr.R")

#' Title
#'
#' @param xs
#'
#' @return
#' @export
#'
#' @examples
name2smiles_pubchem <- function(xs) {
  p <- progressr::progressor(along = xs)
  future.apply::future_lapply(
    future.seed = TRUE,
    X = xs,
    FUN = function(x) {
      tryCatch(
        {
          p(sprintf("x=%g", x))
          cpd <- dataForPubchem[x, "nameCleaned"]
          url <-
            paste0(
              "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/",
              cpd,
              "/property/IsomericSMILES/TXT"
            )
          url <- gsub(
            pattern = "[[:space:]]",
            replacement = "%20",
            x = url
          )
          rvest::read_html(url) |>
            rvest::html_text2()
        },
        error = function(e) {
          NA
        }
      )
    }
  )
}
