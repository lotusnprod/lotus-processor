source("paths.R")
library(future)
library(future.apply)
library(jsonlite)
library(progressr)
library(rvest)
library(purrr)

source("r/progressr.R")

#' Title
#'
#' @param xs
#'
#' @return
#' @export
#'
#' @examples
getClass <- function(xs) {
  p <- progressr::progressor(along = xs)
  future.apply::future_lapply(
    future.seed = TRUE,
    X = xs,
    FUN = function(x) {
      tryCatch(
        {
          p(sprintf("x=%g", x))
          url_id <- paste0(url, order, queries[x])
          result <- rvest::read_html(url_id) |>
            rvest::html_text() |>
            jsonlite::fromJSON()

          df <- data.frame(
            structure_smiles_2D = new$structure_smiles_2D[[x]],
            query = new$query[[x]],
            pathway = if (!is_empty(result$pathway_results)) {
              result$pathway_results
            } else {
              NA
            },
            superclass = if (!is_empty(result$superclass_results)) {
              result$superclass_results
            } else {
              NA
            },
            class = if (!is_empty(result$class_results)) {
              result$class_results
            } else {
              NA
            },
            glycoside = if (!is_empty(result$isglycoside)) {
              result$isglycoside
            } else {
              NA
            }
          )

          return(df)
        },
        error = function(e) {
          df <- data.frame(
            structure_smiles_2D = new$structure_smiles_2D[[x]],
            query = new$query[[x]],
            pathway = NA,
            superclass = NA,
            class = NA,
            glycoside = NA
          )
          return(df)
        }
      )
    }
  )
}
