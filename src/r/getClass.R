source("paths.R")
library(jsonlite)
library(rvest)
library(purrr)

#' Title
#'
#' @param X
#'
#' @return
#' @export
#'
#' @examples
getClass <- function(X) {
  tryCatch(
    {
      url_id <- paste0(url, order, queries[X])
      result <- read_html(url_id) %>%
        html_text() %>%
        fromJSON()

      df <- data.frame(
        structure_smiles_2D = new$structure_smiles_2D[[X]],
        query = new$query[[X]],
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
        structure_smiles_2D = new$structure_smiles_2D[[X]],
        query = new$query[[X]],
        pathway = NA,
        superclass = NA,
        class = NA,
        glycoside = NA
      )
      return(df)
    }
  )
}
