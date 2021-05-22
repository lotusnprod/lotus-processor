source("paths.R")
library(jsonlite)
library(tidyverse)
library(rvest)

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
        pathway = ifelse(
          test = is_empty(result$pathway_results),
          yes = NA,
          no = result$pathway_results
        ),
        superclass = ifelse(
          test = is_empty(result$superclass_results),
          yes = NA,
          no = result$superclass_results
        ),
        class = ifelse(
          test = is_empty(result$class_results),
          yes = NA,
          no = result$class_results
        ),
        glycoside = ifelse(
          test = is_empty(result$isglycoside),
          yes = NA,
          no = result$isglycoside
        )
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
