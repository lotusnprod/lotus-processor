source("paths.R")
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
name2inchi_cactus <- function(xs) {
  p <- progressr::progressor(along = xs)
  future.apply::future_lapply(
    future.seed = TRUE,
    X = xs,
    FUN = function(x) {
      tryCatch(
        {
          p(sprintf("x=%g", x))
          cpd <- dataForCactus[x, "nameCleaned"]
          url <-
            paste0(
              "https://cactus.nci.nih.gov/chemical/structure/",
              cpd,
              "/stdinchi"
            )
          url <- gsub(
            pattern = "\\s",
            replacement = "%20",
            x = url
          )
          rvest::read_html(url) |>
            rvest::html_text()
        },
        error = function(e) {
          NA
        }
      )
    }
  )
}

#' Title
#'
#' @param xs
#'
#' @return
#' @export
#'
#' @examples
name2smiles_cactus <- function(xs) {
  p <- progressr::progressor(along = xs)
  future.apply::future_lapply(
    future.seed = TRUE,
    X = xs,
    FUN = function(x) {
      tryCatch(
        {
          p(sprintf("x=%g", x))
          cpd <- dataForCactus[x, "nameCleaned"]
          url <-
            paste0(
              "https://cactus.nci.nih.gov/chemical/structure/",
              cpd,
              "/smiles"
            )
          url <- gsub(
            pattern = "\\s",
            replacement = "%20",
            x = url
          )
          rvest::read_html(url) |>
            rvest::html_text()
        },
        error = function(e) {
          NA
        }
      )
    }
  )
}
