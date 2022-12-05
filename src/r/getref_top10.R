source("paths.R")
library(future)
library(future.apply)
library(progressr)
library(rcrossref)

source("r/progressr.R")

#' Title
#'
#' @param xs
#'
#' @return
#' @export
#'
#' @examples
getref_top10 <- function(xs) {
  p <- progressr::progressor(along = xs)
  future.apply::future_lapply(
    future.seed = TRUE,
    X = xs,
    FUN = function(x) {
      tryCatch(
        {
          p(sprintf("x=%g", as.numeric(x))) ## little hack
          rcrossref::cr_works(
            flq = c(query.bibliographic = x),
            sort = "score",
            order = "desc",
            limit = 10
          )
        },
        error = function(e) {
          NA
        }
      )
    }
  )
}
