source("paths.R")
library(groundhog)
groundhog.library(rcrossref, date = groundhog.day)

#' Title
#'
#' @param X
#'
#' @return
#' @export
#'
#' @examples
getrefDoi <- function(X) {
  tryCatch(
    {
      cr_works(dois = X)
    },
    error = function(e) {
      NA
    }
  )
}
