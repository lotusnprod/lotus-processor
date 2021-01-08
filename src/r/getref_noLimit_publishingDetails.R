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
getref_noLimit_publishingDetails <- function(X) {
  tryCatch(
    {
      cr_works(
        query = X,
        sort = "score",
        order = "desc"
      )
    },
    error = function(e) {
      NA
    }
  )
}
