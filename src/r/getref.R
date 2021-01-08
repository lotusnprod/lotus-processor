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
getref <- function(X) {
  tryCatch(
    {
      cr_works(
        flq = c(query.bibliographic = X),
        sort = "score",
        order = "desc",
        limit = 1
      )
    },
    error = function(e) {
      NA
    }
  )
}
