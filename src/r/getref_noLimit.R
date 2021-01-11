source("paths.R")
library(rcrossref)

#' Title
#'
#' @param X
#'
#' @return
#' @export
#'
#' @examples
getref_noLimit <- function(X) {
  tryCatch(
    {
      cr_works(
        flq = c(query.bibliographic = X),
        sort = "score",
        order = "desc"
      )
    },
    error = function(e) {
      NA
    }
  )
}
