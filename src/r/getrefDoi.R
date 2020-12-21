library(rcrossref)

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
