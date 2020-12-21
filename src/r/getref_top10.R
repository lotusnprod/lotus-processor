library(rcrossref)

#' Title
#'
#' @param X
#'
#' @return
#' @export
#'
#' @examples
getref_top10 <- function(X) {
  tryCatch(
    {
      cr_works(
        flq = c(query.bibliographic = X),
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
