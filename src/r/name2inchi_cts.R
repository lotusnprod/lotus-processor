source("paths.R")
library(webchem)

#' Title
#'
#' @param i
#'
#' @return
#' @export
#'
#' @examples
name2inchi_cts <- function(i) {
  tryCatch(
    {
      x <- cts_convert(
        query = dataForCTS[i, "nameCleaned"],
        from = "Chemical Name",
        to = "InChI Code",
        verbose = FALSE,
        match = "first"
      )[[1]]
      return(x)
    },
    error = function(e) {
      NA
    }
  )
}
