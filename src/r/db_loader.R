#' Title
#'
#' @param path_to_db
#'
#' @return
#' @export
#'
#' @examples
db_loader <- function(path_to_db) {
  db <- read_delim(
    file = gzfile(path_to_db),
    col_types = cols(.default = "c"),
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  )
  return(db)()
}