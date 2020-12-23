#' Title
#'
#' @param path_to_db
#'
#' @return
#' @export
#'
#' @examples
db_loader <- function(path_to_db) {
  db <- vroom(
    file = gzfile(path_to_db),
    col_types = cols(.default = "c"),
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE,
    quote = "",
    num_threads = 1
  )
  return(db)()
}