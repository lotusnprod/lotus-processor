source("paths.R")
library(vroom)

#' Title
#'
#' @param path
#'
#' @return
#' @export
#'
#' @examples
vroom_read_safe <- function(path) {
  vroom(
    file = path,
    delim = "\t",
    escape_double = TRUE,
    col_types = cols(.default = "c"),
    num_threads = 1
  )
}

#' Title
#'
#' @param x
#' @param path
#'
#' @return
#' @export
#'
#' @examples
vroom_write_safe <- function(x, path) {
  vroom_write(
    x = x,
    path = gzfile(
      description = path,
      compression = 9,
      encoding = "UTF-8"
    ),
    num_threads = 1,
    bom = TRUE,
    quote = "all",
    escape = "double",
    delim = "\t",
    col_names = TRUE,
    progress = TRUE,
    append = FALSE
  )
}
