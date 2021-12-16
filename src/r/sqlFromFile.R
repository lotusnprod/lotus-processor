# taken from https://stackoverflow.com/questions/18914283/how-to-execute-more-than-one-rsqlite-statement-at-once-or-how-to-dump-a-whole-fi

#' Title
#'
#' @param file
#'
#' @return
#' @export
#'
#' @examples
sqlFromFile <- function(file) {
  require(stringr)
  sql <- suppressWarnings(readLines(file))
  sql <- unlist(str_split(paste(sql, collapse = " "), ";"))
  sql <- sql[grep("^ *$", sql, invert = TRUE)]
  sql
}
