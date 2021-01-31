# taken from https://stackoverflow.com/questions/18914283/how-to-execute-more-than-one-rsqlite-statement-at-once-or-how-to-dump-a-whole-fi

#' Title
#'
#' @param conn
#' @param sql
#'
#' @return
#' @export
#'
#' @examples
dbSendQueries <- function(conn, sql) {
  dummyfunction <- function(sql, conn) {
    dbSendQuery(conn, sql)
  }
  lapply(sql, dummyfunction, conn)
}
