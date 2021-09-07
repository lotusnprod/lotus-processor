## courtesy of the Hmisc package, jusst not loading the whole package
## and its long dependencies for it

#' Title
#'
#' @param string
#'
#' @return
#' @export
#'
#' @examples
capitalize <- function(string) {
  capped <- grep("^[A-Z]", string, invert = TRUE)
  substr(string[capped], 1, 1) <-
    toupper(substr(string[capped], 1, 1))
  return(string)
}
