require(stringr)
source("r/y_as_na.R")

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
preparing_name <- function(x) {
  x$nameCleaned <- x$structureOriginal_nominal
  x$nameCleaned <-
    gsub(
      pattern = "\\u03b1",
      replacement = "alpha",
      x = x$nameCleaned,
      fixed = TRUE
    )
  x$nameCleaned <-
    gsub(
      "\\u03b2",
      replacement = "beta",
      x = x$nameCleaned,
      fixed = TRUE
    )
  x$nameCleaned <-
    gsub(
      "\\u03b3",
      replacement = "gamma",
      x = x$nameCleaned,
      fixed = TRUE
    )
  x$nameCleaned <-
    gsub(
      "\\u03b4",
      replacement = "delta",
      x = x$nameCleaned,
      fixed = TRUE
    )
  x$nameCleaned <-
    gsub(
      pattern = "\\u03b5",
      replacement = "epsilon",
      x = x$nameCleaned,
      fixed = TRUE
    )
  x$nameCleaned <-
    gsub(
      pattern = "\\u03c9",
      replacement = "omega",
      x = x$nameCleaned,
      fixed = TRUE
    )
  x$nameCleaned <-
    gsub(
      pattern = "(\\u00b1)-",
      replacement = "",
      x = x$nameCleaned,
      fixed = TRUE
    )
  x$nameCleaned <-
    gsub(
      pattern = "\\u00f6",
      replacement = "ö",
      x = x$nameCleaned,
      fixed = TRUE
    )
  x$nameCleaned <-
    gsub(
      pattern = "\\u2192",
      replacement = "->",
      x = x$nameCleaned,
      fixed = TRUE
    )
  x$nameCleaned <-
    gsub(
      pattern = "α",
      replacement = "alpha",
      x = x$nameCleaned
    )
  x$nameCleaned <-
    gsub(
      pattern = "Α",
      replacement = "alpha",
      x = x$nameCleaned
    )
  x$nameCleaned <-
    gsub(
      pattern = "β",
      replacement = "beta",
      x = x$nameCleaned
    )
  x$nameCleaned <-
    gsub(
      pattern = "Β",
      replacement = "beta",
      x = x$nameCleaned
    )
  x$nameCleaned <-
    gsub(
      pattern = "γ",
      replacement = "gamma",
      x = x$nameCleaned
    )
  x$nameCleaned <-
    gsub(
      pattern = "Γ",
      replacement = "gamma",
      x = x$nameCleaned
    )
  x$nameCleaned <-
    gsub(
      pattern = "δ",
      replacement = "delta",
      x = x$nameCleaned
    )
  x$nameCleaned <-
    gsub(
      pattern = "Δ",
      replacement = "delta",
      x = x$nameCleaned
    )
  x$nameCleaned <-
    gsub(
      pattern = "ε",
      replacement = "epsilon",
      x = x$nameCleaned
    )
  x$nameCleaned <-
    gsub(
      pattern = "Ε",
      replacement = "epsilon",
      x = x$nameCleaned
    )
  x$nameCleaned <-
    gsub(
      pattern = "- ",
      replacement = "-",
      x = x$nameCleaned
    )
  x$nameCleaned <-
    gsub(
      pattern = "–",
      replacement = "-",
      x = x$nameCleaned
    )
  x$nameCleaned <-
    gsub(
      pattern = "\\) ",
      replacement = "\\)",
      x = x$nameCleaned
    )
  x$nameCleaned <-
    gsub(
      pattern = "\"\"\"",
      replacement = "",
      x = x$nameCleaned
    )

  x$nameCleaned <- trimws(x = x$nameCleaned)
  x$nameCleaned <- tolower(x = x$nameCleaned)
  if (mode != "custom") {
    x$nameCleaned <-
      gsub(
        pattern = "-NA$",
        replacement = "",
        x = gsub(
          pattern = ",$",
          replacement = "",
          x = (paste0(
            gsub(
              pattern = ".*,([^.]+)\\:.*",
              replacement = "\\1",
              x = x$nameCleaned
            ),
            "-",
            str_extract(
              pattern = ".*,",
              string = x$nameCleaned
            )
          ))
        )
      )
  }
  x$nameCleaned <- y_as_na(x = x$nameCleaned, y = "NA")

  return(x)
}
