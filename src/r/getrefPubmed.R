source("paths.R")
library(rentrez)

#' Title
#'
#' @param X
#'
#' @return
#' @export
#'
#' @examples
getrefPubmed <- function(X) {
  tryCatch(
    {
      df <- entrez_summary(db = "pubmed", id = X)

      translatedDoi <-
        ifelse(test = "doi" %in% df[["articleids"]][, 1],
          yes = trimws(df[["articleids"]][["value"]][[which(df[["articleids"]] == "doi")]]),
          no = NA
        )
      translatedJournal <- df[["fulljournalname"]]
      translatedTitle <- df[["title"]]
      translatedAuthor <- df[["sortfirstauthor"]]
      translatedDate <- df[["pubdate"]]

      ids <-
        data.frame(
          translatedDoi,
          translatedJournal,
          translatedTitle,
          translatedAuthor,
          translatedDate
        )
      return(ids)
    },
    error = function(e) {
      data.frame(
        "translatedDoi" = NA,
        "translatedJournal" = NA,
        "translatedTitle" = NA,
        "translatedAuthor" = NA,
        "translatedDate" = NA
      )
    }
  )
}
