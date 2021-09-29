source("paths.R")
library(stringdist)
library(tidyr)

#' Title
#'
#' @param data
#' @param referenceType
#' @param method
#'
#' @return
#' @export
#'
#' @examples
getBestReference <- function(data, referenceType, method = "osa") {
  referenceColumnName <- paste("referenceOriginal",
    referenceType,
    sep = "_"
  )

  tableIntial <- data %>%
    mutate(level = row.names(.))

  dataList <- list()
  for (i in seq_along(reflist)) {
    ifelse(
      test = !is.na(reflist[[i]]),
      yes = ifelse(
        test = !is.null(reflist[[i]][["data"]]),
        yes = dataList[[i]] <- data.frame(reflist[[i]][["data"]]),
        no = dataList[[i]] <- data.frame(NA)
      ),
      no = dataList[[i]] <- data.frame(NA)
    )
  }

  names(dataList) <- seq_along(dataList)

  bound <- bind_rows(dataList[!is.na(dataList)], .id = "level")

  tableInterim <- left_join(tableIntial, bound)

  for (i in seq_len(nrow(tableInterim))) {
    tableInterim[i, "distScore"] <-
      stringdist(
        a = as.character(tolower(tableInterim[i, 1])),
        # method is case sensitive
        b = as.character(tolower(tableInterim[i, "title"])),
        # method is case sensitive
        method = method
      )
  }

  tableFinal <- tableInterim %>%
    select(
      all_of(referenceColumnName),
      referenceTranslatedDoi = doi,
      referenceTranslatedJournal = container.title,
      referenceTranslatedTitle = title,
      referenceTranslatedDate = issued,
      referenceTranslatedAuthor = author,
      referenceTranslationScoreCrossref = score,
      referenceTranslationScoreDistance = distScore
    ) %>%
    unnest(
      cols = referenceTranslatedAuthor,
      keep_empty = TRUE
    ) %>%
    distinct_at(.,
      .vars = (c(
        1,
        4
      )),
      .keep_all = TRUE
    ) %>%
    select(
      all_of(referenceColumnName),
      referenceTranslatedDoi,
      referenceTranslatedJournal,
      referenceTranslatedTitle,
      referenceTranslatedDate,
      referenceTranslatedAuthor = family,
      referenceTranslationScoreCrossref,
      referenceTranslationScoreDistance
    ) %>%
    group_by_at(1) %>%
    arrange(referenceTranslationScoreDistance) %>%
    distinct_at(.,
      .vars = (1),
      .keep_all = TRUE
    ) %>%
    ungroup() %>%
    mutate_all(as.character)

  tableFinal[] <-
    lapply(tableFinal, function(x) {
      gsub("\r\n", " ", x)
    })
  tableFinal[] <-
    lapply(tableFinal, function(x) {
      gsub("\r", " ", x)
    })
  tableFinal[] <-
    lapply(tableFinal, function(x) {
      gsub("\n", " ", x)
    })
  tableFinal[] <-
    lapply(tableFinal, function(x) {
      gsub("\t", " ", x)
    })

  return(tableFinal)
}
