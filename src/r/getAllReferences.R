source("paths.R")
library(groundhog)
groundhog.library(stringdist, date = groundhog.day)

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
getAllReferences <- function(data, referenceType, method = "osa") {
  referenceColumnName <- paste("referenceOriginal",
    referenceType,
    sep = "_"
  )

  tableIntial <- data %>%
    mutate(level = row.names(.))

  dataList <- list()
  for (i in 1:length(reflist)) {
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

  names(dataList) <- 1:length(dataList)

  bound <- bind_rows(dataList[!is.na(dataList)], .id = "level")

  tableInterim <- left_join(tableIntial, bound)

  for (i in 1:nrow(tableInterim)) {
    tableInterim[i, "distScore"] <-
      ifelse(
        test = "title" %in% colnames(tableInterim),
        yes =
          stringdist(
            a = as.character(tolower(tableInterim[i, 1])),
            # method is case sensitive
            b = as.character(tolower(tableInterim[i, "title"])),
            # method is case sensitive
            method = method
          ),
        no = NA
      )
  }


  if (!"title" %in% colnames(tableInterim)) {
    tableFinal <- tableInterim %>%
      mutate(
        referenceTranslatedDoi = NA,
        referenceTranslatedJournal = NA,
        referenceTranslatedTitle = NA,
        referenceTranslatedDate = NA,
        referenceTranslatedAuthor = NA,
        referenceTranslationScoreCrossref = NA,
        referenceTranslationScoreDistance = NA,
      ) %>%
      select(
        all_of(referenceColumnName),
        referenceTranslatedDoi,
        referenceTranslatedJournal,
        referenceTranslatedTitle,
        referenceTranslatedDate,
        referenceTranslatedAuthor,
        referenceTranslationScoreCrossref,
        referenceTranslationScoreDistance
      ) %>%
      mutate_all(as.character)
  }

  if ("title" %in% colnames(tableInterim)) {
    tableFinal <- tableInterim %>%
      select(
        all_of(referenceColumnName),
        referenceTranslatedDoi = doi,
        referenceTranslatedJournal = container.title,
        referenceTranslatedTitle = title,
        referenceTranslatedDate = issued,
        referenceTranslatedAuthor = author,
        referenceTranslationScoreCrossref = score,
        referenceTranslationScoreDistance = distScore,
      ) %>%
      unnest(
        cols = c(referenceTranslatedAuthor),
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
      mutate_all(as.character)
  }

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
