#######################################################
####################   Functions   ####################
#######################################################

library(data.table)
library(parallel)
library(pbmcapply)
library(rcrossref)
library(rentrez)
library(splitstackshape)
library(stringdist)
library(stringi)
library(tidyverse)

#######################################################
#######################################################

getref <- function(X) {
  tryCatch({
    cr_works(
      flq = c(query.bibliographic = X),
      sort = 'score',
      order = "desc",
      limit = 1
    )
  },
  error = function(e) {
    NA
  })
}

#######################################################
#######################################################

getref_noLimit <- function(X) {
  tryCatch({
    cr_works(
      flq = c(query.bibliographic = X),
      sort = 'score',
      order = "desc"
    )
  },
  error = function(e) {
    NA
  })
}

#######################################################
#######################################################

getrefPubmed <- function(X) {
  tryCatch({
    df <- entrez_summary(db = "pubmed", id = X)
    
    translatedDoi <-
      ifelse(test = "doi" %in% df[["articleids"]][, 1],
             yes = trimws(df[["articleids"]][["value"]][[which(df[["articleids"]] == "doi")]]),
             no = NA)
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
  })
}

#######################################################
#######################################################

getrefDoi <- function(X) {
  tryCatch({
    cr_works(dois = X)
  },
  error = function(e) {
    NA
  })
}

#######################################################
#######################################################

getBestReference <- function(data, referenceType, method = "osa") {
  referenceColumnName <- paste("referenceOriginal",
                               referenceType,
                               sep = "_")
  
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
      stringdist(
        a = as.character(tableInterim[i, 1]),
        b = as.character(tableInterim[i, "title"]),
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
    unnest(cols = c(referenceTranslatedAuthor),
           keep_empty = TRUE) %>%
    distinct_at(., .vars = (c(1,
                              4)),
                .keep_all = TRUE) %>%
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
    distinct_at(., .vars = (c(1)),
                .keep_all = TRUE) %>%
    ungroup() %>%
    mutate_all(as.character)
  
  tableFinal[] <-
    lapply(tableFinal, function(x)
      gsub("\r\n", " ", x))
  tableFinal[] <-
    lapply(tableFinal, function(x)
      gsub("\r", " ", x))
  tableFinal[] <-
    lapply(tableFinal, function(x)
      gsub("\n", " ", x))
  tableFinal[] <-
    lapply(tableFinal, function(x)
      gsub("\t", " ", x))
  
  return(tableFinal)
}

#######################################################
#######################################################

getAllReferences <- function(data, referenceType, method = "osa") {
  referenceColumnName <- paste("referenceOriginal",
                               referenceType,
                               sep = "_")
  
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
            a = as.character(tableInterim[i, 1]),
            b = as.character(tableInterim[i, "title"]),
            method = method
          ),
        no = NA
      )
  }
  
  
  if (!"title" %in% colnames(tableInterim))
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
  
  if ("title" %in% colnames(tableInterim))
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
    unnest(cols = c(referenceTranslatedAuthor),
           keep_empty = TRUE) %>%
    distinct_at(., .vars = (c(1,
                              4)),
                .keep_all = TRUE) %>%
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
  
  tableFinal[] <-
    lapply(tableFinal, function(x)
      gsub("\r\n", " ", x))
  tableFinal[] <-
    lapply(tableFinal, function(x)
      gsub("\r", " ", x))
  tableFinal[] <-
    lapply(tableFinal, function(x)
      gsub("\n", " ", x))
  tableFinal[] <-
    lapply(tableFinal, function(x)
      gsub("\t", " ", x))
  
  return(tableFinal)
}
#######################################################
#######################################################

getref_top10 <- function(X) {
  tryCatch({
    cr_works(
      flq = c(query.bibliographic = X),
      sort = 'score',
      order = "desc",
      limit = 10
    )
  },
  error = function(e) {
    NA
  })
}

#######################################################
#######################################################