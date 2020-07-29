#######################################################
####################   Functions   ####################
#######################################################

library(Hmisc)
library(data.table)
library(dplyr)
library(parallel)
library(pbmcapply)
library(rcrossref)
library(readr)
library(rentrez)
library(splitstackshape)
library(stringdist)
library(stringr)
library(stringi)
library(tidyr)

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
  
  names(reflist) <- 1:length(reflist)
  
  bound <- bind_rows(reflist[!is.na(reflist)], .id = "level")
  
  joined <- left_join(tableIntial, bound)
  
  joinedData <- joined$data
  
  joinedSelected <- joined %>%
    select(all_of(referenceColumnName),
           level)
  
  tableInterim <- cbind(joinedSelected, joinedData)
  
  for (i in 1:nrow(tableInterim)) {
    tableInterim[i, "distScore"] <-
      stringdist(a = tableInterim[i, 1],
                 b = tableInterim[i, "title"],
                 method = method)
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
  
  names(reflist) <- 1:length(reflist)
  
  bound <- bind_rows(reflist[!is.na(reflist)], .id = "level")
  
  joined <- left_join(tableIntial, bound)
  
  joinedData <- joined$data
  
  joinedSelected <- joined %>%
    select(all_of(referenceColumnName),
           level)
  
  tableInterim <- cbind(joinedSelected, joinedData)
  
  for (i in 1:nrow(tableInterim)) {
    tableInterim[i, "distScore"] <-
      stringdist(a = tableInterim[i, 1],
                 b = tableInterim[i, "title"],
                 method = method)
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
