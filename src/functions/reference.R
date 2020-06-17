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
