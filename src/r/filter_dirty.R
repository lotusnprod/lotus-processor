#' Title
#'
#' @param dataframe
#'
#' @return
#' @export
#'
#' @examples
filter_dirty <- function(dataframe) {
  dfWiki <- dataframe %>%
    filter(database == "wikidata")

  dfRest <- dataframe %>%
    filter(database != "wikidata")

  dfDoi <- dfRest %>%
    filter(referenceType == "doi") %>%
    filter(
      database != "coconut" |
        as.numeric(referenceCleaned_score_titleOrganism) == 1
    )

  dfOriginal <- dfRest %>%
    filter(referenceType == "original") %>%
    filter(
      (as.numeric(referenceCleaned_score_complementTotal) >= 2 |
        grepl(
          pattern = "Phytochemistry",
          x = referenceCleaned_journal,
          ignore.case = TRUE
        ) |
        grepl(
          pattern = "Journal of Agricultural and Food Chemistry",
          x = referenceCleaned_journal,
          ignore.case = TRUE
        ) |
        grepl(
          pattern = "Journal of Natural Products",
          x = referenceCleaned_journal,
          ignore.case = TRUE
        )) &
        as.numeric(referenceCleaned_score_titleOrganism) == 1
    )

  dfPublishingDetails <- dfRest %>%
    filter(referenceType == "publishingDetails") %>%
    filter(as.numeric(referenceCleaned_score_titleOrganism) == 1)

  dfPubmed <- dfRest %>%
    filter(referenceType == "pubmed") %>%
    filter(as.numeric(referenceCleaned_score_titleOrganism) == 1)

  dfSplit <- dfRest %>%
    filter(referenceType == "split") %>%
    filter(
      (as.numeric(referenceCleaned_score_complementTotal) >= 2 |
        grepl(
          pattern = "Phytochemistry",
          x = referenceCleaned_journal,
          ignore.case = TRUE
        ) |
        grepl(
          pattern = "Journal of Agricultural and Food Chemistry",
          x = referenceCleaned_journal,
          ignore.case = TRUE
        ) |
        grepl(
          pattern = "Journal of Natural Products",
          x = referenceCleaned_journal,
          ignore.case = TRUE
        )) &
        as.numeric(referenceCleaned_score_titleOrganism) == 1
    )

  #' additional coconut cleaning
  if (nrow(dfSplit) != 0) {
    dfSplit_1 <- dfSplit %>%
      filter(database == "coconut")

    dfSplit_1_1 <- dfSplit %>%
      filter(database == "coconut") %>%
      rowwise() %>%
      filter(grepl(pattern = organismValue, x = referenceCleanedTitle))

    dfSplit_1_2 <- dfSplit %>%
      filter(database == "coconut") %>%
      anti_join(
        dfSplit_1_1,
        by = c(
          "structureValue" = "structureValue",
          "referenceValue" = "referenceValue"
        )
      )

    dfSplit_2 <- bind_rows(dfSplit_1_1, dfSplit_1_2)

    dfSplit_3 <- dfSplit %>%
      filter(database != "coconut")
  } else {
    dfSplit_2 <- dfSplit
    dfSplit_3 <- dfSplit
  }

  dfTitle <- dfRest %>%
    filter(referenceType == "title") %>%
    filter(as.numeric(referenceCleaned_score_distance) <= 10)

  cleanDataframe <-
    bind_rows(
      dfWiki,
      dfDoi,
      dfPubmed,
      dfTitle,
      dfOriginal,
      dfSplit_2,
      dfSplit_3,
      dfPublishingDetails
    )

  return(cleanDataframe)
}

###############################################################################
