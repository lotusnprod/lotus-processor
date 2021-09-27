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
    filter(database != "coconut" |
      as.numeric(referenceCleaned_score_titleOrganism) == 1)

  dfOriginal <- dfRest %>%
    filter(referenceType == "original") %>%
    filter(
      (
        as.numeric(referenceCleaned_score_complementTotal) >= 2 |
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
          )
      ) &
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
      (
        as.numeric(referenceCleaned_score_complementTotal) >= 2 |
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
          )
      ) &
        as.numeric(referenceCleaned_score_titleOrganism) == 1
    )

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
      dfSplit,
      dfPublishingDetails
    )

  return(cleanDataframe)
}

###############################################################################
