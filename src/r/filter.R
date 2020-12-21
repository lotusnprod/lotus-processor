################################   Functions   ################################

#' Title
#'
#' @param dataframe
#'
#' @return
#' @export
#'
#' @examples
filter_dirty <- function(dataframe) {
  dfDoi <- dataframe %>%
    filter(referenceType == "doi") %>%
    filter(database != "coc_1" |
      as.numeric(referenceCleaned_score_titleOrganism) == 1)

  dfOriginal <- dataframe %>%
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

  dfPublishingDetails <- dataframe %>%
    filter(referenceType == "publishingDetails") %>%
    filter(as.numeric(referenceCleaned_score_titleOrganism) == 1)

  dfPubmed <- dataframe %>%
    filter(referenceType == "pubmed") %>%
    filter(as.numeric(referenceCleaned_score_titleOrganism) == 1)

  dfSplit <- dataframe %>%
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

  dfTitle <- dataframe %>%
    filter(referenceType == "title") %>%
    filter(as.numeric(referenceCleaned_score_distance) <= 10)

  cleanDataframe <-
    bind_rows(
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
