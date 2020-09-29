###############################################################################
################################   Functions   ################################
###############################################################################

filter_dirty <- function(dataframe) {
  cleanDataframe <- dataframe %>%
    filter(
      as.numeric(referenceCleaned_score_crossref) == 1 |
        as.numeric(referenceCleaned_score_distance) <= 10 |
        as.numeric(referenceCleaned_score_complementTotal) >= 2 |
        as.numeric(referenceCleaned_score_titleOrganism) == 1 |
        referenceCleaned_journal == "Phytochemistry" |
        referenceCleaned_journal == "Journal of Agricultural and Food Chemistry" |
        referenceCleaned_journal == "Journal of Natural Products"
    ) %>%
    filter(!(
      referenceType == "doi" &
        database == "coc_1" &
        as.numeric(referenceCleaned_score_titleOrganism) != 1
    )) %>%
    filter(!(
      referenceType == "original" &
        as.numeric(referenceCleaned_score_titleOrganism) != 1
    ))  %>%
    filter(!(
      referenceType == "publishingDetails" &
        as.numeric(referenceCleaned_score_titleOrganism) != 1
    )) %>%
    filter(!(
      referenceType == "pubmed" &
        as.numeric(referenceCleaned_score_titleOrganism) != 1
    )) %>%
    filter(!(referenceType == "split" & (
      as.numeric(referenceCleaned_score_titleOrganism) != 1 |
        as.numeric(referenceCleaned_score_complementTotal) < 2
    ))) %>%
    filter(!(
      referenceType == "title" &
        as.numeric(referenceCleaned_score_distance) > 10
    ))
  
  # %>%
  #   filter((((referenceType == "title" |
  #               referenceType == "original") &
  #              (
  #                referenceCleaned_journal == "Phytochemistry" |
  #                  referenceCleaned_journal == "Journal of Agricultural and Food Chemistry" |
  #                  referenceCleaned_journal == "Journal of Natural Products" |
  #                  as.numeric(referenceCleaned_score_distance) <= 10
  #              )
  #   ) |
  #     (
  #       referenceType == "doi" & database != "coc_1"
  #     )) )
  
  return(cleanDataframe)
}

###############################################################################
###############################################################################
