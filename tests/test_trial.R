library(testthat)
library(tidyverse)
setwd("~/GitLab/opennaturalproductsdb/src")
source("paths.R")

#' Check if taxon is attributed to the right kingdom
#'
#' @param table_organisms The table contaning the taxon/taxa you want to test
#' @param taxon_detected The taxon/taxa you want to test
#' @param kingdom_cleaned The kingdom(s) you want them to belong to
#'
#' @return If expectation is met or not
#' @export
#'
#' @examples
#' expect_taxon(organismTableFull, c("Iris"), c("Plantae", "Viridiplantae", NA))
expect_taxon <-
  function(table_organisms,
           taxon_detected,
           kingdom_cleaned) {
    act <- quasi_label(enquo(taxon_detected), arg = "object")
    act$kingdom <- table_organisms %>%
      filter(word(string = taxon_detected, end = 1) == organismDetected) %>%
      distinct(organismDetected, organismCleaned_dbTaxo_1kingdom)
    expect(
      ok = all(
        sort(act$kingdom$organismCleaned_dbTaxo_1kingdom) %in% sort(kingdom_cleaned)
      ),
      failure_message = sprintf(
        "%s belongs to %s, not  %s.",
        act$lab,
        sort(act$kingdom$organismCleaned_dbTaxo_1kingdom),
        sort(kingdom_cleaned)
      )
    )
    invisible(act$val)
  }

organismTableFull <- read_delim(
  file = gzfile(description = pathDataInterimTablesCleanedOrganismFinal),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  select(
    organismOriginal,
    organismDetected,
    organismCleaned,
    organismCleaned_dbTaxo = organismDbTaxo,
    organismCleaned_dbTaxoTaxonIds = organismTaxonIds,
    organismCleaned_dbTaxoTaxonRanks = organismTaxonRanks,
    organismCleaned_dbTaxoTaxonomy = organismTaxonomy,
    organismCleaned_dbTaxo_1kingdom = organism_1_kingdom,
    organismCleaned_dbTaxo_2phylum = organism_2_phylum,
    organismCleaned_dbTaxo_3class = organism_3_class,
    organismCleaned_dbTaxo_4order = organism_4_order,
    organismCleaned_dbTaxo_5family = organism_5_family,
    organismCleaned_dbTaxo_6genus = organism_6_genus,
    organismCleaned_dbTaxo_7species = organism_7_species,
    organismCleaned_dbTaxo_8variety = organism_8_variety
  )

test_that("taxon is correct", {
  expect_success(expect_taxon(
    organismTableFull,
    c("Iris"),
    c("Plantae", "Viridiplantae", NA)
  ))
  expect_failure(expect_length(1, 2), "has length 1, not length 2.")
  expect_success(expect_length(1:10, 10))
  expect_success(expect_length(letters[1:5], 5))
})
