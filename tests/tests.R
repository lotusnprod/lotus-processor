# setwd(dir = "~/gitlab/opennaturalproductsdb/src")

library(testthat)
library(tidyverse)
source("paths.R")
source("r/log.R")

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

log_debug("loading ...")
log_debug("... initial table")
originalTable <- read_delim(
  file = gzfile(description = pathDataInterimTablesOriginalTable),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  tibble()

log_debug("... results ...")
log_debug("... organisms")
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
  ) %>%
  tibble()

log_debug("... structures ...")
log_debug("... translated")
translatedStructureTable <- read_delim(
  file = gzfile(description = pathDataInterimTablesTranslatedStructureFinal),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  tibble()

log_debug("... cleaned")
cleanedStructureTableFull <- read_delim(
  file = gzfile(description = pathDataInterimTablesCleanedStructureNamed),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  select(
    structureTranslated,
    structureCleanedSmiles = smilesSanitized,
    structureCleanedInchi = inchiSanitized,
    structureCleanedInchikey3D = inchikeySanitized,
    structureCleaned_inchikey2D = shortikSanitized,
    structureCleaned_molecularFormula = formulaSanitized,
    structureCleaned_exactMass = exactmassSanitized,
    structureCleaned_xlogp = xlogpSanitized,
    structureCleaned_stereocenters_unspecified = count_unspecified_atomic_stereocenters,
    structureCleaned_stereocenters_total = count_atomic_stereocenters,
    structureCleaned_nameIupac,
    structureCleaned_nameTraditional
  ) %>%
  tibble()

## no classification (classyFire, NP-classifier) test at the moment

log_debug("... references")
referenceTableFull <- read_delim(
  file = gzfile(description = pathDataInterimTablesCleanedReferenceFile),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate(referenceCleanedDoi = toupper(referenceCleanedDoi)) %>%
  select(-referenceCleaned_score_crossref) %>%
  tibble()

log_debug("joining structures")
structureFull <-
  left_join(translatedStructureTable, cleanedStructureTableFull) %>%
  select(-structureTranslated) %>%
  tibble()

log_debug("... expectations ...")
log_debug("... organisms")
organismTableFullExpectation <- read_delim(
  file = pathTestsOrganisms,
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  tibble()

log_debug("... structures")
structureFullExpectation <- read_delim(
  file = pathTestsStructures,
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  select(-structureCleaned_validatorLog) %>%
  tibble()

log_debug("... references")
referenceTableFullExpectation <- read_delim(
  file = pathTestsReferences,
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  select(-referenceCleaned_score_crossref) %>%
  tibble()

log_debug("testing ...")
log_debug("... organisms")
test_that(
  desc = "organisms",
  code = expect_equal(
    object = organismTableFull,
    expected = organismTableFullExpectation
  )
)

log_debug("... structures")
test_that(
  desc = "structures",
  code = expect_equal(
    object = structureFull,
    expected = structureFullExpectation
  )
)

log_debug("... references")
test_that(
  desc = "references",
  code = expect_equal(
    object = referenceTableFull,
    expected = referenceTableFullExpectation
  )
)

# write.table(
#   x = organismTableFull,
#   file = pathTestsOrganisms,
#   row.names = FALSE,
#   quote = FALSE,
#   sep = "\t",
#   fileEncoding = "UTF-8"
# )
#
# write.table(
#   x = structureFull,
#   file = pathTestsStructures,
#   row.names = FALSE,
#   quote = FALSE,
#   sep = "\t",
#   fileEncoding = "UTF-8"
# )
#
# write.table(
#   x = referenceTableFull,
#   file = pathTestsReferencces,
#   row.names = FALSE,
#   quote = FALSE,
#   sep = "\t",
#   fileEncoding = "UTF-8"
# )


test_that(desc = "taxon", {
  expect_success(expect_taxon(
    organismTableFull,
    c("Iris"),
    c("Plantae", "Viridiplantae", NA)
  ))
  expect_failure(expect_length(1, 2), "has length 1, not length 2.")
  expect_success(expect_length(1:10, 10))
  expect_success(expect_length(letters[1:5], 5))
})
