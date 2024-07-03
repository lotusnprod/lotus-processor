# setwd(dir = "~/git/lotus-processor/src")

mode_test <- TRUE

library(dplyr)
library(readr)
library(testthat)
source("paths.R")
source("r/log_debug.R")

log_debug("loading ...")
log_debug("... initial table")
originalTable <- read_delim(
  file = gzfile(description = pathDataInterimTablesOriginalTable),
  delim = "\t",
  col_types = cols(.default = "c"),
  locale = locales,
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  tibble()

log_debug("... results ...")
log_debug("... organisms")
organismTableFull <- read_delim(
  file = gzfile(description = pathDataInterimTablesProcessedOrganismFinal),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
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
  file = gzfile(description = pathDataInterimTablesProcessedStructureNamed),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  select(-validatorLog) %>%
  tibble()

## no classification (classyFire, NP-classifier) test at the moment

log_debug("... references")
referenceTableFull <- read_delim(
  file = gzfile(description = pathDataInterimTablesProcessedReferenceFile),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  select(-referenceCleaned_score_crossref) %>%
  tibble()

log_debug("validated table ")
validatedTable <- read_delim(
  file = gzfile(description = pathDataInterimTablesAnalyzedPlatinum),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
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
  tibble()

log_debug("... references")
referenceTableFullExpectation <- read_delim(
  file = pathTestsReferences,
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  tibble()

log_debug("... validated table")
validatedTableExpectation <- read_delim(
  file = pathTestsPlatinum,
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
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

log_debug("... not testing unvalidated references since they vary because of Crossref")
# test_that(
#   desc = "references",
#   code = expect_equal(object = referenceTableFull,
#                       expected = referenceTableFullExpectation)
# )

# log_debug("... validated referenced pairs")
# test_that(
#   desc = "validated referenced pairs",
#   code = expect_equal(object = validatedTable,
#                       expected = validatedTableExpectation)
# )

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
#   file = pathTestsReferences,
#   row.names = FALSE,
#   quote = FALSE,
#   sep = "\t",
#   fileEncoding = "UTF-8"
# )
#
# write.table(
#   x = validatedTable,
#   file = pathTestsPlatinum,
#   row.names = FALSE,
#   quote = FALSE,
#   sep = "\t",
#   fileEncoding = "UTF-8"
# )
