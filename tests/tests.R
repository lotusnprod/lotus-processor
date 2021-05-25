# setwd(dir = "~/gitlab/opennaturalproductsdb/src")

mode_test <- TRUE

library(testthat)
library(tidyverse)
source("paths.R")
source("r/log_debug.R")

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

# write.table(
#   x = structureFull,
#   file = pathTestsStructures,
#   row.names = FALSE,
#   quote = FALSE,
#   sep = "\t",
#   fileEncoding = "UTF-8"
# )

# write.table(
#   x = referenceTableFull,
#   file = pathTestsReferences,
#   row.names = FALSE,
#   quote = FALSE,
#   sep = "\t",
#   fileEncoding = "UTF-8"
# )
