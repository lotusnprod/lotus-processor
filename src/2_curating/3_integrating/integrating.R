# title: "integrating bio chemo ref"

# loading
##paths
source("paths.R")

## libraries
library(tidyverse)
library(dplyr)
library(readr)

## functions
source("functions/analysis.R")
source("functions/helpers.R")

## files
### original table
originalTable <- read_delim(
  file = gzfile(description = pathDataInterimTablesOriginalTable),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
)

### organism
organismTableFull <- read_delim(
  file = gzfile(description = pathDataInterimTablesCleanedOrganismFinal),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  select(
    organismOriginal,
    organismCleaned = organismLowestTaxon,
    organismCleaned_dbTaxo = organismDbTaxo,
    organismCleaned_dbTaxoQuality = organismDbTaxoQuality,
    organismCleaned_dbTaxoTaxonId = organismTaxonId,
    organismCleaned_dbTaxo_1kingdom = organism_1_kingdom,
    organismCleaned_dbTaxo_2phylum = organism_2_phylum,
    organismCleaned_dbTaxo_3class = organism_3_class,
    organismCleaned_dbTaxo_4order = organism_4_order,
    organismCleaned_dbTaxo_5family = organism_5_family,
    organismCleaned_dbTaxo_6genus = organism_6_genus,
    organismCleaned_dbTaxo_7species = organism_7_species,
    organismCleaned_dbTaxo_8variety = organism_8_variety
  )

### structure
#### translated
translatedStructureTable <- read_delim(
  file = gzfile(description = pathDataInterimTablesTranslatedStructureFinal),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
)

#### cleaned
cleanedStructureTableFull <- read_delim(
  file = gzfile(description = pathDataInterimTablesCleanedStructureFile),
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
    structureCleaned_validatorLog = validatorLog,
    structureCleaned_molecularFormula = formulaSanitized,
    structureCleaned_exactMass = exactmassSanitized,
    structureCleaned_xlogp = xlogpSanitized
  )

### reference table
referenceTableFull <- read_delim(
  file = gzfile(description = pathDataInterimTablesCleanedReferenceFile),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
)

# splitting minimal and metadata
## structure
structureFull <-
  left_join(translatedStructureTable, cleanedStructureTableFull) %>%
  select(-structureTranslated)

structureMinimal <- structureFull %>%
  distinct(
    structureType,
    structureValue,
    structureCleanedInchi,
    structureCleanedInchikey3D,
    structureCleanedSmiles
  )

structureMetadata <- structureFull %>%
  distinct(
    structureCleanedInchi,
    structureCleanedInchikey3D,
    structureCleanedSmiles,
    structureCleaned_inchikey2D,
    structureCleaned_validatorLog,
    structureCleaned_molecularFormula,
    structureCleaned_exactMass,
    structureCleaned_xlogp,
  )

## organism
organismMinimal <- organismTableFull %>%
  distinct(
    organismOriginal,
    organismCleaned,
    organismCleaned_dbTaxo,
    organismCleaned_dbTaxoTaxonId
  )

organismMetadata <- organismTableFull %>%
  distinct(
    organismCleaned,
    organismCleaned_dbTaxo,
    organismCleaned_dbTaxoTaxonId,
    organismCleaned_dbTaxoQuality,
    organismCleaned_dbTaxo_1kingdom,
    organismCleaned_dbTaxo_2phylum,
    organismCleaned_dbTaxo_3class,
    organismCleaned_dbTaxo_4order,
    organismCleaned_dbTaxo_5family,
    organismCleaned_dbTaxo_6genus,
    organismCleaned_dbTaxo_7species,
    organismCleaned_dbTaxo_8variety
  )

## reference
referenceMinimal <- referenceTableFull %>%
  distinct(
    referenceType,
    referenceValue,
    referenceCleanedDoi,
    referenceCleanedPmcid,
    referenceCleanedPmid
  )

referenceMetadata <- referenceTableFull %>%
  distinct(
    referenceCleanedDoi,
    referenceCleanedPmcid,
    referenceCleanedPmid,
    referenceCleaned_title,
    referenceCleaned_journal,
    referenceCleaned_date,
    referenceCleaned_author,
    referenceCleaned_score_crossref,
    referenceCleaned_score_distance,
    referenceCleaned_score_titleOrganism
  )

# cleaning memory
gc(verbose = TRUE,
   reset = TRUE,
   full = TRUE)

# integrating
## structure
inhouseDbMinimal <-
  left_join(originalTable, structureMinimal)

## reference
inhouseDbMinimal <-
  left_join(inhouseDbMinimal, referenceMinimal)

## organism
inhouseDbMinimal <-
  left_join(inhouseDbMinimal, organismMinimal)

# export
## creating directories if they do not exist
ifelse(
  !dir.exists(pathDataInterimTablesCurated),
  dir.create(pathDataInterimTablesCurated),
  FALSE
)

print(x = "writing the monster table, if running fullmode, this may take a while")

## table
write.table(
  x = inhouseDbMinimal,
  file = gzfile(
    description = pathDataInterimTablesCuratedTable,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

## metadata
## structure
write.table(
  x = structureMetadata,
  file = gzfile(
    description = pathDataInterimTablesCuratedStructureMetadata,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

## organism
write.table(
  x = organismMetadata,
  file = gzfile(
    description = pathDataInterimTablesCuratedOrganismMetadata,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

### reference
write.table(
  x = referenceMetadata,
  file = gzfile(
    description = pathDataInterimTablesCuratedReferenceMetadata,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)