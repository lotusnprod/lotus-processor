# title: "integrating bio chemo ref"

# loading
##paths
source("paths.R")

## libraries
library(tidyverse)

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

### dictionaries
#### structure
structureDictionary <- read_delim(
  file = gzfile(description = pathDataInterimDictionariesStructureDictionary),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
)

#### organism
organismDictionary <- read_delim(
  file = gzfile(description = pathDataInterimDictionariesOrganismDictionary),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
)

#### reference
referenceOrganismDictionary <- read_delim(
  file = gzfile(description = pathDataInterimDictionariesReferenceOrganismDictionary),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
)

### metadata
#### structure
structureMetadata <- read_delim(
  file = gzfile(description = pathDataInterimDictionariesStructureMetadata),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
)

#### organism
organismMetadata <- read_delim(
  file = gzfile(description = pathDataInterimDictionariesOrganismMetadata),
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

# joining previous dictionaries with metadata
## organism
organismOld <-
  left_join(organismDictionary, organismMetadata)

## structure
structureOld <-
  left_join(structureDictionary, structureMetadata)

# joining previous results with new ones
## organism
organismTableFull <- bind_rows(organismTableFull, organismOld) %>%
  distinct()

## structure
structureFull <-
  left_join(translatedStructureTable, cleanedStructureTableFull) %>%
  select(-structureTranslated)

structureFull <- bind_rows(structureFull, structureOld) %>%
  distinct()

## reference
referenceTableFull <-
  bind_rows(referenceTableFull, referenceOrganismDictionary) %>%
  distinct()

# splitting minimal and metadata
## structure
structureMinimal <- structureFull %>%
  filter(!is.na(structureCleanedInchi)) %>%
  distinct(
    structureType,
    structureValue,
    structureCleanedInchi,
    structureCleanedInchikey3D,
    structureCleanedSmiles
  )

structureMetadata <- structureFull %>%
  filter(!is.na(structureCleanedInchi)) %>%
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
  filter(!is.na(organismCleaned)) %>%
  distinct(
    organismOriginal,
    organismCleaned,
    organismCleaned_dbTaxo,
    organismCleaned_dbTaxoTaxonIds,
    organismCleaned_dbTaxoTaxonRanks,
    organismCleaned_dbTaxoTaxonomy
  )

organismMetadata <- organismTableFull %>%
  filter(!is.na(organismCleaned)) %>%
  distinct(
    organismCleaned,
    organismCleaned_dbTaxo,
    organismCleaned_dbTaxoTaxonIds,
    organismCleaned_dbTaxoTaxonRanks,
    organismCleaned_dbTaxoTaxonomy,
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
  filter(
    !is.na(referenceCleanedDoi) |
      !is.na(referenceCleanedPmcid) |
      !is.na(referenceCleanedPmid)
  ) %>%
  distinct(
    organismCleaned,
    referenceType,
    referenceValue,
    referenceCleanedDoi,
    referenceCleanedPmcid,
    referenceCleanedPmid
  )

referenceMetadata <- referenceTableFull %>%
  filter(
    !is.na(referenceCleanedDoi) |
      !is.na(referenceCleanedPmcid) |
      !is.na(referenceCleanedPmid)
  ) %>%
  distinct(
    organismCleaned,
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
  left_join(originalTable, structureMinimal) %>%
  filter(!is.na(structureCleanedInchikey3D))

## organism
inhouseDbMinimal <-
  left_join(inhouseDbMinimal, organismMinimal) %>%
  filter(!is.na(organismCleaned))

## reference
inhouseDbMinimal <-
  left_join(inhouseDbMinimal, referenceMinimal) %>%
  filter(
    !is.na(referenceCleanedDoi) |
      !is.na(referenceCleanedPmcid) |
      !is.na(referenceCleanedPmid) |
      database == "dnp_1"
  )

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

## dictionaries
### structure
ifelse(
  !dir.exists(pathDataInterimDictionariesStructure),
  dir.create(pathDataInterimDictionariesStructure),
  FALSE
)

write.table(
  x = structureMinimal,
  file = gzfile(
    description = pathDataInterimDictionariesStructureDictionary,
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
  x = referenceTableFull,
  file = gzfile(
    description = pathDataInterimDictionariesReferenceOrganismDictionary,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

## metadata
### structure
write.table(
  x = structureMetadata,
  file = gzfile(
    description = pathDataInterimDictionariesStructureMetadata,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

### organism
write.table(
  x = organismMetadata,
  file = gzfile(
    description = pathDataInterimDictionariesOrganismMetadata,
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
    description = pathDataInterimDictionariesReferenceMetadata,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
