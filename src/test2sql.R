cat("This script is the first attempt to create the tables for sql use \n")

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")

cat("... libraries \n")
library(tidyverse)
library(data.table)
source("r/vroom_safe.R")

dbList <- lapply(pathDataInterimDbDir, vroom_read_safe)

dbTable <- rbindlist(l = dbList, fill = TRUE) %>%
  select(
    database,
    organismOriginal = biologicalsource,
    structureOriginal_inchi = inchi,
    structureOriginal_nominal = name,
    structureOriginal_smiles = smiles,
    referenceOriginal_authors = reference_authors,
    referenceOriginal_doi = reference_doi,
    referenceOriginal_external = reference_external,
    referenceOriginal_isbn = reference_isbn,
    referenceOriginal_journal = reference_journal,
    referenceOriginal_original = reference_original,
    referenceOriginal_pubmed = reference_pubmed,
    referenceOriginal_publishingDetails = reference_publishingDetails,
    referenceOriginal_split = reference_split,
    referenceOriginal_title = reference_title,
  ) %>%
  tibble()

source_data <- dbTable %>%
  mutate(curationStateId = 0) %>%
  mutate(id = row_number()) %>%
  group_by(database) %>%
  mutate(sourceDatabaseId = group_indices()) %>%
  ungroup() %>%
  select(
    id,
    sourceDatabaseId,
    database,
    structureName = structureOriginal_nominal,
    structureInchi = structureOriginal_inchi,
    structureSmiles = structureOriginal_smiles,
    organism = organismOriginal,
    referenceDoi = referenceOriginal_doi,
    referencePubmed = referenceOriginal_pubmed,
    referenceOriginal = referenceOriginal_original,
    referenceJournal = referenceOriginal_journal,
    referenceTitle = referenceOriginal_title,
    referenceAuthors = referenceOriginal_authors,
    referenceIsbn = referenceOriginal_isbn,
    referenceSplit = referenceOriginal_split,
    referenceExternal = referenceOriginal_external,
    referencePublishingDetails = referenceOriginal_publishingDetails,
    curationStateId
  )

rm(dbTable)

source_databases <- source_data %>%
  select(
    id = sourceDatabaseId,
    name = database
  ) %>%
  distinct()

source_data <- source_data %>%
  select(-database)

structureDictionary <-
  vroom_read_safe(path = pathDataInterimDictionariesStructureDictionary)

organismDictionary <-
  vroom_read_safe(path = pathDataInterimDictionariesOrganismDictionary)

referenceOrganismDictionary <-
  vroom_read_safe(path = pathDataInterimDictionariesReferenceOrganismDictionary)

structureMetadata <-
  vroom_read_safe(path = pathDataInterimDictionariesStructureMetadata)

organismMetadata <-
  vroom_read_safe(path = pathDataInterimDictionariesOrganismMetadata)

structureOld <-
  left_join(
    structureDictionary,
    structureMetadata
  )

organismOld <-
  left_join(
    organismDictionary,
    organismMetadata
  )

inhouseDbMinimal <-
  vroom_read_safe(path = pathDataInterimTablesCuratedTable)

detected_organisms <- organismOld %>%
  distinct(
    organismDetected,
    organismCleaned
  ) %>%
  mutate(id = row_number()) %>%
  select(
    id,
    organismDetected,
    organismCleaned
  )

cleaned_organisms <- organismOld %>%
  distinct(organismCleaned) %>%
  mutate(id = row_number()) %>%
  select(id,
    name = organismCleaned
  )

detected_organisms <- detected_organisms %>%
  left_join(.,
    cleaned_organisms,
    by = c("organismCleaned" = "name")
  ) %>%
  select(
    id = id.x,
    organismDetected,
    organismCleanedId = id.y
  )

taxonomic_databases <- organismOld %>%
  distinct(organismCleaned_dbTaxo) %>%
  group_by(organismCleaned_dbTaxo) %>%
  mutate(id = group_indices()) %>%
  select(id,
    name = organismCleaned_dbTaxo
  )

taxonomic_information <- organismOld %>%
  left_join(.,
    cleaned_organisms,
    by = c("organismCleaned" = "name")
  ) %>%
  select(
    cleanedOrganismId = id,
    everything()
  ) %>%
  left_join(.,
    taxonomic_databases,
    by = c("organismCleaned_dbTaxo" = "name")
  ) %>%
  select(
    taxonomicDatabaseId = id,
    everything()
  ) %>%
  distinct(cleanedOrganismId,
    taxonomicDatabaseId,
    .keep_all = TRUE
  ) %>%
  mutate(id = row_number()) %>%
  select(
    id,
    cleanedOrganismId,
    taxonomicDatabaseId,
    taxonomic_id = organismCleaned_id,
    ranks = organismCleaned_dbTaxoTaxonRanks,
    taxonomy = organismCleaned_dbTaxoTaxonomy,
    rank = organismCleaned_rank
  )

cleaned_references <- referenceOrganismDictionary %>%
  distinct(
    referenceCleanedDoi,
    referenceCleanedPmcid,
    referenceCleanedPmid,
    referenceCleanedTitle
  ) %>%
  mutate(id = row_number()) %>%
  select(
    id,
    doi = referenceCleanedDoi,
    pmcid = referenceCleanedPmcid,
    pmid = referenceCleanedPmid,
    title = referenceCleanedTitle
  )

detected_references <- referenceOrganismDictionary %>%
  distinct(organismDetected,
    referenceType,
    referenceValue,
    .keep_all = TRUE
  ) %>%
  mutate(id = row_number()) %>%
  select(
    id.x = id,
    organismDetected,
    referenceType,
    referenceValue,
    doi = referenceCleanedDoi,
    pmcid = referenceCleanedPmcid,
    pmid = referenceCleanedPmid,
    title = referenceCleanedTitle,
    everything()
  ) %>%
  left_join(., cleaned_references) %>%
  select(
    id = id.x,
    referenceCleanedId = id,
    organismDetected,
    referenceType,
    referenceValue,
    journal = referenceCleaned_journal,
    date = referenceCleaned_date,
    scoreCrossref = referenceCleaned_score_crossref,
    scoreDistance = referenceCleaned_score_distance,
    scoreTitleOrganism = referenceCleaned_score_titleOrganism,
    scoreComplementDate = referenceCleaned_score_complementDate,
    scoreComplementAuthor = referenceCleaned_score_complementAuthor,
    scoreComplementJournal = referenceCleaned_score_complementJournal,
    scoreComplementTotal = referenceCleaned_score_complementTotal
  ) %>%
  left_join(.,
    detected_organisms,
    by = c("organismDetected" = "organismDetected")
  ) %>%
  select(
    id = id.x,
    organismDetectedId = id.y,
    everything(),
    -organismDetected,
    -organismCleanedId
  )

cleaned_structures <- structureOld %>%
  distinct(
    structureCleanedSmiles,
    structureCleanedInchi,
    structureCleanedInchikey3D
  ) %>%
  mutate(id = row_number()) %>%
  select(
    id,
    structureCleanedSmiles,
    structureCleanedInchi,
    structureCleanedInchikey3D
  ) %>%
  left_join(
    .,
    structureOld %>% distinct(
      structureCleanedSmiles,
      structureCleanedInchi,
      structureCleanedInchikey3D,
      .keep_all = TRUE
    )
  ) %>%
  distinct() %>%
  select(
    id,
    traditionalName = structureCleaned_nameTraditional,
    iupacName = structureCleaned_nameIupac,
    inchikey = structureCleanedInchikey3D,
    inchi = structureCleanedInchi,
    smiles = structureCleanedSmiles,
    shortInchikey = structureCleaned_inchikey2D,
    molecularFormula = structureCleaned_molecularFormula,
    exactMass = structureCleaned_exactMass,
    xlogp = structureCleaned_xlogp,
    stereocentersTotal = structureCleaned_stereocenters_total,
    stereocentersUnspecified = structureCleaned_stereocenters_unspecified
  )

processed_data <- inhouseDbMinimal %>%
  left_join(
    .,
    organismDictionary %>% distinct(
      organismOriginal,
      organismDetected,
      organismCleaned
    ),
    by = c(
      "organismOriginal" = "organismOriginal",
      "organismCleaned" = "organismCleaned"
    )
  ) %>%
  left_join(.,
    detected_organisms,
    by = c("organismDetected" = "organismDetected")
  ) %>%
  select(
    organismDetectedId = id,
    everything(),
    -database,
    -organismOriginal,
    -organismCleaned
  ) %>%
  left_join(
    .,
    cleaned_structures %>%
      distinct(
        id,
        inchikey,
        inchi,
        smiles
      ),
    by = c(
      "structureCleanedSmiles" = "smiles",
      "structureCleanedInchi" = "inchi",
      "structureCleanedInchikey3D" = "inchikey"
    )
  ) %>%
  select(
    organismCleanedId,
    structureCleanedId = id,
    organismDetectedId,
    referenceType,
    referenceValue
  ) %>%
  left_join(
    .,
    detected_references
  ) %>%
  distinct(
    organismCleanedId,
    structureCleanedId,
    referenceCleanedId
  ) %>%
  mutate(id = row_number()) %>%
  select(
    id,
    organismCleanedId,
    structureCleanedId,
    referenceCleanedId
  )

rm(
  dbList,
  inhouseDbMinimal,
  organismDictionary,
  organismMetadata,
  organismOld,
  referenceOrganismDictionary,
  structureDictionary,
  structureMetadata,
  structureOld,
  temp
)
