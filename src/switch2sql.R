cat("This script is the first attempt to create the tables for sql use \n")

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")

cat("... libraries \n")
library(data.table)
library(DBI)
library(RSQLite)
library(tidyverse)
source("r/vroom_safe.R")
source("r/sqlFromFile.R")
source("r/dbSendQueries.R")

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

data_source <- dbTable %>%
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
    referencePublishingDetails = referenceOriginal_publishingDetails
  )

rm(dbTable)

source_databases <- data_source %>%
  select(id = sourceDatabaseId,
         name = database) %>%
  distinct()

data_source <- data_source %>%
  select(-database)

structureDictionary <-
  vroom_read_safe(path = pathDataInterimDictionariesStructureDictionary)

organismDictionary <-
  vroom_read_safe(path = pathDataInterimDictionariesOrganismDictionary)

referenceOrganismDictionary <-
  vroom_read_safe(path = pathDataInterimDictionariesReferenceOrganismDictionary) %>%
  select(-organismDetected) %>%
  left_join(.,
            organismDictionary %>%
              distinct(organismOriginal, organismDetected))
## else the organismDetected in the dictionary only corresponds to the word(organismDetected,1)

structureMetadata <-
  vroom_read_safe(path = pathDataInterimDictionariesStructureMetadata)

organismMetadata <-
  vroom_read_safe(path = pathDataInterimDictionariesOrganismMetadata)

structureOld <-
  left_join(structureDictionary,
            structureMetadata)

organismOld <-
  left_join(organismDictionary,
            organismMetadata)

inhouseDbMinimal <-
  vroom_read_safe(path = pathDataInterimTablesCuratedTable)

manuallyValidated <-
  vroom_read_safe(path = "../data/validation/manuallyValidated.tsv.gz") %>%
  select(colnames(inhouseDbMinimal)) %>%
  mutate(curationStateId_1 = 1)

automaticallyValidated <-
  vroom_read_safe(path = pathDataInterimTablesAnalysedPlatinum) %>%
  select(colnames(inhouseDbMinimal)) %>%
  mutate(curationStateId_2 = 2)

manuallyRemoved <-
  vroom_read_safe(path = "../data/validation/manuallyRemoved.tsv.gz") %>%
  select(colnames(inhouseDbMinimal)) %>%
  mutate(curationStateId_3 = 3)

curation_states <-
  data.frame(
    id = c(1, 2, 3, 4),
    name = c(
      "manual_validation",
      "automatic_validation",
      "manual_removal",
      "automatic_removal"
    )
  )

organisms_detected <- organismOld %>%
  distinct(organismDetected,
           organismCleaned) %>%
  mutate(id = row_number()) %>%
  select(id,
         organismDetected,
         organismCleaned)

organisms_cleaned <- organismOld %>%
  distinct(organismCleaned) %>%
  mutate(id = row_number()) %>%
  select(id,
         name = organismCleaned)

organisms_synonyms <- organisms_detected %>%
  left_join(.,
            organisms_cleaned,
            by = c("organismCleaned" = "name")) %>%
  select(id = id.x,
         name = organismDetected,
         organismCleanedId = id.y)

taxonomic_databases <- organismOld %>%
  distinct(organismCleaned_dbTaxo) %>%
  group_by(organismCleaned_dbTaxo) %>%
  mutate(id = group_indices()) %>%
  select(id,
         name = organismCleaned_dbTaxo)

taxonomic_information <- organismOld %>%
  left_join(.,
            organisms_cleaned,
            by = c("organismCleaned" = "name")) %>%
  select(cleanedOrganismId = id,
         everything()) %>%
  left_join(.,
            taxonomic_databases,
            by = c("organismCleaned_dbTaxo" = "name")) %>%
  select(taxonomicDatabaseId = id,
         everything()) %>%
  distinct(cleanedOrganismId,
           taxonomicDatabaseId,
           .keep_all = TRUE) %>%
  mutate(id = row_number()) %>%
  select(
    id,
    cleanedOrganismId,
    taxonomicDatabaseId,
    taxonomicId = organismCleaned_id,
    ranks = organismCleaned_dbTaxoTaxonRanks,
    taxonomy = organismCleaned_dbTaxoTaxonomy,
    rank = organismCleaned_rank
  )

references_cleaned <- referenceOrganismDictionary %>%
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

references_detected <- referenceOrganismDictionary %>%
  distinct(organismDetected,
           referenceType,
           referenceValue,
           .keep_all = TRUE) %>%
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
  left_join(., references_cleaned) %>%
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
            organisms_detected,
            by = c("organismDetected" = "organismDetected")) %>%
  select(id = id.x,
         organismDetectedId = id.y,
         everything(),
         -organismDetected)

structures_cleaned <- structureOld %>%
  distinct(structureCleanedSmiles,
           structureCleanedInchi,
           structureCleanedInchikey3D) %>%
  mutate(id = row_number()) %>%
  select(id,
         structureCleanedSmiles,
         structureCleanedInchi,
         structureCleanedInchikey3D) %>%
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
    shortInchikey = structureCleaned_inchikey2D,
    inchi = structureCleanedInchi,
    smiles = structureCleanedSmiles,
    stereocentersTotal = structureCleaned_stereocenters_total,
    stereocentersUnspecified = structureCleaned_stereocenters_unspecified,
    molecularFormula = structureCleaned_molecularFormula,
    exactMass = structureCleaned_exactMass,
    xlogp = structureCleaned_xlogp
  )

inhouseDbMinimal <- inhouseDbMinimal %>%
  mutate(curationStateId_4 = 4) %>%
  left_join(., manuallyValidated) %>%
  left_join(., manuallyRemoved) %>%
  left_join(., automaticallyValidated) %>%
  pivot_longer(cols = (ncol(.) - 3):ncol(.)) %>%
  arrange(value) %>%
  distinct(
    database,
    organismOriginal,
    referenceType,
    referenceValue,
    structureType,
    structureValue,
    structureCleanedSmiles,
    structureCleanedInchi,
    structureCleanedInchikey3D,
    organismCleaned,
    referenceCleanedTitle,
    .keep_all = TRUE
  ) %>%
  select(everything(),
         -name,
         curationStateId = value)

data_processed_temp <- inhouseDbMinimal %>%
  left_join(
    .,
    organismDictionary %>% distinct(organismOriginal,
                                    organismDetected,
                                    organismCleaned),
    by = c("organismOriginal" = "organismOriginal",
           "organismCleaned" = "organismCleaned")
  ) %>%
  left_join(.,
            organisms_synonyms,
            by = c("organismDetected" = "name")) %>%
  select(organismDetectedId = id,
         everything()) %>%
  left_join(
    .,
    structures_cleaned %>%
      distinct(id,
               inchikey,
               inchi,
               smiles),
    by = c(
      "structureCleanedSmiles" = "smiles",
      "structureCleanedInchi" = "inchi",
      "structureCleanedInchikey3D" = "inchikey"
    )
  ) %>%
  select(structureCleanedId = id,
         everything()) %>%
  left_join(.,
            references_detected) %>%
  distinct(
    database,
    organismOriginal,
    referenceType,
    referenceValue,
    structureType,
    structureValue,
    organismCleanedId,
    structureCleanedId,
    referenceCleanedId,
    curationStateId
  ) %>%
  select(organism = organismOriginal, everything())

data_processed <- data_processed_temp %>%
  distinct(organismCleanedId,
           structureCleanedId,
           referenceCleanedId,
           curationStateId) %>%
  mutate(id = row_number()) %>%
  select(id,
         structureCleanedId,
         organismCleanedId,
         referenceCleanedId,
         curationStateId)

data_processed__data_source <- data_source %>%
  left_join(source_databases, by = c("sourceDatabaseId" = "id")) %>%
  select(
    dataSourceId = id,
    structureOriginal_nominal = structureName,
    structureOriginal_inchi = structureInchi,
    structureOriginal_smiles = structureSmiles,
    organism,
    referenceOriginal_doi = referenceDoi,
    referenceOriginal_pubmed = referencePubmed,
    referenceOriginal_original = referenceOriginal,
    referenceOriginal_journal = referenceJournal,
    referenceOriginal_title = referenceTitle,
    referenceOriginal_authors = referenceAuthors,
    referenceOriginal_isbn = referenceIsbn,
    referenceOriginal_split = referenceSplit,
    referenceOriginal_external = referenceExternal,
    referenceOriginal_publishingDetails = referencePublishingDetails,
    database = name,
  ) %>%
  pivot_longer(
    cols = 2:4,
    names_to = c("drop", "structureType"),
    names_sep = "_",
    values_to = "structureValue",
    values_drop_na = TRUE
  ) %>%
  pivot_longer(
    cols = 3:12,
    names_to = c("drop2", "referenceType"),
    names_sep = "_",
    values_to = "referenceValue",
    values_drop_na = TRUE
  ) %>%
  select(-drop, -drop2) %>%
  left_join(., data_processed_temp) %>%
  left_join(., data_processed) %>%
  select(dataSourceId, dataProcessedId = id) %>%
  filter(!is.na(dataProcessedId)) %>%
  mutate(id = row_number())

rm(
  automaticallyValidated,
  manuallyRemoved,
  manuallyValidated,
  dbList,
  inhouseDbMinimal,
  organismDictionary,
  organismMetadata,
  organisms_detected,
  organismOld,
  referenceOrganismDictionary,
  references_detected,
  structureDictionary,
  structureMetadata,
  structureOld,
  temp
)

drv <- SQLite()

## TEMP
file.create("../data/processed/lotusNew.sqlite")

## TEMP
db <- dbConnect(drv = drv,
                dbname = "../data/processed/lotusNew.sqlite")

## TEMP
dbSendQueries(conn = db,
              sqlFromFile("schema_db/0000_create_initial_tables.sql"))

dbListObjects(db)

dbListFields(db, "curation_states")
dbListFields(db, "data_processed")
dbListFields(db, "data_processed__data_source")
dbListFields(db, "data_source")
dbListFields(db, "databases_source")
dbListFields(db, "databases_types")
dbListFields(db, "organisms_cleaned")
dbListFields(db, "organisms_source")
dbListFields(db, "organisms_synonyms")
dbListFields(db, "organisms_types")
dbListFields(db, "references_cleaned")
dbListFields(db, "references_source")
dbListFields(db, "references_types")
dbListFields(db, "structures_cleaned")
dbListFields(db, "structures_source")
dbListFields(db, "structures_types")
dbListFields(db, "taxonomic_databases")
dbListFields(db, "taxonomic_information")

dbWriteTable(
  conn = db,
  name = "curation_states",
  value = curation_states,
  row.names = FALSE,
  append = TRUE
)

dbWriteTable(
  conn = db,
  name = "data_processed",
  value = data_processed,
  row.names = FALSE,
  append = TRUE
)

dbWriteTable(
  conn = db,
  name = "data_processed__data_source",
  value = data_processed__data_source,
  row.names = FALSE,
  append = TRUE
)

dbWriteTable(
  conn = db,
  name = "data_source",
  value = data_source,
  row.names = FALSE,
  append = TRUE
)

dbWriteTable(
  conn = db,
  name = "organisms_cleaned",
  value = organisms_cleaned,
  row.names = FALSE,
  append = TRUE
)

dbWriteTable(
  conn = db,
  name = "organisms_synonyms",
  value = organisms_synonyms,
  row.names = FALSE,
  append = TRUE
)

dbWriteTable(
  conn = db,
  name = "references_cleaned",
  value = references_cleaned,
  row.names = FALSE,
  append = TRUE
)

dbWriteTable(
  conn = db,
  name = "source_databases",
  value = source_databases,
  row.names = FALSE,
  append = TRUE
)

dbWriteTable(
  conn = db,
  name = "structures_cleaned",
  value = structures_cleaned,
  row.names = FALSE,
  append = TRUE
)

dbWriteTable(
  conn = db,
  name = "taxonomic_databases",
  value = taxonomic_databases,
  row.names = FALSE,
  append = TRUE
)

dbWriteTable(
  conn = db,
  name = "taxonomic_information",
  value = taxonomic_information,
  row.names = FALSE,
  append = TRUE
)

dbDisconnect(db)
