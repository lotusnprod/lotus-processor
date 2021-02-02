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

dbTypes <- read_delim(file = "../docs/dataset.tsv",
                      delim = "\t") %>%
  select(database,
         type)

dbList <- lapply(pathDataInterimDbDir, vroom_read_safe)

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

structureOld <-
  left_join(structureDictionary,
            structureMetadata)

organismOld <-
  left_join(organismDictionary,
            organismMetadata)

dbTable <- rbindlist(l = dbList, fill = TRUE) %>%
  select(
    database,
    organismOriginal_clean = organism_clean,
    organismOriginal_dirty = organism_dirty,
    structureOriginal_inchi = structure_inchi,
    structureOriginal_nominal = structure_name,
    structureOriginal_smiles = structure_smiles,
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

originalTable <- dbTable %>%
  pivot_longer(
    cols = 7:ncol(.),
    names_to = c("drop", "referenceType"),
    names_sep = "_",
    values_to = "referenceValue",
    values_drop_na = TRUE
  ) %>%
  pivot_longer(
    cols = 4:6,
    names_to = c("drop2", "structureType"),
    names_sep = "_",
    values_to = "structureValue",
    values_drop_na = TRUE
  ) %>%
  pivot_longer(
    cols = 2:3,
    names_to = c("drop3", "organismType"),
    names_sep = "_",
    values_to = "organismValue",
    values_drop_na = TRUE
  ) %>%
  select(-drop, -drop2, -drop3) %>%
  distinct()

databases_source <- originalTable %>%
  distinct(database) %>%
  mutate(databaseSourceId = row_number()) %>%
  left_join(., dbTypes) %>%
  group_by(type) %>%
  mutate(typeId = group_indices()) %>%
  ungroup()

organisms_source <- originalTable %>%
  distinct(organismType, organismValue) %>%
  mutate(organismSourceId = row_number()) %>%
  group_by(organismType) %>%
  mutate(typeId = group_indices()) %>%
  ungroup()

structures_source <- originalTable %>%
  distinct(structureType, structureValue) %>%
  mutate(structureSourceId = row_number()) %>%
  group_by(structureType) %>%
  mutate(typeId = group_indices()) %>%
  ungroup()

references_source <- originalTable %>%
  distinct(referenceType, referenceValue) %>%
  mutate(referenceSourceId = row_number()) %>%
  group_by(referenceType) %>%
  mutate(typeId = group_indices()) %>%
  ungroup()

databases_types <- databases_source %>%
  select(id = typeId,
         name = type) %>%
  distinct() %>%
  arrange(id)

organisms_types <- organisms_source %>%
  select(id = typeId,
         name = organismType) %>%
  distinct() %>%
  arrange(id)

structures_types <- structures_source %>%
  select(id = typeId,
         name = structureType) %>%
  distinct() %>%
  arrange(id)

references_types <- references_source %>%
  select(id = typeId,
         name = referenceType) %>%
  distinct() %>%
  arrange(id)

data_source <- originalTable %>%
  left_join(., databases_source %>% select(-typeId)) %>%
  left_join(., organisms_source %>% select(-typeId)) %>%
  left_join(., structures_source %>% select(-typeId)) %>%
  left_join(., references_source %>% select(-typeId)) %>%
  mutate(id = row_number())

databases_source <- databases_source %>%
  select(id = databaseSourceId,
         name = database,
         typeId)

organisms_source <- organisms_source %>%
  select(id = organismSourceId,
         value = organismValue,
         typeId)

structures_source <- structures_source %>%
  select(id = structureSourceId,
         value = structureValue,
         typeId)

references_source <- references_source %>%
  select(id = referenceSourceId,
         value = referenceValue,
         typeId)

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
  ungroup() %>%
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
  left_join(., referenceOrganismDictionary) %>%
  select(
    id,
    doi = referenceCleanedDoi,
    pmcid = referenceCleanedPmcid,
    pmid = referenceCleanedPmid,
    title = referenceCleanedTitle,
    everything()
  )

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

inhouseDbMinimal_complemented <- inhouseDbMinimal %>%
  mutate(curationStateId_4 = 4) %>%
  left_join(., manuallyValidated) %>%
  left_join(., manuallyRemoved) %>%
  left_join(., automaticallyValidated) %>%
  pivot_longer(cols = (ncol(.) - 3):ncol(.)) %>%
  arrange(value) %>%
  distinct(
    database,
    organismType,
    organismValue,
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

data_processed_temp <- inhouseDbMinimal_complemented %>%
  left_join(.,
            organisms_cleaned,
            by = c("organismCleaned" = "name")) %>%
  select(organismCleanedId = id,
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
  left_join(
    .,
    references_cleaned,
    by = c(
      "organismType" = "organismType",
      "organismValue" = "organismValue",
      "referenceType" = "referenceType",
      "referenceValue" = "referenceValue",
      "referenceCleanedTitle" = "title"
    )
  ) %>%
  select(referenceCleanedId = id,
         everything()) %>%
  distinct(
    database,
    organismType,
    organismValue,
    referenceType,
    referenceValue,
    structureType,
    structureValue,
    organismCleanedId,
    structureCleanedId,
    referenceCleanedId,
    curationStateId
  )

data_processed <- data_processed_temp %>%
  distinct(organismCleanedId,
           structureCleanedId,
           referenceCleanedId,
           curationStateId) %>%
  mutate(id = row_number()) %>%
  left_join(., data_processed_temp)

data_processed__data_source <- data_processed_temp %>%
  left_join(., data_source %>%
              select(dataSourceId = id,
                     everything(),
                     -type)) %>%
  left_join(., data_processed %>%
              select(dataProcessedId = id,
                     everything())) %>%
  filter(!is.na(dataProcessedId)) %>%
  distinct(dataSourceId, dataProcessedId) %>%
  mutate(id = row_number()) %>%
  select(id, dataSourceId, dataProcessedId)

data_processed <- data_processed %>%
  select(id,
         structureCleanedId,
         organismCleanedId,
         referenceCleanedId,
         curationStateId)

data_source <- data_source %>%
  select(id,
         databaseSourceId,
         organismSourceId,
         structureSourceId,
         referenceSourceId)

references_cleaned <- references_cleaned %>%
  select(id,
         doi,
         pmcid,
         pmid,
         title)

rm(
  automaticallyValidated,
  dbList,
  dbTable,
  dbTypes,
  inhouseDbMinimal,
  inhouseDbMinimal_complemented,
  manuallyRemoved,
  manuallyValidated,
  organismDictionary,
  organismMetadata,
  organisms_detected,
  organismOld,
  originalTable,
  referenceOrganismDictionary,
  structureDictionary,
  structureMetadata,
  structureOld
)

drv <- SQLite()

## TEMP
file.create(lotusDB)

## TEMP
db <- dbConnect(drv = drv,
                dbname = lotusDB)
## TEMP
dbSendQueries(conn = db,
              sqlFromFile("schema_db/0000_create_initial_tables.sql"))

dbListObjects(db)

dbListFields(db, "curation_states")
colnames(curation_states)
dbListFields(db, "data_processed")
colnames(data_processed)
dbListFields(db, "data_processed__data_source")
colnames(data_processed__data_source)
dbListFields(db, "data_source")
colnames(data_source)
dbListFields(db, "databases_source")
colnames(databases_source)
dbListFields(db, "databases_types")
colnames(databases_types)
dbListFields(db, "organisms_cleaned")
colnames(organisms_cleaned)
dbListFields(db, "organisms_source")
colnames(organisms_source)
dbListFields(db, "organisms_synonyms")
colnames(organisms_synonyms)
dbListFields(db, "organisms_types")
colnames(organisms_types)
dbListFields(db, "references_cleaned")
colnames(references_cleaned)
dbListFields(db, "references_source")
colnames(references_source)
dbListFields(db, "references_types")
colnames(references_types)
dbListFields(db, "structures_cleaned")
colnames(structures_cleaned)
dbListFields(db, "structures_source")
colnames(structures_source)
dbListFields(db, "structures_types")
colnames(structures_types)
dbListFields(db, "taxonomic_databases")
colnames(taxonomic_databases)
dbListFields(db, "taxonomic_information")
colnames(taxonomic_information)

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
  name = "databases_source",
  value = databases_source,
  row.names = FALSE,
  append = TRUE
)

dbWriteTable(
  conn = db,
  name = "databases_types",
  value = databases_types,
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
  name = "organisms_source",
  value = organisms_source,
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
  name = "organisms_types",
  value = organisms_types,
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
  name = "references_source",
  value = references_source,
  row.names = FALSE,
  append = TRUE
)

dbWriteTable(
  conn = db,
  name = "references_types",
  value = references_types,
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
  name = "structures_source",
  value = structures_source,
  row.names = FALSE,
  append = TRUE
)

dbWriteTable(
  conn = db,
  name = "structures_types",
  value = structures_types,
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
