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

## chemical support tables to discuss and do
## references support tables to discuss and do

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
  mutate(curationTypeId_1 = 1)

automaticallyValidated <-
  vroom_read_safe(path = pathDataInterimTablesAnalysedPlatinum) %>%
  select(colnames(inhouseDbMinimal)) %>%
  mutate(curationTypeId_2 = 2)

manuallyRemoved <-
  vroom_read_safe(path = "../data/validation/manuallyRemoved.tsv.gz") %>%
  select(colnames(inhouseDbMinimal)) %>%
  mutate(curationTypeId_3 = 3)

curation_type <-
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

database_source <- originalTable %>%
  distinct(database) %>%
  mutate(databaseSourceId = row_number()) %>%
  left_join(., dbTypes) %>%
  group_by(type) %>%
  mutate(databaseTypeId = cur_group_id()) %>%
  ungroup()

organism_source <- originalTable %>%
  distinct(organismType, organismValue) %>%
  mutate(organismSourceId = row_number()) %>%
  group_by(organismType) %>%
  mutate(organismTypeId = cur_group_id()) %>%
  ungroup()

structure_source <- originalTable %>%
  distinct(structureType, structureValue) %>%
  mutate(structureSourceId = row_number()) %>%
  group_by(structureType) %>%
  mutate(structureTypeId = cur_group_id()) %>%
  ungroup()

reference_source <- originalTable %>%
  distinct(referenceType, referenceValue) %>%
  mutate(referenceSourceId = row_number()) %>%
  group_by(referenceType) %>%
  mutate(referenceTypeId = cur_group_id()) %>%
  ungroup()

database_type <- database_source %>%
  select(id = databaseTypeId,
         name = type) %>%
  distinct() %>%
  arrange(id)

organism_type <- organism_source %>%
  select(id = organismTypeId,
         name = organismType) %>%
  distinct() %>%
  arrange(id)

structure_type <- structure_source %>%
  select(id = structureTypeId,
         name = structureType) %>%
  distinct() %>%
  arrange(id)

reference_type <- reference_source %>%
  select(id = referenceTypeId,
         name = referenceType) %>%
  distinct() %>%
  arrange(id)

data_source <- originalTable %>%
  left_join(., database_source %>% select(-databaseTypeId)) %>%
  left_join(., organism_source %>% select(-organismTypeId)) %>%
  left_join(., structure_source %>% select(-structureTypeId)) %>%
  left_join(., reference_source %>% select(-referenceTypeId)) %>%
  mutate(id = row_number())

database_source <- database_source %>%
  select(id = databaseSourceId,
         name = database,
         databaseTypeId)

organism_source <- organism_source %>%
  select(id = organismSourceId,
         value = organismValue,
         organismTypeId)

structure_source <- structure_source %>%
  select(id = structureSourceId,
         value = structureValue,
         structureTypeId)

reference_source <- reference_source %>%
  select(id = referenceSourceId,
         value = referenceValue,
         referenceTypeId)

organism_detected <- organismOld %>%
  distinct(organismDetected,
           organismCleaned) %>%
  mutate(id = row_number()) %>%
  select(id,
         organismDetected,
         organismCleaned)

organism_cleaned <- organismOld %>%
  distinct(organismCleaned) %>%
  mutate(id = row_number()) %>%
  select(id,
         name = organismCleaned)

organism_synonym <- organism_detected %>%
  left_join(.,
            organism_cleaned,
            by = c("organismCleaned" = "name")) %>%
  select(id = id.x,
         name = organismDetected,
         organismCleanedId = id.y)

organism_database <- organismOld %>%
  distinct(organismCleaned_dbTaxo) %>%
  group_by(organismCleaned_dbTaxo) %>%
  mutate(id = cur_group_id()) %>%
  ungroup() %>%
  select(id,
         name = organismCleaned_dbTaxo)

organism_information <- organismOld %>%
  left_join(.,
            organism_cleaned,
            by = c("organismCleaned" = "name")) %>%
  select(organismCleanedId = id,
         everything()) %>%
  left_join(.,
            organism_database,
            by = c("organismCleaned_dbTaxo" = "name")) %>%
  select(organismDatabaseId = id,
         everything()) %>%
  distinct(organismCleanedId,
           organismDatabaseId,
           .keep_all = TRUE) %>%
  mutate(id = row_number()) %>%
  select(
    id,
    organismCleanedId,
    organismDatabaseId,
    taxonId = organismCleaned_id,
    ranks = organismCleaned_dbTaxoTaxonRanks,
    taxonomy = organismCleaned_dbTaxoTaxonomy,
    rank = organismCleaned_rank
  )

reference_cleaned <- referenceOrganismDictionary %>%
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

structure_cleaned <- structureOld %>%
  distinct(structureCleanedSmiles,
           structureCleanedInchi,
           structureCleanedInchikey) %>%
  mutate(id = row_number()) %>%
  select(id,
         structureCleanedSmiles,
         structureCleanedInchi,
         structureCleanedInchikey) %>%
  left_join(
    .,
    structureOld %>% distinct(
      structureCleanedSmiles,
      structureCleanedInchi,
      structureCleanedInchikey,
      .keep_all = TRUE
    )
  ) %>%
  distinct() %>%
  select(
    id,
    traditionalName = structureCleaned_nameTraditional,
    iupacName = structureCleaned_nameIupac,
    inchikey = structureCleanedInchikey,
    inchikey2D = structureCleaned_inchikey2D,
    inchi = structureCleanedInchi,
    inchi2D = structureCleaned_inchi2D,
    smiles = structureCleanedSmiles,
    smiles2D = structureCleaned_smiles2D,
    stereocentersTotal = structureCleaned_stereocenters_total,
    stereocentersUnspecified = structureCleaned_stereocenters_unspecified,
    molecularFormula = structureCleaned_molecularFormula,
    exactMass = structureCleaned_exactMass,
    xlogp = structureCleaned_xlogp
  )

inhouseDbMinimal_complemented <- inhouseDbMinimal %>%
  mutate(curationTypeId_4 = 4) %>%
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
    structureCleanedInchikey,
    organismCleaned,
    referenceCleanedTitle,
    .keep_all = TRUE
  ) %>%
  select(everything(),
         -name,
         curationTypeId = value)

data_cleaned_temp <- inhouseDbMinimal_complemented %>%
  left_join(.,
            organism_cleaned,
            by = c("organismCleaned" = "name")) %>%
  select(organismCleanedId = id,
         everything()) %>%
  left_join(
    .,
    structure_cleaned %>%
      distinct(id,
               inchikey,
               inchi,
               smiles),
    by = c(
      "structureCleanedSmiles" = "smiles",
      "structureCleanedInchi" = "inchi",
      "structureCleanedInchikey" = "inchikey"
    )
  ) %>%
  select(structureCleanedId = id,
         everything()) %>%
  left_join(
    .,
    reference_cleaned,
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
    curationTypeId
  )

data_cleaned <- data_cleaned_temp %>%
  distinct(organismCleanedId,
           structureCleanedId,
           referenceCleanedId,
           curationTypeId) %>%
  mutate(id = row_number()) %>%
  left_join(., data_cleaned_temp)

data_source__data_cleaned <- data_cleaned_temp %>%
  left_join(., data_source %>%
              select(dataSourceId = id,
                     everything(),
                     -type)) %>%
  left_join(., data_cleaned %>%
              select(dataCleanedId = id,
                     everything())) %>%
  filter(!is.na(dataCleanedId)) %>%
  distinct(dataSourceId, dataCleanedId) %>%
  mutate(id = row_number()) %>%
  select(id,
         dataSourceId,
         dataCleanedId)

data_cleaned <- data_cleaned %>%
  distinct(id,
           structureCleanedId,
           organismCleanedId,
           referenceCleanedId,
           curationTypeId) %>%
  select(id,
         structureCleanedId,
         organismCleanedId,
         referenceCleanedId,
         curationTypeId)

data_source <- data_source %>%
  select(id,
         databaseSourceId,
         organismSourceId,
         structureSourceId,
         referenceSourceId)

reference_cleaned <- reference_cleaned %>%
  distinct(id,
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
  organism_detected,
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
              sql = sqlFromFile("schema_db/0000_create_initial_tables.sql"))

dbListObjects(db)

cat("checking if tables match \n")
identical(dbListFields(db, "curation_type"), colnames(curation_type))
identical(dbListFields(db, "data_cleaned"), colnames(data_cleaned))
identical(dbListFields(db, "data_source"),
          colnames(data_source))
identical(
  dbListFields(db, "data_source__data_cleaned"),
  colnames(data_source__data_cleaned)
)
identical(dbListFields(db, "database_source"),
          colnames(database_source))
identical(dbListFields(db, "database_type"),
          colnames(database_type))
identical(dbListFields(db, "organism_cleaned"),
          colnames(organism_cleaned))
identical(dbListFields(db, "organism_database"),
          colnames(organism_database))
identical(dbListFields(db, "organism_information"),
         colnames(organism_information))
identical(dbListFields(db, "organism_source"),
          colnames(organism_source))
identical(dbListFields(db, "organism_synonym"),
          colnames(organism_synonym))
identical(dbListFields(db, "organism_type"),
          colnames(organism_type))
identical(dbListFields(db, "reference_cleaned"),
          colnames(reference_cleaned))
# identical(dbListFields(db, "reference_database"),
#           colnames(reference_database)) ## TO DO
# identica(dbListFields(db, "reference_information"),
#          colnames(reference_information)) ## TO DO
identical(dbListFields(db, "reference_source"),
          colnames(reference_source))
identical(dbListFields(db, "reference_type"),
          colnames(reference_type))
identical(dbListFields(db, "structure_cleaned"),
          colnames(structure_cleaned))
# identical(
#   dbListFields(db, "structure_cleaned__structure_information"),
#   colnames(structure_cleaned__structure_information)
# ) ## TO DO
# identical(dbListFields(db, "structure_information"),
#           colnames(structure_information)) ## TO DO
identical(dbListFields(db, "structure_source"),
         colnames(structure_source))
identical(dbListFields(db, "structure_type"),
          colnames(structure_type))

dbWriteTable(
  conn = db,
  name = "curation_type",
  value = curation_type,
  row.names = FALSE,
  append = TRUE
)

dbWriteTable(
  conn = db,
  name = "data_cleaned",
  value = data_cleaned,
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
  name = "data_source__data_cleaned",
  value = data_source__data_cleaned,
  row.names = FALSE,
  append = TRUE
)

dbWriteTable(
  conn = db,
  name = "database_source",
  value = database_source,
  row.names = FALSE,
  append = TRUE
)

dbWriteTable(
  conn = db,
  name = "database_type",
  value = database_type,
  row.names = FALSE,
  append = TRUE
)

dbWriteTable(
  conn = db,
  name = "organism_cleaned",
  value = organism_cleaned,
  row.names = FALSE,
  append = TRUE
)

dbWriteTable(
  conn = db,
  name = "organism_database",
  value = organism_database,
  row.names = FALSE,
  append = TRUE
)

dbWriteTable(
  conn = db,
  name = "organism_information",
  value = organism_information,
  row.names = FALSE,
  append = TRUE
)

dbWriteTable(
  conn = db,
  name = "organism_source",
  value = organism_source,
  row.names = FALSE,
  append = TRUE
)

dbWriteTable(
  conn = db,
  name = "organism_synonym",
  value = organism_synonym,
  row.names = FALSE,
  append = TRUE
)

dbWriteTable(
  conn = db,
  name = "organism_type",
  value = organism_type,
  row.names = FALSE,
  append = TRUE
)

dbWriteTable(
  conn = db,
  name = "reference_cleaned",
  value = reference_cleaned,
  row.names = FALSE,
  append = TRUE
)

# dbWriteTable(
#   conn = db,
#   name = "reference_database",
#   value = reference_database,
#   row.names = FALSE,
#   append = TRUE
# )

# dbWriteTable(
#   conn = db,
#   name = "reference_information",
#   value = reference_information,
#   row.names = FALSE,
#   append = TRUE
# )

dbWriteTable(
  conn = db,
  name = "reference_source",
  value = reference_source,
  row.names = FALSE,
  append = TRUE
)

dbWriteTable(
  conn = db,
  name = "reference_type",
  value = reference_type,
  row.names = FALSE,
  append = TRUE
)

dbWriteTable(
  conn = db,
  name = "structure_cleaned",
  value = structure_cleaned,
  row.names = FALSE,
  append = TRUE
)

# dbWriteTable(
#   conn = db,
#   name = "structure_cleaned__structure_information",
#   value = structure_cleaned__structure_information,
#   row.names = FALSE,
#   append = TRUE
# )

# dbWriteTable(
#   conn = db,
#   name = "structure_information",
#   value = structure_information,
#   row.names = FALSE,
#   append = TRUE
# )

dbWriteTable(
  conn = db,
  name = "structure_source",
  value = structure_source,
  row.names = FALSE,
  append = TRUE
)

dbWriteTable(
  conn = db,
  name = "structure_type",
  value = structure_type,
  row.names = FALSE,
  append = TRUE
)

dbDisconnect(db)

end <- Sys.time()

cat("Script finished in", format(end - start), "\n")
