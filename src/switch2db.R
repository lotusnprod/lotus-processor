source("r/log_debug.R")
log_debug("This script is the first attempt to create the tables for sql use")

start <- Sys.time()

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")

log_debug("... libraries")
library(data.table)
library(DBI)
library(RPostgreSQL)
library(tidyverse)
source("r/vroom_safe.R")
source("r/sqlFromFile.R")
source("r/dbSendQueries.R")

## chemical support tables to discuss and do
## references support tables to discuss and do

drv <- PostgreSQL()

log_debug("... connecting to the database")
# db <- dbConnect(
#   drv = drv,
#   dbname = "lotus",
#   user = "rutza",
#   host = "localhost",
# )

db <- dbConnect(
  drv = drv,
  dbname = "lotus",
  user = "lotusadmin",
  host = "127.0.0.1",
  port = 9999,
  password = Sys.getenv("LOTUS_DB_PWD")
)

log_debug("... listing objects")
dbListObjects(db)

log_debug("... loading files ...")

dbTypes <- read_delim(
  file = "../docs/dataset.csv",
  delim = ","
) %>%
  select(
    database,
    type
  )

structureDictionary <-
  vroom_read_safe(path = pathDataInterimDictionariesStructureDictionary)

organismDictionary <-
  vroom_read_safe(path = pathDataInterimDictionariesOrganismDictionary)

## temp path
ott_taxonomy <-
  vroom(file = "../data/external/taxonomySource/organism/taxonomy.tsv")

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

## fixed for the moment, will have to change
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
  left_join(
    structureDictionary,
    structureMetadata
  )

organismOld <-
  left_join(
    organismDictionary,
    organismMetadata
  )

originalTable <-
  vroom_read_safe(path = pathDataInterimTablesOriginalTable)

database_type <- dbTypes %>%
  distinct(type)

database_type_old <- dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_database_type.sql")
)

database_type_new <- database_type_old %>%
  full_join(., database_type) %>%
  distinct()

database_type <- database_type_new %>%
  anti_join(., database_type_old) %>%
  select(name = type)

log_debug(
  "writing",
  nrow(database_type),
  "new rows to 'database_type' table"
)

dbWriteTable(
  conn = db,
  name = "database_type",
  value = database_type,
  row.names = FALSE,
  append = TRUE
)

database_source <- originalTable %>%
  distinct(database) %>%
  left_join(., dbTypes)

database_type <- dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_database_type.sql")
)

database_source <- database_source %>%
  left_join(., database_type) %>%
  select(
    name = database,
    database_type_id
  )

database_source_old <- dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_database_source.sql")
)

database_source <- database_source %>%
  anti_join(., database_source_old %>% select(name)) %>%
  distinct()

log_debug(
  "writing",
  nrow(database_source),
  "new rows to 'database_source' table"
)

dbWriteTable(
  conn = db,
  name = "database_source",
  value = database_source,
  row.names = FALSE,
  append = TRUE
)

organism_type <- originalTable %>%
  distinct(organismType) %>%
  select(organism_type = organismType)

organism_type_old <- dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_organism_type.sql")
)

organism_type_new <- organism_type_old %>%
  full_join(., organism_type) %>%
  distinct()

organism_type <- organism_type_new %>%
  anti_join(., organism_type_old) %>%
  select(name = organism_type) %>%
  distinct()

log_debug(
  "writing",
  nrow(organism_type),
  "new rows to 'organism_type' table"
)

dbWriteTable(
  conn = db,
  name = "organism_type",
  value = organism_type,
  row.names = FALSE,
  append = TRUE
)

organism_type <- dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_organism_type.sql")
)

organism_source_old <- dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_organism_source.sql")
)

organism_source <- originalTable %>%
  distinct(organismType, organismValue) %>%
  select(
    organism_type = organismType,
    organism_value = organismValue
  ) %>%
  left_join(., organism_type) %>%
  select(
    value = organism_value,
    organism_type_id
  ) %>%
  anti_join(., organism_source_old) %>%
  distinct()

log_debug(
  "writing",
  nrow(organism_source),
  "new rows to 'organism_source' table"
)

dbWriteTable(
  conn = db,
  name = "organism_source",
  value = organism_source,
  row.names = FALSE,
  append = TRUE
)

structure_type <- originalTable %>%
  distinct(structureType) %>%
  select(structure_type = structureType)

structure_type_old <- dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_structure_type.sql")
)

structure_type_new <- structure_type_old %>%
  full_join(., structure_type) %>%
  distinct()

structure_type <- structure_type_new %>%
  anti_join(., structure_type_old) %>%
  select(name = structure_type) %>%
  distinct()

log_debug(
  "writing",
  nrow(structure_type),
  "new rows to 'structure_type' table"
)

dbWriteTable(
  conn = db,
  name = "structure_type",
  value = structure_type,
  row.names = FALSE,
  append = TRUE
)

structure_type <- dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_structure_type.sql")
)

structure_source_old <- dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_structure_source.sql")
)

structure_source <- originalTable %>%
  distinct(structureType, structureValue) %>%
  select(
    structure_type = structureType,
    structure_value = structureValue
  ) %>%
  left_join(., structure_type) %>%
  select(
    value = structure_value,
    structure_type_id
  ) %>%
  anti_join(., structure_source_old) %>%
  distinct()

log_debug(
  "writing",
  nrow(structure_source),
  "new rows to 'structure_source' table"
)

dbWriteTable(
  conn = db,
  name = "structure_source",
  value = structure_source,
  row.names = FALSE,
  append = TRUE
)

reference_type <- originalTable %>%
  distinct(referenceType) %>%
  select(reference_type = referenceType)

reference_type_old <- dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_reference_type.sql")
)

reference_type_new <- reference_type_old %>%
  full_join(., reference_type) %>%
  distinct()

reference_type <- reference_type_new %>%
  anti_join(., reference_type_old) %>%
  select(name = reference_type) %>%
  distinct()

log_debug(
  "writing",
  nrow(reference_type),
  "new rows to 'reference_type' table"
)

dbWriteTable(
  conn = db,
  name = "reference_type",
  value = reference_type,
  row.names = FALSE,
  append = TRUE
)

reference_type <- dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_reference_type.sql")
)

reference_source_old <- dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_reference_source.sql")
)

reference_source <- originalTable %>%
  distinct(referenceType, referenceValue) %>%
  select(
    reference_type = referenceType,
    reference_value = referenceValue
  ) %>%
  left_join(., reference_type) %>%
  select(
    value = reference_value,
    reference_type_id
  ) %>%
  anti_join(., reference_source_old) %>%
  distinct()

log_debug(
  "writing",
  nrow(reference_source),
  "new rows to 'reference_source' table"
)

dbWriteTable(
  conn = db,
  name = "reference_source",
  value = reference_source,
  row.names = FALSE,
  append = TRUE
)

database_source <- dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_database_source.sql")
)

organism_source <- dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_organism_source.sql")
)

structure_source <- dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_structure_source.sql")
)

reference_source <- dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_reference_source.sql")
)

data_source <- originalTable %>%
  left_join(
    .,
    database_source %>% select(
      database = name,
      databaseSourceId = id,
      everything()
    )
  ) %>%
  left_join(
    .,
    organism_source %>% select(
      organismValue = value,
      organismSourceId = id,
      everything()
    )
  ) %>%
  left_join(
    .,
    structure_source %>% select(
      structureValue = value,
      structureSourceId = id,
      everything()
    )
  ) %>%
  left_join(
    .,
    reference_source %>% select(
      referenceValue = value,
      referenceSourceId = id,
      everything()
    )
  ) %>%
  distinct(
    databaseSourceId,
    organismSourceId,
    structureSourceId,
    referenceSourceId
  ) %>%
  select(
    database_source_id = databaseSourceId,
    organism_source_id = organismSourceId,
    structure_source_id = structureSourceId,
    reference_source_id = referenceSourceId
  )

data_source_old <- dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_data_source_full.sql")
)

data_source <- data_source %>%
  anti_join(., data_source_old) %>%
  distinct()

log_debug(
  "writing",
  nrow(data_source),
  "new rows to 'data_source' table"
)

dbWriteTable(
  conn = db,
  name = "data_source",
  value = data_source,
  row.names = FALSE,
  append = TRUE
)

data_source <- dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_data_source.sql")
)

organism_detected <- organismOld %>%
  distinct(
    organismDetected,
    organismCleaned
  ) %>%
  mutate(id = row_number()) %>%
  select(id,
    organism_detected = organismDetected,
    organism_cleaned = organismCleaned
  )

organism_cleaned_old <- dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_organism_cleaned.sql")
)

organism_cleaned <- organismOld %>%
  select(name = organismCleaned) %>%
  anti_join(., organism_cleaned_old) %>%
  distinct()

log_debug(
  "writing",
  nrow(organism_cleaned),
  "new rows to 'organism_cleaned' table"
)

dbWriteTable(
  conn = db,
  name = "organism_cleaned",
  value = organism_cleaned,
  row.names = FALSE,
  append = TRUE
)

organism_cleaned <- dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_organism_cleaned.sql")
)

organism_synonym <- organism_detected %>%
  left_join(.,
    organism_cleaned,
    by = c("organism_cleaned" = "name")
  ) %>%
  select(
    id = id.x,
    name = organism_detected,
    organism_cleaned_id = id.y
  )

organism_synonym_old <- dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_organism_synonym.sql")
)

organism_synonym <- organism_synonym %>%
  select(
    name,
    organism_cleaned_id
  ) %>%
  anti_join(., organism_synonym_old) %>%
  distinct()

log_debug(
  "writing",
  nrow(organism_synonym),
  "new rows to 'organism_synonym' table"
)

dbWriteTable(
  conn = db,
  name = "organism_synonym",
  value = organism_synonym,
  row.names = FALSE,
  append = TRUE
)

organism_database_old <- dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_organism_database.sql")
)

organism_database <- organismOld %>%
  distinct(organismCleaned_dbTaxo) %>%
  select(name = organismCleaned_dbTaxo) %>%
  anti_join(., organism_database_old) %>%
  distinct()

log_debug(
  "writing",
  nrow(organism_database),
  "new rows to 'organism_database' table"
)

dbWriteTable(
  conn = db,
  name = "organism_database",
  value = organism_database,
  row.names = FALSE,
  append = TRUE
)

organism_database <- dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_organism_database.sql")
)

organism_information_old <- dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_organism_information.sql")
)

organism_information <- organismOld %>%
  left_join(.,
    organism_cleaned,
    by = c("organismCleaned" = "name")
  ) %>%
  select(
    organism_cleaned_id = id,
    everything()
  ) %>%
  left_join(.,
    organism_database,
    by = c("organismCleaned_dbTaxo" = "name")
  ) %>%
  select(
    organism_database_id = id,
    everything()
  ) %>%
  distinct(organism_cleaned_id,
    organism_database_id,
    .keep_all = TRUE
  ) %>%
  anti_join(organism_information_old) %>%
  select(
    organism_cleaned_id,
    organism_database_id,
    taxon_id = organismCleaned_id,
    ranks = organismCleaned_dbTaxoTaxonRanks,
    taxonomy = organismCleaned_dbTaxoTaxonomy,
    rank = organismCleaned_rank
  ) %>%
  distinct()

log_debug(
  "writing",
  nrow(organism_information),
  "new rows to 'organism_information' table"
)

dbWriteTable(
  conn = db,
  name = "organism_information",
  value = organism_information,
  row.names = FALSE,
  append = TRUE
)

reference_cleaned_old <- dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_reference_cleaned.sql")
)

reference_cleaned <- referenceOrganismDictionary %>%
  filter(!is.na(referenceCleanedTitle)) %>%
  distinct(
    referenceCleanedDoi,
    referenceCleanedPmcid,
    referenceCleanedPmid,
    referenceCleanedTitle
  ) %>%
  left_join(., referenceOrganismDictionary) %>%
  select(
    doi = referenceCleanedDoi,
    pmcid = referenceCleanedPmcid,
    pmid = referenceCleanedPmid,
    title = referenceCleanedTitle
  ) %>%
  anti_join(., reference_cleaned_old) %>%
  distinct()

log_debug(
  "writing",
  nrow(reference_cleaned),
  "new rows to 'reference_cleaned' table"
)

dbWriteTable(
  conn = db,
  name = "reference_cleaned",
  value = reference_cleaned,
  row.names = FALSE,
  append = TRUE
)

reference_cleaned <- dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_reference_cleaned.sql")
)

structure_cleaned_old <- dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_structure_cleaned.sql")
)

structure_cleaned <- structureOld %>%
  distinct(
    structureCleanedSmiles,
    structureCleanedInchi,
    structureCleanedInchikey
  ) %>%
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
  mutate(
    stereocentersTotal = as.integer(structureCleaned_stereocenters_total),
    stereocentersUnspecified = as.integer(structureCleaned_stereocenters_unspecified),
    exactMass = as.numeric(structureCleaned_exactMass),
    xlogp = as.numeric(structureCleaned_xlogp)
  ) %>%
  select(
    traditional_name = structureCleaned_nameTraditional,
    iupac_name = structureCleaned_nameIupac,
    inchikey = structureCleanedInchikey,
    inchikey2d = structureCleaned_inchikey2D,
    inchi = structureCleanedInchi,
    inchi2d = structureCleaned_inchi2D,
    smiles = structureCleanedSmiles,
    smiles2d = structureCleaned_smiles2D,
    stereocenters_total = stereocentersTotal,
    stereocenters_unspecified = stereocentersUnspecified,
    molecular_formula = structureCleaned_molecularFormula,
    exact_mass = exactMass,
    xlogp
  ) %>%
  anti_join(., structure_cleaned_old) %>%
  distinct()

log_debug(
  "writing",
  nrow(structure_cleaned),
  "new rows to 'structure_cleaned' table"
)

dbWriteTable(
  conn = db,
  name = "structure_cleaned",
  value = structure_cleaned,
  row.names = FALSE,
  append = TRUE
)

structure_cleaned <- dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_structure_cleaned.sql")
)

inhouseDbMinimal_complemented <-
  semi_join(inhouseDbMinimal, originalTable) %>%
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
    curationTypeId = value
  )

data_cleaned_temp <- inhouseDbMinimal_complemented %>%
  left_join(.,
    organism_cleaned,
    by = c("organismCleaned" = "name")
  ) %>%
  select(
    organismCleanedId = id,
    everything()
  ) %>%
  left_join(
    .,
    structure_cleaned %>%
      distinct(
        id,
        inchikey,
        inchi,
        smiles
      ),
    by = c(
      "structureCleanedSmiles" = "smiles",
      "structureCleanedInchi" = "inchi",
      "structureCleanedInchikey" = "inchikey"
    )
  ) %>%
  select(
    structureCleanedId = id,
    everything()
  ) %>%
  left_join(.,
    reference_cleaned %>% distinct(id, title),
    by = c("referenceCleanedTitle" = "title")
  ) %>%
  select(
    referenceCleanedId = id,
    everything()
  ) %>%
  distinct(
    database,
    organism_type = organismType,
    organism_value = organismValue,
    reference_type = referenceType,
    reference_vaue = referenceValue,
    structure_type = structureType,
    structure_value = structureValue,
    organism_cleaned_id = organismCleanedId,
    structure_cleaned_id = structureCleanedId,
    reference_cleaned_id = referenceCleanedId,
    curation_type_id = curationTypeId
  )

data_cleaned <- data_cleaned_temp %>%
  distinct(
    organism_cleaned_id,
    structure_cleaned_id,
    reference_cleaned_id,
    curation_type_id
  )

data_cleaned_old <- dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_data_cleaned.sql")
)

data_cleaned <- data_cleaned %>%
  anti_join(., data_cleaned_old) %>%
  distinct()

log_debug(
  "writing",
  nrow(data_cleaned),
  "new rows to 'data_cleaned' table"
)

dbWriteTable(
  conn = db,
  name = "data_cleaned",
  value = data_cleaned,
  row.names = FALSE,
  append = TRUE
)

data_cleaned <- dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_data_cleaned.sql")
)

data_source__data_cleaned_old <- dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_data_source__data_cleaned.sql")
)

data_source <- dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_data_source.sql")
)

data_source__data_cleaned <- data_cleaned_temp %>%
  left_join(
    .,
    data_source %>%
      select(
        data_source_id = id,
        everything()
      )
  ) %>%
  left_join(., data_cleaned %>%
    select(
      data_cleaned_id = id,
      everything()
    )) %>%
  filter(!is.na(data_cleaned_id)) %>%
  distinct(data_source_id, data_cleaned_id)

data_source__data_cleaned <- data_source__data_cleaned %>%
  anti_join(., data_source__data_cleaned_old) %>%
  distinct()

log_debug(
  "writing",
  nrow(data_source__data_cleaned),
  "new rows to 'data_source__data_cleaned' table"
)

dbWriteTable(
  conn = db,
  name = "data_source__data_cleaned",
  value = data_source__data_cleaned,
  row.names = FALSE,
  append = TRUE
)

dbDisconnect(db)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
