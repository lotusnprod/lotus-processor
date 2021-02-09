cat("This script is the first attempt to create the tables for postgres use \n")

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")

cat("... libraries \n")
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

cat("... connecting to the database \n")

db <- dbConnect(
  drv = drv,
  dbname = "lotus",
  user = "adafede",
  host = "localhost"
)

cat("... listing objects \n")
dbListObjects(db)

cat("... loading files ... \n")

dbTypes <- read_delim(file = "../docs/dataset.tsv",
                      delim = "\t") %>%
  select(database,
         type)

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

originalTable <-
  vroom_read_safe(path = pathDataInterimTablesOriginalTable)

structureOld <-
  left_join(structureDictionary,
            structureMetadata)

organismOld <-
  left_join(organismDictionary,
            organismMetadata)

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

curation_type_old <-
  dbGetQuery(conn = db,
             statement = "
             SELECT
             *
             FROM curation_type")

curation_type <- curation_type %>%
  anti_join(., curation_type_old) %>%
  distinct(name)

dbWriteTable(
  conn = db,
  name = "curation_type",
  value = curation_type,
  row.names = FALSE,
  append = TRUE
)

curation_type <-
  dbGetQuery(conn = db,
             statement = "
             SELECT
             *
             FROM curation_type")

database_type <- dbTypes %>%
  distinct(type)

database_type_old <- dbGetQuery(
  conn = db,
  statement = "
  SELECT
  database_type.id AS \"databaseTypeId\",
  database_type.name AS type
  FROM database_type"
)

database_type <- database_type_old %>%
  full_join(., database_type) %>%
  distinct() %>%
  anti_join(., database_type_old) %>%
  select(name = type)

dbWriteTable(
  conn = db,
  name = "database_type",
  value = database_type,
  row.names = FALSE,
  append = TRUE
)

database_type <- dbGetQuery(
  conn = db,
  statement = "
  SELECT
  database_type.id AS \"databaseTypeId\",
  database_type.name AS type
  FROM database_type"
)

database_source <- originalTable %>%
  distinct(database) %>%
  left_join(., dbTypes) %>%
  left_join(., database_type) %>%
  select(name = database,
         databaseTypeId)

database_source_old <-
  dbGetQuery(conn = db,
             statement = "
  SELECT
  name
  FROM database_source")

database_source <- database_source %>%
  anti_join(., database_source_old) %>%
  distinct()

dbWriteTable(
  conn = db,
  name = "database_source",
  value = database_source,
  row.names = FALSE,
  append = TRUE
)

organism_type <- originalTable %>%
  distinct(organismType)

organism_type_old <- dbGetQuery(
  conn = db,
  statement = "
  SELECT
  organism_type.id AS \"organismTypeId\",
  organism_type.name AS \"organismType\"
  FROM organism_type"
)

organism_type <- organism_type_old %>%
  full_join(., organism_type) %>%
  distinct() %>%
  anti_join(., organism_type_old) %>%
  select(name = organismType) %>%
  distinct()

dbWriteTable(
  conn = db,
  name = "organism_type",
  value = organism_type,
  row.names = FALSE,
  append = TRUE
)

organism_type <- dbGetQuery(
  conn = db,
  statement = "
  SELECT
  organism_type.id AS \"organismTypeId\",
  organism_type.name AS \"organismType\"
  FROM organism_type"
)

organism_source_old <- dbGetQuery(conn = db,
                                  statement = "
  SELECT
  *
  FROM organism_source")

organism_source <- originalTable %>%
  distinct(organismType, organismValue) %>%
  left_join(., organism_type) %>%
  select(value = organismValue,
         organismTypeId) %>%
  anti_join(., organism_source_old) %>%
  distinct()

dbWriteTable(
  conn = db,
  name = "organism_source",
  value = organism_source,
  row.names = FALSE,
  append = TRUE
)

structure_type <- originalTable %>%
  distinct(structureType)

structure_type_old <- dbGetQuery(
  conn = db,
  statement = "
  SELECT
  structure_type.id AS \"structureTypeId\",
  structure_type.name AS \"structureType\"
  FROM structure_type"
)

structure_type <- structure_type_old %>%
  full_join(., structure_type) %>%
  distinct() %>%
  anti_join(., structure_type_old) %>%
  select(name = structureType) %>%
  distinct()

dbWriteTable(
  conn = db,
  name = "structure_type",
  value = structure_type,
  row.names = FALSE,
  append = TRUE
)

structure_type <- dbGetQuery(
  conn = db,
  statement = "
  SELECT
  structure_type.id AS \"structureTypeId\",
    structure_type.name AS \"structureType\"
    FROM structure_type"
)

structure_source_old <- dbGetQuery(conn = db,
                                   statement = "
  SELECT
  *
  FROM structure_source")

structure_source <- originalTable %>%
  distinct(structureType, structureValue) %>%
  left_join(., structure_type) %>%
  select(value = structureValue,
         structureTypeId) %>%
  anti_join(., structure_source_old) %>%
  distinct()

dbWriteTable(
  conn = db,
  name = "structure_source",
  value = structure_source,
  row.names = FALSE,
  append = TRUE
)

reference_type <- originalTable %>%
  distinct(referenceType)

reference_type_old <- dbGetQuery(
  conn = db,
  statement = "
  SELECT
  reference_type.id AS \"referenceTypeId\",
  reference_type.name AS \"referenceType\"
  FROM reference_type"
)

reference_type <- reference_type_old %>%
  full_join(., reference_type) %>%
  distinct() %>%
  anti_join(., reference_type_old) %>%
  select(name = referenceType) %>%
  distinct()

dbWriteTable(
  conn = db,
  name = "reference_type",
  value = reference_type,
  row.names = FALSE,
  append = TRUE
)

reference_type <- dbGetQuery(
  conn = db,
  statement = "
  SELECT
  reference_type.id AS \"referenceTypeId\",
  reference_type.name AS \"referenceType\"
  FROM reference_type"
)

reference_source_old <- dbGetQuery(conn = db,
                                   statement = "
  SELECT
  *
  FROM reference_source")

reference_source <- originalTable %>%
  distinct(referenceType, referenceValue) %>%
  left_join(., reference_type) %>%
  select(value = referenceValue,
         referenceTypeId) %>%
  anti_join(., reference_source_old) %>%
  distinct()

dbWriteTable(
  conn = db,
  name = "reference_source",
  value = reference_source,
  row.names = FALSE,
  append = TRUE
)

database_source <- dbGetQuery(conn = db,
                              statement = "
                              SELECT
                              *
                              FROM database_source")

organism_source <- dbGetQuery(conn = db,
                              statement = "
                              SELECT
                              *
                              FROM organism_source")

structure_source <- dbGetQuery(conn = db,
                               statement = "
                               SELECT
                               *
                               FROM structure_source")

reference_source <- dbGetQuery(conn = db,
                               statement = "
                               SELECT
                               *
                               FROM reference_source")

data_source <- originalTable %>%
  left_join(.,
            database_source %>% select(
              database = name,
              databaseSourceId = id,
              everything()
            )) %>%
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
  distinct(databaseSourceId,
           organismSourceId,
           structureSourceId,
           referenceSourceId)

data_source_old <- dbGetQuery(conn = db,
                              statement = "
                              SELECT
                              *
                              FROM data_source")

data_source <- data_source %>%
  anti_join(., data_source_old) %>%
  distinct()

dbWriteTable(
  conn = db,
  name = "data_source",
  value = data_source,
  row.names = FALSE,
  append = TRUE
)

data_source <- dbGetQuery(conn = db,
                          statement = "
                          SELECT
                          *
                          FROM data_source")

organism_detected <- organismOld %>%
  distinct(organismDetected,
           organismCleaned) %>%
  mutate(id = row_number()) %>%
  select(id,
         organismDetected,
         organismCleaned)

organism_cleaned_old <- dbGetQuery(conn = db,
                                   statement = "
                                   SELECT
                                   name
                                   FROM organism_cleaned")

organism_cleaned <- organismOld %>%
  select(name = organismCleaned) %>%
  anti_join(., organism_cleaned_old) %>%
  distinct()

dbWriteTable(
  conn = db,
  name = "organism_cleaned",
  value = organism_cleaned,
  row.names = FALSE,
  append = TRUE
)

organism_cleaned <- dbGetQuery(conn = db,
                               statement = "
                               SELECT
                               id,
                               name
                               FROM organism_cleaned")

organism_synonym <- organism_detected %>%
  left_join(.,
            organism_cleaned,
            by = c("organismCleaned" = "name")) %>%
  select(id = id.x,
         name = organismDetected,
         organismCleanedId = id.y)

organism_synonym_old <- dbGetQuery(
  conn = db,
  statement = "
  SELECT
  name,
  \"organismCleanedId\"
  FROM organism_synonym"
)

organism_synonym <- organism_synonym %>%
  select(name, organismCleanedId) %>%
  anti_join(., organism_synonym_old) %>%
  distinct()

dbWriteTable(
  conn = db,
  name = "organism_synonym",
  value = organism_synonym,
  row.names = FALSE,
  append = TRUE
)

organism_database_old <- dbGetQuery(conn = db,
                                    statement = "
                                    SELECT
                                    name
                                    FROM organism_database")

organism_database <- organismOld %>%
  distinct(organismCleaned_dbTaxo) %>%
  select(name = organismCleaned_dbTaxo) %>%
  anti_join(., organism_database_old) %>%
  distinct()

dbWriteTable(
  conn = db,
  name = "organism_database",
  value = organism_database,
  row.names = FALSE,
  append = TRUE
)

organism_database <- dbGetQuery(conn = db,
                                statement = "
                                SELECT
                                *
                                FROM organism_database")

organism_information_old <- dbGetQuery(
  conn = db,
  statement = "
  SELECT
  \"organismCleanedId\",
  \"organismDatabaseId\",
  \"taxonId\",
  ranks,
  taxonomy,
  rank
  FROM organism_information"
)

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
  anti_join(organism_information_old) %>%
  select(
    organismCleanedId,
    organismDatabaseId,
    taxonId = organismCleaned_id,
    ranks = organismCleaned_dbTaxoTaxonRanks,
    taxonomy = organismCleaned_dbTaxoTaxonomy,
    rank = organismCleaned_rank
  ) %>%
  distinct()

dbWriteTable(
  conn = db,
  name = "organism_information",
  value = organism_information,
  row.names = FALSE,
  append = TRUE
)

reference_cleaned_old <- dbGetQuery(conn = db,
                                    statement = "
                                    SELECT
                                    *
                                    FROM reference_cleaned")

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

dbWriteTable(
  conn = db,
  name = "reference_cleaned",
  value = reference_cleaned,
  row.names = FALSE,
  append = TRUE
)

reference_cleaned <- dbGetQuery(conn = db,
                                statement = "
                                SELECT
                                *
                                FROM reference_cleaned")

structure_cleaned_old <- dbGetQuery(conn = db,
                                    statement = "
                                    SELECT
                                    *
                                    FROM structure_cleaned")

structure_cleaned <- structureOld %>%
  distinct(structureCleanedSmiles,
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
  mutate(
    stereocentersTotal = as.integer(structureCleaned_stereocenters_total),
    stereocentersUnspecified = as.integer(structureCleaned_stereocenters_unspecified),
    exactMass = as.numeric(structureCleaned_exactMass),
    xlogp = as.numeric(structureCleaned_xlogp)
  ) %>%
  select(
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
  ) %>%
  anti_join(., structure_cleaned_old) %>%
  distinct()

## DECIMALS PROBLEM WITH PGLOADER

# dbWriteTable(
#   conn = db,
#   name = "structure_cleaned",
#   value = structure_cleaned,
#   row.names = FALSE,
#   append = TRUE
# )

structure_cleaned <- dbGetQuery(conn = db,
                                statement = "
                                SELECT
                                *
                                FROM structure_cleaned")

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
  left_join(.,
            reference_cleaned %>% distinct(id, title),
            by = c("referenceCleanedTitle" = "title")) %>%
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
           curationTypeId)

data_cleaned_old <- dbGetQuery(conn = db,
                               statement = "
                               SELECT
                               *
                               FROM data_cleaned")

data_cleaned <- data_cleaned %>%
  anti_join(., data_cleaned_old) %>%
  distinct()

dbWriteTable(
  conn = db,
  name = "data_cleaned",
  value = data_cleaned,
  row.names = FALSE,
  append = TRUE
)

data_cleaned <- dbGetQuery(conn = db,
                           statement = "
                           SELECT
                           *
                           FROM data_cleaned")

data_source__data_cleaned_old <-
  dbGetQuery(conn = db,
             statement = "
             SELECT
             *
             FROM data_source__data_cleaned")

data_source <- dbGetQuery(
  conn = db,
  statement = "
  SELECT data_source.id,
  database_source.name   AS database,
  organism_source.value  AS \"organismValue\",
  organism_type.name     AS \"organismType\",
  reference_source.value AS \"referenceValue\",
  reference_type.name    AS \"referenceType\",
  structure_source.value AS \"structureValue\",
  structure_type.name    AS \"structureType\"
  FROM data_source
  LEFT JOIN database_source
  ON data_source.\"databaseSourceId\" = database_source.id
  LEFT JOIN organism_source
  ON data_source.\"organismSourceId\" = organism_source.id
  LEFT JOIN organism_type
  ON organism_source.\"organismTypeId\" = organism_type.id
  LEFT JOIN reference_source
  ON data_source.\"referenceSourceId\" = reference_source.id
  LEFT JOIN reference_type
  ON reference_source.\"referenceTypeId\" = reference_type.id
  LEFT JOIN structure_source
  ON data_source.\"structureSourceId\" = structure_source.id
  LEFT JOIN structure_type
  ON structure_source.\"structureTypeId\" = structure_type.id"
)

data_source__data_cleaned <- data_cleaned_temp %>%
  left_join(.,
            data_source %>%
              select(dataSourceId = id,
                     everything())) %>%
  left_join(., data_cleaned %>%
              select(dataCleanedId = id,
                     everything())) %>%
  filter(!is.na(dataCleanedId)) %>%
  distinct(dataSourceId, dataCleanedId) %>%
  anti_join(., data_source__data_cleaned_old) %>%
  distinct()

dbWriteTable(
  conn = db,
  name = "data_source__data_cleaned",
  value = data_source__data_cleaned,
  row.names = FALSE,
  append = TRUE
)

dbDisconnect(db)

end <- Sys.time()

cat("Script finished in", format(end - start), "\n")
