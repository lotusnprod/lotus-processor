source("r/log_debug.R")
log_debug("This script is the first attempt to create the tables for sql use")

start <- Sys.time()

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")

log_debug("... libraries")
library(DBI)
library(dplyr)
library(readr)
library(RPostgreSQL)
library(stringr)
library(tidyr)
source("r/sqlFromFile.R")
source("r/dbSendQueries.R")

## chemical support tables to discuss and do
## references support tables to discuss and do

drv <- RPostgres::Postgres()

log_debug("... connecting to the database")
# db <- dbConnect(
#   drv = drv,
#   dbname = "lotus",
#   user = "rutza",
#   host = "localhost",
# )

db <- DBI::dbConnect(
  drv = drv,
  dbname = dbname,
  user = user,
  host = host,
  port = port,
  password = password
)

log_debug("... listing objects")
DBI::dbListObjects(db)

log_debug("... loading files ...")

dbTypes <- readr::read_delim(
  file = "../docs/dataset.csv",
  delim = ","
) |>
  dplyr::select(
    database,
    type
  )

structureDictionary <-
  readr::read_delim(
    file = pathDataInterimDictionariesStructureDictionary,
    # n_max = 1000,
    col_types = "c"
  )

organismDictionary <-
  readr::read_delim(
    file = pathDataInterimDictionariesOrganismDictionary,
    # n_max = 1000,
    col_types = "c"
  )

## temp path
ott_taxonomy <-
  readr::read_delim(
    file = "../data/external/taxonomySource/organism/taxonomy.tsv",
    # n_max = 1000,
    col_types = "c"
  )

referenceOrganismDictionary <-
  readr::read_delim(
    file = pathDataInterimDictionariesReferenceOrganismDictionary,
    # n_max = 1000,
    col_types = cols(referenceCleanedPmid = "c")
  )

structureMetadata <-
  readr::read_delim(
    file = pathDataInterimDictionariesStructureMetadata,
    # n_max = 1000,
    col_types = "c"
  )

organismMetadata <-
  readr::read_delim(
    file = pathDataInterimDictionariesOrganismMetadata,
    # n_max = 1000,
    col_types = "c"
  )

inhouseDbMinimal <-
  readr::read_delim(
    file = pathDataInterimTablesCuratedTable,
    # n_max = 1000,
    col_types = "c"
  )

manuallyValidated <-
  readr::read_delim(
    file = "../data/validation/manuallyValidated.tsv.gz",
    col_select = any_of(colnames(inhouseDbMinimal))
  ) |>
  dplyr::mutate(curationTypeId_1 = 1)

automaticallyValidated <-
  readr::read_delim(
    file = pathDataInterimTablesAnalyzedPlatinum,
    # n_max = 1000,
    col_select = colnames(inhouseDbMinimal)
  ) |>
  dplyr::mutate(curationTypeId_2 = 2)

manuallyRemoved <-
  readr::read_delim(
    file = "../data/validation/manuallyRemoved.tsv.gz",
    col_select = any_of(colnames(inhouseDbMinimal))
  ) |>
  dplyr::mutate(curationTypeId_3 = 3)

originalTable <-
  readr::read_delim(
    file = pathDataInterimTablesOriginalTable,
    # n_max = 1000,
    col_types = "c"
  )

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
  dplyr::left_join(
    structureDictionary,
    structureMetadata
  )

organismOld <-
  dplyr::left_join(
    organismDictionary,
    organismMetadata
  )

database_type <- dbTypes |>
  dplyr::distinct(type)

database_type_old <- DBI::dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_database_type.sql")
)

database_type <- database_type_old |>
  dplyr::full_join(database_type) |>
  dplyr::distinct() |>
  dplyr::anti_join(database_type_old)

database_type <- database_type |>
  dplyr::select(name = type)

log_debug(
  "writing",
  nrow(database_type),
  "new rows to 'database_type' table"
)

DBI::dbWriteTable(
  conn = db,
  name = "database_type",
  value = database_type,
  append = TRUE,
  # overwrite = TRUE,
  row.names = FALSE
)

database_source <- originalTable |>
  dplyr::distinct(database) |>
  dplyr::left_join(dbTypes)

database_type <- DBI::dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_database_type.sql")
)

database_source <- database_source |>
  dplyr::left_join(database_type) |>
  dplyr::select(
    name = database,
    database_type_id
  )

database_source_old <- DBI::dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_database_source.sql")
)
database_source <- database_source |>
  dplyr::anti_join(database_source_old |>
    dplyr::select(name)) |>
  dplyr::distinct()

log_debug(
  "writing",
  nrow(database_source),
  "new rows to 'database_source' table"
)

DBI::dbWriteTable(
  conn = db,
  name = "database_source",
  value = database_source,
  append = TRUE,
  # overwrite = TRUE,
  row.names = FALSE
)

organism_type <- originalTable |>
  dplyr::distinct(organismType) |>
  dplyr::select(organism_type = organismType)

organism_type_old <- DBI::dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_organism_type.sql")
)
organism_type <- organism_type_old |>
  dplyr::full_join(organism_type) |>
  dplyr::distinct() |>
  dplyr::anti_join(organism_type_old)

organism_type <- organism_type |>
  dplyr::select(name = organism_type) |>
  dplyr::distinct()

log_debug(
  "writing",
  nrow(organism_type),
  "new rows to 'organism_type' table"
)

DBI::dbWriteTable(
  conn = db,
  name = "organism_type",
  value = organism_type,
  append = TRUE,
  # overwrite = TRUE,
  row.names = FALSE
)

organism_type <- DBI::dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_organism_type.sql")
)

organism_source <- originalTable |>
  dplyr::distinct(organismType, organismValue) |>
  dplyr::select(
    organism_type = organismType,
    organism_value = organismValue
  ) |>
  dplyr::left_join(organism_type) |>
  dplyr::select(
    value = organism_value,
    organism_type_id
  )

organism_source_old <- DBI::dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_organism_source.sql")
)
organism_source <- organism_source |>
  dplyr::anti_join(organism_source_old) |>
  dplyr::distinct()

log_debug(
  "writing",
  nrow(organism_source),
  "new rows to 'organism_source' table"
)

DBI::dbWriteTable(
  conn = db,
  name = "organism_source",
  value = organism_source,
  append = TRUE,
  # overwrite = TRUE,
  row.names = FALSE
)

structure_type <- originalTable |>
  dplyr::distinct(structureType) |>
  dplyr::select(structure_type = structureType)

structure_type_old <- DBI::dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_structure_type.sql")
)
structure_type <- structure_type_old |>
  dplyr::full_join(structure_type) |>
  dplyr::distinct() |>
  dplyr::anti_join(structure_type_old)

structure_type <- structure_type |>
  dplyr::select(name = structure_type) |>
  dplyr::distinct()

log_debug(
  "writing",
  nrow(structure_type),
  "new rows to 'structure_type' table"
)

DBI::dbWriteTable(
  conn = db,
  name = "structure_type",
  value = structure_type,
  append = TRUE,
  # overwrite = TRUE,
  row.names = FALSE
)

structure_type <- DBI::dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_structure_type.sql")
)

structure_source <- originalTable |>
  dplyr::distinct(structureType, structureValue) |>
  dplyr::select(
    structure_type = structureType,
    structure_value = structureValue
  ) |>
  dplyr::left_join(structure_type) |>
  dplyr::select(
    value = structure_value,
    structure_type_id
  )

structure_source_old <- DBI::dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_structure_source.sql")
)
structure_source <- structure_source |>
  dplyr::anti_join(structure_source_old) |>
  dplyr::distinct()

log_debug(
  "writing",
  nrow(structure_source),
  "new rows to 'structure_source' table"
)

DBI::dbWriteTable(
  conn = db,
  name = "structure_source",
  value = structure_source,
  append = TRUE,
  # overwrite = TRUE,
  row.names = FALSE
)

reference_type <- originalTable |>
  dplyr::distinct(referenceType) |>
  dplyr::select(reference_type = referenceType)

reference_type_old <- DBI::dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_reference_type.sql")
)
reference_type <- reference_type_old |>
  dplyr::full_join(reference_type) |>
  dplyr::distinct() |>
  dplyr::anti_join(reference_type_old)

reference_type <- reference_type |>
  dplyr::select(name = reference_type) |>
  dplyr::distinct()

log_debug(
  "writing",
  nrow(reference_type),
  "new rows to 'reference_type' table"
)

DBI::dbWriteTable(
  conn = db,
  name = "reference_type",
  value = reference_type,
  append = TRUE,
  # overwrite = TRUE,
  row.names = FALSE
)

reference_type <- DBI::dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_reference_type.sql")
)

reference_source <- originalTable |>
  dplyr::distinct(referenceType, referenceValue) |>
  dplyr::select(
    reference_type = referenceType,
    reference_value = referenceValue
  ) |>
  dplyr::left_join(reference_type) |>
  dplyr::select(
    value = reference_value,
    reference_type_id
  )

reference_source_old <- DBI::dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_reference_source.sql")
)
reference_source <- reference_source |>
  dplyr::anti_join(reference_source_old) |>
  dplyr::distinct()

log_debug(
  "writing",
  nrow(reference_source),
  "new rows to 'reference_source' table"
)

DBI::dbWriteTable(
  conn = db,
  name = "reference_source",
  value = reference_source,
  append = TRUE,
  # overwrite = TRUE,
  row.names = FALSE
)

database_source <- DBI::dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_database_source.sql")
)

organism_source <- DBI::dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_organism_source.sql")
)

structure_source <- DBI::dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_structure_source.sql")
)

reference_source <- DBI::dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_reference_source.sql")
)

data_source <- originalTable |>
  dplyr::left_join(
    database_source |>
      dplyr::select(
        database = name,
        databaseSourceId = id,
        dplyr::everything()
      )
  ) |>
  left_join(
    organism_source |>
      dplyr::select(
        organismValue = value,
        organismSourceId = id,
        dplyr::everything()
      )
  ) |>
  left_join(
    structure_source |>
      dplyr::select(
        structureValue = value,
        structureSourceId = id,
        dplyr::everything()
      )
  ) |>
  left_join(
    reference_source |>
      dplyr::select(
        referenceValue = value,
        referenceSourceId = id,
        dplyr::everything()
      )
  ) |>
  dplyr::distinct(
    databaseSourceId,
    organismSourceId,
    structureSourceId,
    referenceSourceId
  ) |>
  dplyr::select(
    database_source_id = databaseSourceId,
    organism_source_id = organismSourceId,
    structure_source_id = structureSourceId,
    reference_source_id = referenceSourceId
  )

data_source_old <- DBI::dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_data_source_full.sql")
)
data_source <- data_source |>
  dplyr::anti_join(data_source_old) |>
  dplyr::distinct()

log_debug(
  "writing",
  nrow(data_source),
  "new rows to 'data_source' table"
)

DBI::dbWriteTable(
  conn = db,
  name = "data_source",
  value = data_source,
  append = TRUE,
  # overwrite = TRUE,
  row.names = FALSE
)

data_source <- DBI::dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_data_source.sql")
)

organism_detected <- organismOld |>
  dplyr::distinct(
    organismDetected,
    organismCleaned
  ) |>
  dplyr::mutate(id = row_number()) |>
  dplyr::select(id,
    organism_detected = organismDetected,
    organism_cleaned = organismCleaned
  )

organism_cleaned <- organismOld |>
  dplyr::select(name = organismCleaned)

organism_cleaned_old <- DBI::dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_organism_cleaned.sql")
)
organism_cleaned <- organism_cleaned |>
  dplyr::anti_join(organism_cleaned_old) |>
  dplyr::distinct()

log_debug(
  "writing",
  nrow(organism_cleaned),
  "new rows to 'organism_cleaned' table"
)

DBI::dbWriteTable(
  conn = db,
  name = "organism_cleaned",
  value = organism_cleaned,
  append = TRUE,
  # overwrite = TRUE,
  row.names = FALSE
)

organism_cleaned <- DBI::dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_organism_cleaned.sql")
)

organism_synonym <- organism_detected |>
  dplyr::left_join(organism_cleaned,
    by = c("organism_cleaned" = "name")
  ) |>
  dplyr::select(
    id = id.x,
    name = organism_detected,
    organism_cleaned_id = id.y
  )

organism_synonym <- organism_synonym |>
  dplyr::select(
    name,
    organism_cleaned_id
  )

organism_synonym_old <- DBI::dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_organism_synonym.sql")
)
organism_synonym <- organism_synonym |>
  dplyr::anti_join(organism_synonym_old) |>
  dplyr::distinct()

log_debug(
  "writing",
  nrow(organism_synonym),
  "new rows to 'organism_synonym' table"
)

DBI::dbWriteTable(
  conn = db,
  name = "organism_synonym",
  value = organism_synonym,
  append = TRUE,
  # overwrite = TRUE,
  row.names = FALSE
)

organism_database <- organismOld |>
  dplyr::distinct(organismCleaned_dbTaxo) |>
  dplyr::select(name = organismCleaned_dbTaxo)

organism_database_old <- DBI::dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_organism_database.sql")
)
organism_database <- organism_database |>
  dplyr::anti_join(organism_database_old) |>
  dplyr::distinct()

log_debug(
  "writing",
  nrow(organism_database),
  "new rows to 'organism_database' table"
)

DBI::dbWriteTable(
  conn = db,
  name = "organism_database",
  value = organism_database,
  append = TRUE,
  # overwrite = TRUE,
  row.names = FALSE
)

organism_database <- DBI::dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_organism_database.sql")
)

organism_information <- organismOld |>
  dplyr::left_join(organism_cleaned,
    by = c("organismCleaned" = "name")
  ) |>
  dplyr::select(
    organism_cleaned_id = id,
    dplyr::everything()
  ) |>
  dplyr::left_join(organism_database,
    by = c("organismCleaned_dbTaxo" = "name")
  ) |>
  dplyr::select(
    organism_database_id = id,
    dplyr::everything()
  ) |>
  dplyr::distinct(organism_cleaned_id,
    organism_database_id,
    .keep_all = TRUE
  )

organism_information_old <- DBI::dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_organism_information.sql")
)
organism_information <- organism_information |>
  dplyr::anti_join(organism_information_old)

organism_information <- organism_information |>
  dplyr::select(
    organism_cleaned_id,
    organism_database_id,
    taxon_id = organismCleaned_id,
    ranks = organismCleaned_dbTaxoTaxonRanks,
    taxonomy = organismCleaned_dbTaxoTaxonomy,
    rank = organismCleaned_rank
  ) |>
  dplyr::distinct()

log_debug(
  "writing",
  nrow(organism_information),
  "new rows to 'organism_information' table"
)

DBI::dbWriteTable(
  conn = db,
  name = "organism_information",
  value = organism_information,
  append = TRUE,
  # overwrite = TRUE,
  row.names = FALSE
)

reference_cleaned <- referenceOrganismDictionary |>
  dplyr::filter(!is.na(referenceCleanedTitle)) |>
  dplyr::distinct(
    referenceCleanedDoi,
    referenceCleanedPmcid,
    referenceCleanedPmid,
    referenceCleanedTitle
  ) |>
  dplyr::left_join(referenceOrganismDictionary) |>
  dplyr::select(
    doi = referenceCleanedDoi,
    pmcid = referenceCleanedPmcid,
    pmid = referenceCleanedPmid,
    title = referenceCleanedTitle
  )

reference_cleaned_old <- DBI::dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_reference_cleaned.sql")
)
reference_cleaned <- reference_cleaned |>
  dplyr::anti_join(reference_cleaned_old) |>
  dplyr::distinct()

log_debug(
  "writing",
  nrow(reference_cleaned),
  "new rows to 'reference_cleaned' table"
)

DBI::dbWriteTable(
  conn = db,
  name = "reference_cleaned",
  value = reference_cleaned,
  append = TRUE,
  # overwrite = TRUE,
  row.names = FALSE
)

reference_cleaned <- DBI::dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_reference_cleaned.sql")
)

structure_cleaned <- structureOld |>
  dplyr::distinct(
    structureCleanedSmiles,
    structureCleanedInchi,
    structureCleanedInchikey
  ) |>
  dplyr::left_join(
    structureOld |>
      dplyr::distinct(
        structureCleanedSmiles,
        structureCleanedInchi,
        structureCleanedInchikey,
        .keep_all = TRUE
      )
  ) |>
  dplyr::distinct() |>
  dplyr::mutate(
    stereocentersTotal = as.integer(structureCleaned_stereocenters_total),
    stereocentersUnspecified = as.integer(structureCleaned_stereocenters_unspecified),
    exactMass = as.numeric(structureCleaned_exactMass),
    xlogp = as.numeric(structureCleaned_xlogp)
  ) |>
  dplyr::select(
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
  )

structure_cleaned_old <- DBI::dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_structure_cleaned.sql")
)
structure_cleaned <- structure_cleaned |>
  dplyr::anti_join(structure_cleaned_old) |>
  dplyr::distinct()

log_debug(
  "writing",
  nrow(structure_cleaned),
  "new rows to 'structure_cleaned' table"
)

DBI::dbWriteTable(
  conn = db,
  name = "structure_cleaned",
  value = structure_cleaned,
  append = TRUE,
  # overwrite = TRUE,
  row.names = FALSE
)

structure_cleaned <- DBI::dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_structure_cleaned.sql")
)

rm(
  data_source_old,
  data_source,
  referenceOrganismDictionary,
  ott_taxonomy,
  structureOld,
  structure_cleaned_old,
  structureMetadata,
  structureDictionary,
  structure_source,
  structure_source_old,
  organismOld,
  organismMetadata
)
gc()

inhouseDbMinimal_complemented <-
  dplyr::semi_join(inhouseDbMinimal, originalTable) %>%
  dplyr::mutate(curationTypeId_4 = 4) %>%
  dplyr::left_join(., manuallyValidated) %>%
  dplyr::left_join(., manuallyRemoved) %>%
  dplyr::left_join(., automaticallyValidated) %>%
  tidyr::pivot_longer(cols = (ncol(.) - 3):ncol(.)) %>%
  dplyr::filter(!is.na(value)) %>%
  dplyr::select(dplyr::everything(), -name,
    curationTypeId = value
  ) %>%
  dplyr::arrange(curationTypeId) %>%
  dplyr::distinct(
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
  )

data_cleaned_temp <- inhouseDbMinimal_complemented |>
  dplyr::left_join(organism_cleaned,
    by = c("organismCleaned" = "name")
  ) |>
  dplyr::select(
    organismCleanedId = id,
    dplyr::everything()
  ) |>
  dplyr::left_join(
    structure_cleaned |>
      dplyr::distinct(
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
  ) |>
  dplyr::select(
    structureCleanedId = id,
    dplyr::everything()
  ) %>%
  dplyr::left_join(reference_cleaned %>% distinct(id, title),
    by = c("referenceCleanedTitle" = "title")
  ) |>
  dplyr::select(
    referenceCleanedId = id,
    dplyr::everything()
  ) |>
  dplyr::distinct(
    database,
    organism_type = organismType,
    organism_value = organismValue,
    reference_type = referenceType,
    reference_value = referenceValue,
    structure_type = structureType,
    structure_value = structureValue,
    organism_cleaned_id = organismCleanedId,
    structure_cleaned_id = structureCleanedId,
    reference_cleaned_id = referenceCleanedId,
    curation_type_id = curationTypeId
  )

data_cleaned <- data_cleaned_temp |>
  dplyr::distinct(
    organism_cleaned_id,
    structure_cleaned_id,
    reference_cleaned_id,
    curation_type_id
  )

data_cleaned_old <- DBI::dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_data_cleaned.sql")
)
data_cleaned <- data_cleaned |>
  dplyr::anti_join(data_cleaned_old) |>
  dplyr::distinct()

log_debug(
  "writing",
  nrow(data_cleaned),
  "new rows to 'data_cleaned' table"
)

DBI::dbWriteTable(
  conn = db,
  name = "data_cleaned",
  value = data_cleaned,
  append = TRUE,
  # overwrite = TRUE,
  row.names = FALSE
)

data_cleaned <- DBI::dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_data_cleaned.sql")
)

data_source <- DBI::dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_data_source.sql")
)

rm(
  inhouseDbMinimal,
  inhouseDbMinimal_complemented,
  automaticallyValidated,
  originalTable
)
gc()

data_source__data_cleaned <- data_cleaned_temp |>
  dplyr::left_join(data_source |>
    dplyr::select(
      data_source_id = id,
      dplyr::everything()
    )) |>
  dplyr::left_join(data_cleaned |>
    dplyr::select(
      data_cleaned_id = id,
      dplyr::everything()
    )) |>
  dplyr::filter(!is.na(data_cleaned_id)) |>
  dplyr::distinct(data_source_id, data_cleaned_id)

data_source__data_cleaned_old <- DBI::dbGetQuery(
  conn = db,
  statement = sqlFromFile(file = "queries_db/extract_data_source__data_cleaned.sql")
)
data_source__data_cleaned <- data_source__data_cleaned |>
  dplyr::anti_join(data_source__data_cleaned_old) |>
  dplyr::distinct()

log_debug(
  "writing",
  nrow(data_source__data_cleaned),
  "new rows to 'data_source__data_cleaned' table"
)

DBI::dbWriteTable(
  conn = db,
  name = "data_source__data_cleaned",
  value = data_source__data_cleaned,
  append = TRUE,
  # overwrite = TRUE,
  row.names = FALSE
)

DBI::dbDisconnect(db)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
