source("r/log_debug.R")
log_debug(
  "This script aligns all sources with previous results and outputs files \n",
  "containing new entries for editing"
)

start <- Sys.time()

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")

log_debug("... functions")
source("r/split_data_table_quote.R")
source("r/sqlFromFile.R")
source("r/standardizing_original.R")

log_debug("loading ...")
log_debug("... libraries")
library(data.table)
library(DBI)
library(dplyr)
library(readr)
library(tidyr)

log_debug("... files ...")
log_debug("... DBs")

if (mode == "full" | mode == "custom") {
  if (ssot_access == TRUE) {
    library(RPostgreSQL)

    drv <- RPostgreSQL::PostgreSQL()

    log_debug("... connecting to the database")
    # db <- DBI::dbConnect(
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

    log_debug("... listing remote objects")
    DBI::dbListObjects(db)

    log_debug("... extracting already processed data")
    oldTable <- DBI::dbGetQuery(
      conn = db,
      statement = sqlFromFile("queries_db/extract_data_source.sql")
    )
  } else {
    oldTable <- data.frame() |>
      dplyr::mutate(
        database = NA,
        organism_value = NA,
        organism_type = NA,
        reference_value = NA,
        reference_type = NA,
        structure_value = NA,
        structure_type = NA
      ) |>
      dplyr::mutate_all(as.character)
  }
  log_debug("... list of source databases")
  dbList <- lapply(
    pathDataInterimDbDir,
    readr::read_delim,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  )

  log_debug("... dictionaries ...")
  if (file.exists(pathDataInterimDictionariesStructureDictionary)) {
    log_debug("... structures")
    structureDictionary <-
      readr::read_delim(
        file = pathDataInterimDictionariesStructureDictionary,
        delim = "\t",
        col_types = cols(.default = "c"),
        locale = locales
      )
  }

  if (file.exists(pathDataInterimDictionariesStructureAntiDictionary)) {
    log_debug("... previously unsucessfully querried structures")
    structureAntiDictionary <-
      readr::read_delim(
        file = pathDataInterimDictionariesStructureAntiDictionary,
        delim = "\t",
        col_types = cols(.default = "c"),
        locale = locales
      )
  }

  if (file.exists(pathDataInterimDictionariesOrganismDictionary)) {
    log_debug("... organisms")
    organismDictionary <-
      readr::read_delim(
        file = pathDataInterimDictionariesOrganismDictionary,
        delim = "\t",
        col_types = cols(.default = "c"),
        locale = locales
      )
  }

  if (file.exists(pathDataInterimDictionariesReferenceDictionary)) {
    log_debug("... references")
    referenceDictionary <-
      readr::read_delim(
        file = pathDataInterimDictionariesReferenceDictionary,
        delim = "\t",
        col_types = cols(.default = "c"),
        locale = locales
      ) |>
      data.table::data.table()
  }

  log_debug("renaming and selecting columns")
  dbTable <- data.table::rbindlist(l = dbList, fill = TRUE) |>
    data.frame() |>
    ## see https://github.com/lotusnprod/lotus-processor/issues/58
    ## removing for now
    dplyr::filter(database != "metabolights")

  dbTable[setdiff(
    x = accepted_fields,
    y = names(dbTable)
  )] <- NA

  dbTable <- dbTable |>
    dplyr::select(
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
    ) |>
    dplyr::mutate(
      referenceOriginal_pubmed = as.character(referenceOriginal_pubmed)
    )

  if (mode == "full") {
    log_debug("sampling rows for test mode")
    "%ni%" <- Negate("%in%")
    set.seed(
      seed = 42,
      kind = "Mersenne-Twister",
      normal.kind = "Inversion"
    )
    dbTable_sampled_1 <- dbTable |>
      dplyr::filter(database %ni% forbidden_export) |>
      dplyr::filter(is.na(referenceOriginal_external)) |>
      dplyr::filter(
        !(!is.na(structureOriginal_inchi) &
          !is.na(structureOriginal_nominal)) |
          !(!is.na(structureOriginal_smiles) &
            !is.na(structureOriginal_nominal))
      ) |> ## to avoid too many names (long for CI)
      dplyr::sample_n(size = 490)
    set.seed(
      seed = 42,
      kind = "Mersenne-Twister",
      normal.kind = "Inversion"
    )
    dbTable_sampled_2 <- dbTable |>
      dplyr::filter(database %ni% forbidden_export) |>
      dplyr::filter(!is.na(organismOriginal_dirty)) |>
      dplyr::filter(is.na(referenceOriginal_external)) |>
      dplyr::filter(
        !(!is.na(structureOriginal_inchi) &
          !is.na(structureOriginal_nominal)) |
          !(!is.na(structureOriginal_smiles) &
            !is.na(structureOriginal_nominal))
      ) |> ## to avoid too many names (long for CI)
      dplyr::sample_n(size = 10)
    dbTable_sampled <-
      dplyr::bind_rows(dbTable_sampled_1, dbTable_sampled_2)
    originalTable_sampled <- dbTable_sampled |>
      dplyr::select(database, dplyr::everything()) |>
      tidyr::pivot_longer(
        cols = 7:16,
        names_to = c("drop", "referenceType"),
        names_sep = "_",
        values_to = "referenceValue",
        values_drop_na = TRUE
      ) |>
      tidyr::pivot_longer(
        cols = 4:6,
        names_to = c("drop2", "structureType"),
        names_sep = "_",
        values_to = "structureValue",
        values_drop_na = TRUE
      ) |>
      tidyr::pivot_longer(
        cols = 2:3,
        names_to = c("drop3", "organismType"),
        names_sep = "_",
        values_to = "organismValue",
        values_drop_na = TRUE
      ) |>
      dplyr::select(-drop, -drop2, -drop3) |>
      dplyr::distinct()
  }
} else {
  dbTable <- readr::read_delim(
    file = pathTestsFile,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  ) |>
    tidyr::pivot_wider(
      names_from = "organismType",
      values_from = "organismValue",
      names_prefix = "organismOriginal_"
    ) |>
    tidyr::pivot_wider(
      names_from = "structureType",
      values_from = "structureValue",
      names_prefix = "structureOriginal_"
    ) |>
    tidyr::pivot_wider(
      names_from = "referenceType",
      values_from = "referenceValue",
      names_prefix = "referenceOriginal_"
    ) |>
    tidyr::unnest() |>
    dplyr::mutate(
      referenceOriginal_authors = NA,
      referenceOriginal_external = NA,
      referenceOriginal_isbn = NA,
      referenceOriginal_journal = NA
    ) |>
    dplyr::select(
      database,
      organismOriginal_clean,
      organismOriginal_dirty,
      structureOriginal_inchi,
      structureOriginal_nominal,
      structureOriginal_smiles,
      referenceOriginal_authors,
      referenceOriginal_doi,
      referenceOriginal_external,
      referenceOriginal_isbn,
      referenceOriginal_journal,
      referenceOriginal_original,
      referenceOriginal_pubmed,
      referenceOriginal_publishingDetails,
      referenceOriginal_split,
      referenceOriginal_title,
    )
}

log_debug("... full original table")
originalTable <- dbTable |>
  dplyr::select(database, dplyr::everything()) |>
  tidyr::pivot_longer(
    cols = 7:16,
    names_to = c("drop", "referenceType"),
    names_sep = "_",
    values_to = "referenceValue",
    values_drop_na = TRUE
  ) |>
  dplyr::select(-drop) |>
  tidyr::pivot_longer(
    cols = 4:6,
    names_to = c("drop2", "structureType"),
    names_sep = "_",
    values_to = "structureValue",
    values_drop_na = TRUE
  ) |>
  dplyr::select(-drop2) |>
  tidyr::pivot_longer(
    cols = 2:3,
    names_to = c("drop3", "organismType"),
    names_sep = "_",
    values_to = "organismValue",
    values_drop_na = TRUE
  ) |>
  dplyr::select(-drop3) |>
  dplyr::distinct()

if (mode == "full" | mode == "custom") {
  log_debug("new entries only ...")
  originalTable <- dplyr::anti_join(
    originalTable,
    oldTable |>
      dplyr::select(
        database,
        organismValue = organism_value,
        organismType = organism_type,
        referenceValue = reference_value,
        referenceType = reference_type,
        structureValue = structure_value,
        structureType = structure_type
      )
  )
}

log_debug("ensuring proper encoding ...")
originalTable$organismValue <-
  iconv(
    x = originalTable$organismValue,
    from = "UTF-8",
    to = "UTF-8",
    sub = ""
  )

originalTable$referenceValue <-
  iconv(
    x = originalTable$referenceValue,
    from = "UTF-8",
    to = "UTF-8",
    sub = ""
  )

originalTable$structureValue <-
  iconv(
    x = originalTable$structureValue,
    from = "UTF-8",
    to = "UTF-8",
    sub = ""
  )

log_debug("keeping entries not previously curated only ...")
log_debug("... inchi table")
structureTable_inchi <- originalTable |>
  dplyr::filter(structureType == "inchi") |>
  dplyr::filter(!is.na(structureValue)) |>
  dplyr::distinct(structureValue)

if (mode != "test") {
  if (file.exists(pathDataInterimDictionariesStructureDictionary)) {
    structureTable_inchi <-
      dplyr::anti_join(
        x = structureTable_inchi,
        y = structureDictionary
      )
  }

  if (file.exists(pathDataInterimDictionariesStructureAntiDictionary)) {
    structureTable_inchi <-
      dplyr::anti_join(
        x = structureTable_inchi,
        y = structureAntiDictionary
      )
  }
}

structureTable_inchi <- structureTable_inchi |>
  dplyr::select(structureOriginal_inchi = structureValue)

if (nrow(structureTable_inchi) == 0) {
  structureTable_inchi <- data.frame(structureOriginal_inchi = NA)
}

log_debug("... smiles table")
structureTable_smiles <- originalTable |>
  dplyr::filter(structureType == "smiles") |>
  dplyr::filter(!is.na(structureValue)) |>
  dplyr::distinct(structureValue)

if (mode != "test") {
  if (file.exists(pathDataInterimDictionariesStructureDictionary)) {
    structureTable_smiles <-
      dplyr::anti_join(
        x = structureTable_smiles,
        y = structureDictionary
      )
  }

  if (file.exists(pathDataInterimDictionariesStructureAntiDictionary)) {
    structureTable_smiles <-
      dplyr::anti_join(
        x = structureTable_smiles,
        y = structureAntiDictionary
      )
  }
}

structureTable_smiles <- structureTable_smiles |>
  dplyr::select(structureOriginal_smiles = structureValue)

if (nrow(structureTable_smiles) == 0) {
  structureTable_smiles <- data.frame(structureOriginal_smiles = NA)
}

log_debug("... chemical names table")
structureTable_nominal <- originalTable |>
  dplyr::filter(structureType == "nominal") |>
  dplyr::filter(!is.na(structureValue)) |>
  dplyr::distinct(structureValue)

if (mode != "test") {
  if (file.exists(pathDataInterimDictionariesStructureDictionary)) {
    structureTable_nominal <-
      dplyr::anti_join(
        x = structureTable_nominal,
        y = structureDictionary
      )
  }

  if (file.exists(pathDataInterimDictionariesStructureAntiDictionary)) {
    structureTable_nominal <-
      dplyr::anti_join(
        x = structureTable_nominal,
        y = structureAntiDictionary
      )
  }
}

structureTable_nominal <- structureTable_nominal |>
  dplyr::select(structureOriginal_nominal = structureValue)

if (nrow(structureTable_nominal) == 0) {
  structureTable_nominal <- data.frame(structureOriginal_nominal = NA)
}

log_debug("... structures table")
structureTable_full <-
  dplyr::bind_rows(
    structureTable_inchi |>
      dplyr::mutate(structureType = "inchi") |>
      dplyr::select(structureType, structureValue = structureOriginal_inchi),
    structureTable_smiles |>
      dplyr::mutate(structureType = "smiles") |>
      dplyr::select(structureType, structureValue = structureOriginal_smiles),
    structureTable_nominal |>
      dplyr::mutate(structureType = "nominal") |>
      dplyr::select(structureType, structureValue = structureOriginal_nominal)
  ) |>
  dplyr::distinct()

if (mode == "custom") {
  structureTable_nominal <- data.frame(structureOriginal_nominal = NA)
}

if (nrow(structureTable_full) == 0) {
  structureTable_full[1, ] <- NA
}

log_debug("... organisms tables ...")
log_debug("... clean")
organismTable_clean <- originalTable |>
  dplyr::filter(organismType == "clean") |>
  dplyr::filter(!is.na(organismValue)) |>
  dplyr::distinct(organismValue)

if (nrow(organismTable_clean) == 0) {
  organismTable_clean[1, "organismValue"] <- NA
}

log_debug("... dirty")
organismTable_dirty <- originalTable |>
  dplyr::filter(organismType == "dirty") |>
  dplyr::filter(!is.na(organismValue)) |>
  dplyr::distinct(organismValue)

if (nrow(organismTable_dirty) == 0) {
  organismTable_dirty[1, "organismValue"] <- NA
}

if (mode != "test") {
  if (file.exists(pathDataInterimDictionariesOrganismDictionary)) {
    organismTable_clean <- dplyr::anti_join(
      x = organismTable_clean,
      y = organismDictionary
    )

    organismTable_dirty <- dplyr::anti_join(
      x = organismTable_dirty,
      y = organismDictionary
    )
  }
}

organismTable_clean <- organismTable_clean |>
  dplyr::select(organismOriginal_clean = organismValue)

organismTable_dirty <- organismTable_dirty |>
  dplyr::select(organismOriginal_dirty = organismValue) |>
  data.table::data.table()

log_debug("... structures table")
organismTable_full <- dplyr::bind_rows(
  organismTable_clean |>
    dplyr::mutate(organismType = "clean") |>
    dplyr::select(organismType, organismValue = organismOriginal_clean),
  organismTable_dirty |>
    dplyr::mutate(organismType = "dirty") |>
    dplyr::select(organismType, organismValue = organismOriginal_dirty)
) |>
  dplyr::distinct()

log_debug("... references tables ...")
log_debug("... DOI table")
referenceTable_doi <- dbTable |>
  dplyr::filter(!is.na(referenceOriginal_doi)) |>
  dplyr::distinct(referenceOriginal_doi) |>
  dplyr::select(referenceOriginal = referenceOriginal_doi) |>
  dplyr::mutate_all(as.character)

if (mode != "test") {
  if (file.exists(pathDataInterimDictionariesReferenceDictionary)) {
    referenceTable_doi <- dplyr::anti_join(
      x = referenceTable_doi |>
        mutate(origin = "doi"),
      y = referenceDictionary
    ) |>
      dplyr::select(-origin)
  }
}

referenceTable_doi <- referenceTable_doi |>
  dplyr::select(referenceOriginal_doi = referenceOriginal)

if (nrow(referenceTable_doi) == 0) {
  referenceTable_doi[1, ] <- NA
}

log_debug("... PMID table")
referenceTable_pubmed <- dbTable |>
  dplyr::filter(!is.na(referenceOriginal_pubmed)) |>
  dplyr::distinct(referenceOriginal_pubmed) |>
  dplyr::select(referenceOriginal = referenceOriginal_pubmed) |>
  dplyr::mutate_all(as.character)

if (mode != "test") {
  if (file.exists(pathDataInterimDictionariesReferenceDictionary)) {
    referenceTable_pubmed <- dplyr::anti_join(
      x = referenceTable_pubmed |>
        dplyr::mutate(origin = "pubmed"),
      y = referenceDictionary
    ) |>
      dplyr::select(-origin)
  }
}

referenceTable_pubmed <- referenceTable_pubmed |>
  dplyr::select(referenceOriginal_pubmed = referenceOriginal)

if (nrow(referenceTable_pubmed) == 0) {
  referenceTable_pubmed <- rbind(referenceTable_pubmed, list(NA))
}

if (nrow(referenceTable_pubmed) != 1) {
  row.names(referenceTable_pubmed) <-
    referenceTable_pubmed$referenceOriginal_pubmed
}

log_debug("... reference title table")
referenceTable_title <- dbTable |>
  dplyr::filter(is.na(referenceOriginal_doi)) |>
  dplyr::filter(is.na(referenceOriginal_pubmed)) |>
  dplyr::filter(!is.na(referenceOriginal_title)) |>
  dplyr::distinct(referenceOriginal_title) |>
  dplyr::select(referenceOriginal = referenceOriginal_title) |>
  dplyr::mutate_all(as.character)

if (mode != "test") {
  if (file.exists(pathDataInterimDictionariesReferenceDictionary)) {
    referenceTable_title <- dplyr::anti_join(
      x = referenceTable_title |>
        mutate(origin = "title"),
      y = referenceDictionary
    ) |>
      dplyr::select(-origin)
  }
}

referenceTable_title <- referenceTable_title |>
  dplyr::select(referenceOriginal_title = referenceOriginal)

if (nrow(referenceTable_title) == 0) {
  referenceTable_title[1, ] <- NA
}

referenceTable_title <- referenceTable_title |>
  data.table::data.table()

log_debug(".. reference publishing details table")
referenceTable_publishingDetails <- dbTable |>
  dplyr::filter(is.na(referenceOriginal_doi)) |>
  dplyr::filter(is.na(referenceOriginal_pubmed)) |>
  dplyr::filter(is.na(referenceOriginal_title)) |>
  dplyr::filter(!is.na(referenceOriginal_publishingDetails)) |>
  dplyr::distinct(referenceOriginal_publishingDetails) |>
  dplyr::select(referenceOriginal = referenceOriginal_publishingDetails) |>
  dplyr::mutate_all(as.character)

if (mode != "test") {
  if (file.exists(pathDataInterimDictionariesReferenceDictionary)) {
    referenceTable_publishingDetails <-
      dplyr::anti_join(
        x = referenceTable_publishingDetails |>
          mutate(origin = "publishingDetails"),
        y = referenceDictionary
      ) |>
      dplyr::select(-origin)
  }
}

referenceTable_publishingDetails <-
  referenceTable_publishingDetails |>
  dplyr::select(referenceOriginal_publishingDetails = referenceOriginal)

if (nrow(referenceTable_publishingDetails) == 0) {
  referenceTable_publishingDetails <-
    rbind(referenceTable_publishingDetails, list(NA))
}

referenceTable_publishingDetails <-
  referenceTable_publishingDetails |>
  data.table::data.table()

log_debug("... reference split table")
referenceTable_split <- dbTable %>%
  dplyr::filter(is.na(referenceOriginal_doi)) |>
  dplyr::filter(is.na(referenceOriginal_pubmed)) |>
  dplyr::filter(is.na(referenceOriginal_title)) |>
  dplyr::filter(is.na(referenceOriginal_publishingDetails)) |>
  dplyr::filter(!is.na(referenceOriginal_split)) |>
  dplyr::distinct(referenceOriginal_split) |>
  dplyr::select(referenceOriginal = referenceOriginal_split) |>
  dplyr::mutate_all(as.character)

if (mode != "test") {
  if (file.exists(pathDataInterimDictionariesReferenceDictionary)) {
    referenceTable_split <-
      dplyr::anti_join(
        x = referenceTable_split |>
          mutate(origin = "split"),
        y = referenceDictionary
      ) |>
      dplyr::select(-origin)
  }
}

referenceTable_split <- referenceTable_split |>
  dplyr::select(referenceOriginal_split = referenceOriginal)

if (nrow(referenceTable_split) == 0) {
  referenceTable_split <-
    rbind(referenceTable_split, list(NA))
}

referenceTable_split <- referenceTable_split |>
  data.table::data.table()

log_debug("... original references table")
referenceTable_original <- dbTable |>
  dplyr::filter(is.na(referenceOriginal_doi)) |>
  dplyr::filter(is.na(referenceOriginal_pubmed)) |>
  dplyr::filter(is.na(referenceOriginal_title)) |>
  dplyr::filter(is.na(referenceOriginal_publishingDetails)) |>
  dplyr::filter(is.na(referenceOriginal_split)) |>
  dplyr::filter(!is.na(referenceOriginal_original)) |>
  dplyr::distinct(referenceOriginal_original) |>
  dplyr::select(referenceOriginal = referenceOriginal_original) |>
  dplyr::mutate_all(as.character) |>
  data.table::data.table()

if (mode != "test") {
  if (file.exists(pathDataInterimDictionariesReferenceDictionary)) {
    referenceTable_original <-
      dplyr::anti_join(
        x = referenceTable_original |>
          mutate(origin = "original"),
        y = referenceDictionary
      ) |>
      dplyr::select(-origin)
  }
}

referenceTable_original <- referenceTable_original |>
  dplyr::select(referenceOriginal_original = referenceOriginal)

if (nrow(referenceTable_original) == 0) {
  referenceTable_original <-
    rbind(referenceTable_original, list(NA))
}

log_debug("... full references table")
referenceTable_full <- originalTable |>
  dplyr::select(
    organismType,
    organismValue,
    referenceType,
    referenceValue
  ) |>
  dplyr::distinct()

# data_source <- DBI::dbGetQuery(
#   conn = db,
#   statement = "SELECT * FROM data_source"
# )

log_debug("ensuring directories exist ...")
create_dir(export = pathDataInterimTablesOriginalOrganism)
create_dir(export = pathDataInterimTablesOriginalReference)
create_dir(export = pathDataInterimTablesOriginalStructure)

##### with rm
create_dir_with_rm(export = pathDataInterimTablesOriginalOrganism)
create_dir_with_rm(
  export = pathDataInterimTablesOriginalReferenceOriginalFolder
)
create_dir_with_rm(
  export = pathDataInterimTablesOriginalReferencePublishingDetailsFolder
)
create_dir_with_rm(export = pathDataInterimTablesOriginalReferenceSplitFolder)
create_dir_with_rm(export = pathDataInterimTablesOriginalReferenceTitleFolder)

log_debug("exporting ...")
if (nrow(organismTable_clean) != 0) {
  readr::write_delim(
    x = organismTable_clean,
    delim = "\t",
    file = gzfile(
      description = pathDataInterimTablesOriginalOrganismFile,
      compression = 9,
      encoding = "UTF-8"
    ),
    na = "",
    quote = "none"
  )
  ## because gnverify does not parse quotes
}

log_debug(pathDataInterimTablesOriginal)
readr::write_delim(
  x = organismTable_full,
  delim = "\t",
  file = pathDataInterimTablesOriginalOrganismFull,
  na = ""
)

log_debug(pathDataInterimTablesOriginalReferenceDoi)
readr::write_delim(
  x = referenceTable_doi,
  delim = "\t",
  file = pathDataInterimTablesOriginalReferenceDoi,
  na = ""
)

log_debug(pathDataInterimTablesOriginalReferenceOriginalFolder)
split_data_table_quote(
  x = referenceTable_original,
  no_rows_per_frame = 1000,
  text = "",
  path_to_store = pathDataInterimTablesOriginalReferenceOriginalFolder
)

log_debug(pathDataInterimTablesOriginalReferencePubmed)
readr::write_delim(
  x = referenceTable_pubmed,
  delim = "\t",
  file = pathDataInterimTablesOriginalReferencePubmed,
  na = ""
)

log_debug(pathDataInterimTablesOriginalReferencePublishingDetailsFolder)
split_data_table_quote(
  x = referenceTable_publishingDetails,
  no_rows_per_frame = 1000,
  text = "",
  path_to_store = pathDataInterimTablesOriginalReferencePublishingDetailsFolder
)

log_debug(pathDataInterimTablesOriginalReferenceSplitFolder)
split_data_table_quote(
  x = referenceTable_split,
  no_rows_per_frame = 1000,
  text = "",
  path_to_store = pathDataInterimTablesOriginalReferenceSplitFolder
)

log_debug(pathDataInterimTablesOriginalReferenceTitleFolder)
split_data_table_quote(
  x = referenceTable_title,
  no_rows_per_frame = 1000,
  text = "",
  path_to_store = pathDataInterimTablesOriginalReferenceTitleFolder
)

log_debug(pathDataInterimTablesOriginalReferenceFull)
readr::write_delim(
  x = referenceTable_full,
  delim = "\t",
  file = pathDataInterimTablesOriginalReferenceFull,
  na = ""
)

log_debug(pathDataInterimTablesOriginalStructureInchi)
readr::write_delim(
  x = structureTable_inchi,
  delim = "\t",
  file = pathDataInterimTablesOriginalStructureInchi,
  na = ""
)

log_debug(pathDataInterimTablesOriginalStructureNominal)
readr::write_delim(
  x = structureTable_nominal,
  delim = "\t",
  file = pathDataInterimTablesOriginalStructureNominal,
  na = ""
)

log_debug(pathDataInterimTablesOriginalStructureSmiles)
readr::write_delim(
  x = structureTable_smiles,
  delim = "\t",
  file = pathDataInterimTablesOriginalStructureSmiles,
  na = ""
)

log_debug(pathDataInterimTablesOriginalStructureFull)
readr::write_delim(
  x = structureTable_full,
  delim = "\t",
  file = pathDataInterimTablesOriginalStructureFull,
  na = ""
)

log_debug(pathDataInterimTablesOriginalTable)
readr::write_delim(
  x = originalTable,
  delim = "\t",
  file = pathDataInterimTablesOriginalTable,
  na = ""
)

if (mode == "full") {
  log_debug(pathTestsFile)
  readr::write_delim(
    x = originalTable_sampled,
    delim = "\t",
    file = pathTestsFile,
    na = ""
  )
}

end <- Sys.time()

log_debug("Script finished in", format(end - start))
