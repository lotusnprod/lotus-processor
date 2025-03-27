source("r/log_debug.R")
log_debug("This script outputs some metrics related to the DB")

start <- Sys.time()

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")

log_debug("... libraries")
library(data.table)
library(dplyr)
library(readr)
library(splitstackshape)

log_debug("loading ...")
log_debug("databases list ...")
dataset <- readr::read_delim(
  file = "../docs/dataset.csv",
  delim = ","
) |>
  dplyr::select(
    -`initial retrieved unique entries`,
    -`cleaned referenced structure-organism pairs`,
    -`pairs validated for wikidata export`,
    -`actual pairs on wikidata`,
    -`timestamp`
  ) |>
  data.frame(check.names = FALSE)

log_debug("initial table ...")
dbTable <- lapply(
  pathDataInterimDbDir,
  readr::read_delim,
  delim = "\t",
  col_types = cols(.default = "c"),
  locale = locales
) |>
  data.table::rbindlist(fill = TRUE) |>
  dplyr::filter(database != "custom") |>
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
  )

log_debug("final table ...")
inhouseDbMinimal <-
  readr::read_delim(
    file = pathDataInterimTablesCuratedTable,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales,
    col_select = c(
      "database",
      "organismCleaned",
      "structureCleanedInchi",
      "referenceCleanedDoi"
    )
  ) |>
  dplyr::filter(database != "custom") |>
  dplyr::distinct()

log_debug("validated for export ...")
openDb <-
  readr::read_delim(
    file = pathDataInterimTablesAnalyzedPlatinum,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales,
    col_select = c(
      "database",
      "organismCleaned",
      "organismCleaned_dbTaxoTaxonRanks",
      "structureCleanedInchi",
      "structureCleaned_inchikey2D",
      "referenceCleanedDoi"
    )
  ) |>
  dplyr::filter(database != "custom") |>
  dplyr::distinct()

log_debug("exported ...")
data_organism <-
  readr::read_delim(
    file = wikidataLotusExporterDataOutputTaxaPath,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales,
    col_select = c(
      "organism_wikidata" = "wikidataId",
      "organismCleaned" = "names_pipe_separated"
    )
  ) |>
  splitstackshape::cSplit("organismCleaned", sep = "|", direction = "long") |>
  dplyr::distinct()

data_structures <-
  readr::read_delim(
    file = wikidataLotusExporterDataOutputStructuresPath,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales,
    col_select = c(
      "structure_wikidata" = "wikidataId",
      "structureCleanedInchiKey" = "inchiKey",
      "structureCleanedInchi" = "inchi",
      "canonicalSmiles",
      "isomericSmiles"
    )
  ) |>
  dplyr::distinct() |>
  dplyr::mutate(
    structureValue = if_else(
      condition = !is.na(isomericSmiles),
      true = isomericSmiles,
      false = canonicalSmiles
    )
  )

data_structures_translation <-
  readr::read_delim(
    file = pathDataInterimDictionariesStructureDictionary
  )

data_references <-
  readr::read_delim(
    file = wikidataLotusExporterDataOutputReferencesPath,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales,
    col_select = c(
      "reference_wikidata" = "wikidataId",
      "referenceCleanedDoi" = "dois_pipe_separated",
      "referenceCleanedTitle" = "title",
    )
  ) |>
  splitstackshape::cSplit(
    "referenceCleanedDoi",
    sep = "|",
    direction = "long"
  ) |>
  dplyr::distinct()

wikidata_pairs <-
  readr::read_delim(
    file = pathLastWdExport,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales,
    col_select = c(
      "structureValue" = "structure_smiles",
      "organismCleaned" = "organism_clean",
      "referenceCleanedDoi" = "reference_doi"
    )
  ) |>
  filter(
    !is.na(structureValue) &
      !is.na(organismCleaned) &
      !is.na(referenceCleanedDoi)
  ) |>
  dplyr::left_join(data_structures_translation) |>
  dplyr::left_join(data_organism) |>
  dplyr::left_join(data_structures) |>
  dplyr::left_join(data_references) |>
  dplyr::distinct()

log_debug("... closed db")
closedDb <-
  readr::read_delim(
    file = file.path(pathDataInterimTablesAnalyzed, "closed.tsv.gz"),
    col_types = cols(.default = "c"),
    locale = locales,
    col_select = c(
      "database",
      "organismCleaned",
      "organismCleaned_dbTaxoTaxonRanks",
      "structureCleanedInchi",
      "structureCleaned_inchikey2D",
      "referenceCleanedDoi"
    )
  )

log_debug("... generating timestamps")
system(command = paste("bash", pathTimestampsScript))

timestamps <- readr::read_delim(
  file = pathDataInterimTimestamps,
  delim = "\t",
  locale = locales,
  col_names = FALSE
)

dates <- timestamps$X1 |>
  matrix(ncol = 2, byrow = TRUE) |>
  data.frame() |>
  dplyr::mutate(
    database = gsub(
      pattern = "../data/external/dbSource/",
      replacement = "",
      x = X1,
      fixed = TRUE
    ),
    timestamp = as.character.POSIXt(X2)
  ) |>
  dplyr::mutate(
    database = gsub(
      pattern = "/.*",
      replacement = "",
      x = database
    )
  ) |>
  dplyr::arrange(dplyr::desc(timestamp)) |>
  dplyr::distinct(database, .keep_all = TRUE) |>
  dplyr::select(database, timestamp)

log_debug("performing inner join with uploaded entries")
openDb_u_wiki <- openDb |>
  dplyr::inner_join(wikidata_pairs)

inhouseDb <- dplyr::bind_rows(closedDb, openDb)

initial_stats <- dbTable |>
  dplyr::group_by(database) |>
  dplyr::count(name = "initial retrieved unique entries")

cleaned_stats_3D <- inhouseDbMinimal |>
  dplyr::group_by(database) |>
  dplyr::distinct(
    organismCleaned,
    structureCleanedInchi,
    referenceCleanedDoi
  ) |>
  dplyr::count(name = "cleaned referenced structure-organism pairs")

final_stats_3D <- openDb |>
  dplyr::group_by(database) |>
  dplyr::distinct(
    organismCleaned,
    structureCleanedInchi,
    structureCleaned_inchikey2D,
    referenceCleanedDoi
  ) |>
  dplyr::count(name = "pairs validated for wikidata export")

wiki_stats_3D <- openDb_u_wiki |>
  dplyr::group_by(database) |>
  dplyr::distinct(
    organismCleaned,
    structureCleanedInchi,
    structureCleaned_inchikey2D,
    referenceCleanedDoi
  ) |>
  dplyr::count(name = "actual pairs on wikidata")

stats_table <- dplyr::left_join(initial_stats, cleaned_stats_3D) |>
  dplyr::left_join(final_stats_3D) |>
  dplyr::left_join(wiki_stats_3D)

dataset <- dataset %>%
  dplyr::full_join(dates) %>%
  dplyr::full_join(stats_table) %>%
  dplyr::select(
    database,
    `type`,
    `initial retrieved unique entries`,
    `cleaned referenced structure-organism pairs`,
    `pairs validated for wikidata export`,
    `actual pairs on wikidata`,
    dplyr::everything()
  ) %>%
  dplyr::relocate(timestamp, .after = status) %>%
  dplyr::mutate_all(~ replace(., is.na(.), "-"))

pairsOpenDb_3D <- openDb |>
  dplyr::filter(
    !is.na(organismCleaned) &
      !is.na(structureCleanedInchi)
  ) |>
  dplyr::distinct(structureCleanedInchi, organismCleaned, .keep_all = TRUE)

pairsOpenDb_2D <- openDb |>
  dplyr::filter(
    !is.na(organismCleaned) &
      !is.na(structureCleaned_inchikey2D)
  ) |>
  dplyr::distinct(
    structureCleaned_inchikey2D,
    organismCleaned,
    .keep_all = TRUE
  )

"%ni%" <- Negate("%in%")

pairsOutsideClosed_3D <- inhouseDb |>
  dplyr::filter(
    !is.na(organismCleaned) &
      !is.na(structureCleanedInchi)
  ) |>
  dplyr::distinct(structureCleanedInchi, organismCleaned, .keep_all = TRUE) |>
  dplyr::filter(database %ni% forbidden_export)

pairsOutsideClosed_2D <- inhouseDb |>
  dplyr::filter(
    !is.na(organismCleaned) &
      !is.na(structureCleaned_inchikey2D)
  ) |>
  dplyr::distinct(
    structureCleaned_inchikey2D,
    organismCleaned,
    .keep_all = TRUE
  ) |>
  dplyr::filter(database %ni% forbidden_export)

pairsFull_3D <- inhouseDb |>
  dplyr::filter(
    !is.na(organismCleaned) &
      !is.na(structureCleanedInchi)
  ) |>
  dplyr::distinct(structureCleanedInchi, organismCleaned, .keep_all = TRUE)

pairsFull_2D <- inhouseDb |>
  dplyr::filter(
    !is.na(organismCleaned) &
      !is.na(structureCleaned_inchikey2D)
  ) |>
  dplyr::distinct(
    structureCleaned_inchikey2D,
    organismCleaned,
    .keep_all = TRUE
  )

pairsClosed_3D <- closedDb |>
  dplyr::distinct(structureCleanedInchi, organismCleaned, .keep_all = TRUE)

pairsClosed_2D <- closedDb |>
  dplyr::distinct(
    structureCleaned_inchikey2D,
    organismCleaned,
    .keep_all = TRUE
  )

stats_3D <- pairsOutsideClosed_3D |>
  dplyr::group_by(database) |>
  dplyr::count() |>
  dplyr::arrange(dplyr::desc(n))

stats_2D <- pairsOutsideClosed_2D |>
  dplyr::group_by(database) |>
  dplyr::count() |>
  dplyr::arrange(dplyr::desc(n))

# unique
log_debug("analyzing unique organisms per db")
## biological taxa
### open NP DB
openDbOrganism <- openDb |>
  dplyr::filter(!is.na(organismCleaned)) |>
  dplyr::distinct(organismCleaned)

### inhouseDB
inhouseDbOrganism <- inhouseDb |>
  dplyr::filter(!is.na(organismCleaned)) |>
  dplyr::distinct(organismCleaned)

log_debug(paste(
  "inhouse:",
  nrow(inhouseDbOrganism),
  "distinct organisms",
  sep = " "
))

### Closed
closedDbOrganism <- closedDb |>
  dplyr::filter(!is.na(organismCleaned)) |>
  dplyr::distinct(organismCleaned)

log_debug(paste(
  "closed:",
  nrow(closedDbOrganism),
  "distinct organisms",
  sep = " "
))

## structures
log_debug("analyzing unique structures (2D) per db")
### open NP DB
openDbStructure_3D <- openDb |>
  dplyr::filter(!is.na(structureCleanedInchi)) |>
  dplyr::distinct(structureCleanedInchi, .keep_all = TRUE)

### inhouseDB
inhouseDbStructure_3D <- inhouseDb |>
  dplyr::filter(!is.na(structureCleanedInchi)) |>
  dplyr::distinct(structureCleanedInchi, .keep_all = TRUE)

log_debug(paste(
  "inhouse:",
  nrow(inhouseDbStructure_3D),
  "distinct structures (3D)",
  sep = " "
))

### open NP DB
openDbStructure_2D <- openDb |>
  dplyr::filter(!is.na(structureCleaned_inchikey2D)) |>
  dplyr::distinct(structureCleaned_inchikey2D, .keep_all = TRUE)

### inhouseDB
inhouseDbStructure_2D <- inhouseDb |>
  dplyr::filter(!is.na(structureCleaned_inchikey2D)) |>
  dplyr::distinct(structureCleaned_inchikey2D, .keep_all = TRUE)

log_debug(paste(
  "inhouse:",
  nrow(inhouseDbStructure_2D),
  "distinct structures (2D)",
  sep = " "
))

### closed
closedDbStructure_3D <- closedDb |>
  dplyr::filter(!is.na(structureCleanedInchi)) |>
  dplyr::distinct(structureCleanedInchi, .keep_all = TRUE)

log_debug(paste(
  "closed:",
  nrow(closedDbStructure_3D),
  "distinct structures (3D)",
  sep = " "
))

closedDbStructure_2D <- closedDb |>
  dplyr::filter(!is.na(structureCleaned_inchikey2D)) |>
  dplyr::distinct(structureCleaned_inchikey2D, .keep_all = TRUE)

log_debug(paste(
  "closed:",
  nrow(closedDbStructure_2D),
  "distinct structures (2D)",
  sep = " "
))

structuresPerOrganism_3D <- pairsOpenDb_3D |>
  dplyr::filter(grepl(
    pattern = "species",
    x = organismCleaned_dbTaxoTaxonRanks
  )) |>
  dplyr::filter(grepl(pattern = " ", x = organismCleaned)) |>
  dplyr::distinct(organismCleaned, structureCleanedInchi) |>
  dplyr::group_by(organismCleaned) |>
  dplyr::count()

structuresPerOrganism_2D <- pairsOpenDb_2D |>
  dplyr::filter(grepl(
    pattern = "species",
    x = organismCleaned_dbTaxoTaxonRanks
  )) |>
  dplyr::filter(grepl(pattern = " ", x = organismCleaned)) |>
  dplyr::distinct(organismCleaned, structureCleaned_inchikey2D) |>
  dplyr::group_by(organismCleaned) |>
  dplyr::count()

tableStructures_3D <-
  c(
    "only001_structure" = sum(structuresPerOrganism_3D$n == 1),
    "between001and010_structures" = sum(
      structuresPerOrganism_3D$n > 1 &
        structuresPerOrganism_3D$n <= 10
    ),
    "between010and100_structures" = sum(
      structuresPerOrganism_3D$n > 10 &
        structuresPerOrganism_3D$n <= 100
    ),
    "above100_structures" = sum(structuresPerOrganism_3D$n > 100)
  ) |>
  data.frame()

colnames(tableStructures_3D)[1] <- "organisms"

tableStructures_2D <-
  c(
    "only001_structure" = sum(structuresPerOrganism_2D$n == 1),
    "between001and010_structures" = sum(
      structuresPerOrganism_2D$n > 1 &
        structuresPerOrganism_2D$n <= 10
    ),
    "between010and100_structures" = sum(
      structuresPerOrganism_2D$n > 10 &
        structuresPerOrganism_2D$n <= 100
    ),
    "above100_structures" = sum(structuresPerOrganism_2D$n > 100)
  ) |>
  data.frame()

colnames(tableStructures_2D)[1] <- "organisms"

organismsPerStructure_3D <- pairsOpenDb_3D |>
  dplyr::filter(grepl(
    pattern = "species",
    x = organismCleaned_dbTaxoTaxonRanks
  )) |>
  dplyr::filter(grepl(pattern = " ", x = organismCleaned)) |>
  dplyr::distinct(organismCleaned, structureCleanedInchi) |>
  dplyr::group_by(structureCleanedInchi) |>
  dplyr::count()

organismsPerStructure_2D <- pairsOpenDb_2D %>%
  dplyr::filter(grepl(
    pattern = "species",
    x = organismCleaned_dbTaxoTaxonRanks
  )) |>
  dplyr::filter(grepl(pattern = " ", x = organismCleaned)) |>
  dplyr::distinct(organismCleaned, structureCleaned_inchikey2D) |>
  dplyr::group_by(structureCleaned_inchikey2D) |>
  dplyr::count()

tableOrganisms_3D <-
  c(
    "only001_organism" = sum(organismsPerStructure_3D$n == 1),
    "between001and010_organisms" = sum(
      organismsPerStructure_3D$n > 1 &
        organismsPerStructure_3D$n <= 10
    ),
    "between010and100_organisms" = sum(
      organismsPerStructure_3D$n > 10 &
        organismsPerStructure_3D$n <= 100
    ),
    "above100_organisms" = sum(organismsPerStructure_3D$n > 100)
  ) |>
  data.frame()

colnames(tableOrganisms_3D)[1] <- "structures"

tableOrganisms_2D <-
  c(
    "only001_organism" = sum(organismsPerStructure_2D$n == 1),
    "between001and010_organisms" = sum(
      organismsPerStructure_2D$n > 1 &
        organismsPerStructure_2D$n <= 10
    ),
    "between010and100_organisms" = sum(
      organismsPerStructure_2D$n > 10 &
        organismsPerStructure_2D$n <= 100
    ),
    "above100_organisms" = sum(organismsPerStructure_2D$n > 100)
  ) |>
  data.frame()

colnames(tableOrganisms_2D)[1] <- "structures"

readr::write_delim(
  x = dataset,
  delim = ",",
  file = "../docs/dataset.csv"
)

cat(
  paste(
    nrow(pairsOpenDb_3D),
    "|",
    nrow(pairsOpenDb_2D),
    "(3D|2D)",
    "unique referenced structure-organism pairs. \n",
    "They consist of \n",
    nrow(openDbStructure_3D),
    "|",
    nrow(openDbStructure_2D),
    "(3D|2D)",
    "unique curated structures and \n",
    nrow(openDbOrganism),
    "unique organisms,\n",
    "originating from \n",
    nrow(
      pairsOpenDb_2D |>
        dplyr::distinct(database)
    ),
    "initial open databases. \n",
    "\n",
    "Among 2D structures, \n",
    tableOrganisms_2D[1, 1],
    "are present in only 1 organism, \n",
    tableOrganisms_2D[2, 1],
    "are present in between 1 and 10 organisms, \n",
    tableOrganisms_2D[3, 1],
    "are present in between 10 and 100 organisms, \n",
    tableOrganisms_2D[4, 1],
    "are present in more than 100 organisms. \n",
    "\n",
    "Among organisms, \n",
    tableStructures_2D[1, 1],
    "contain only 1 2D structure, \n",
    tableStructures_2D[2, 1],
    "contain between 1 and 10 2D structures, \n",
    tableStructures_2D[3, 1],
    "contain between 10 and 100 2D structures, \n",
    tableStructures_2D[4, 1],
    "contain more than 100 2D structures. \n",
    sep = " "
  ),
  file = "../docs/metrics.md"
)

create_dir(pathDataProcessed)
readr::write_delim(
  x = structuresPerOrganism_2D,
  delim = "\t",
  file = file.path(
    pathDataProcessed,
    "structures2DPerOrganism.tsv.gz"
  )
)

readr::write_delim(
  x = organismsPerStructure_2D,
  delim = "\t",
  file = file.path(
    pathDataProcessed,
    "organismsPerStructure2D.tsv.gz"
  )
)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
