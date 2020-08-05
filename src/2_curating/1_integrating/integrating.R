# title: "Open NP DB (original) compileR"

# loading paths
source("paths.R")

source("functions/helpers.R")

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

# loading files
dbs <- lapply(pathDataInterimDbDir, function(x) {
  out <- db_loader(x)
  return(out)
})

inhouseDb <- rbindlist(l = dbs, fill = TRUE)

# selecting
inhouseDbSelected <- inhouseDb %>%
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
  )

if (mode == "min")
  set.seed(42)

if (mode == "min")
  inhouseDbSelected <- inhouseDbSelected %>%
  sample_n(2000)

inhouseDbSelected[] <-
  lapply(inhouseDbSelected, function(x)
    gsub("\r\n", " ", x))
inhouseDbSelected[] <-
  lapply(inhouseDbSelected, function(x)
    gsub("\r", " ", x))
inhouseDbSelected[] <-
  lapply(inhouseDbSelected, function(x)
    gsub("\n", " ", x))
inhouseDbSelected[] <-
  lapply(inhouseDbSelected, function(x)
    gsub("\t", " ", x))

## organism
inhouseDbOrganism <- inhouseDbSelected %>%
  filter(!is.na(organismOriginal)) %>%
  distinct(organismOriginal)

## reference
### DOI
inhouseDbReference_doi <- inhouseDbSelected %>%
  filter(!is.na(referenceOriginal_doi)) %>%
  distinct(referenceOriginal_doi)

row.names(inhouseDbReference_doi) <-
  inhouseDbReference_doi$referenceOriginal_doi

### pubmed
inhouseDbReference_pubmed <- inhouseDbSelected %>%
  filter(!is.na(referenceOriginal_pubmed)) %>%
  distinct(referenceOriginal_pubmed)

row.names(inhouseDbReference_pubmed) <-
  inhouseDbReference_pubmed$referenceOriginal_pubmed

### title
inhouseDbReference_title <- inhouseDbSelected %>%
  filter(is.na(referenceOriginal_doi)) %>%
  filter(is.na(referenceOriginal_pubmed)) %>%
  filter(!is.na(referenceOriginal_title)) %>%
  distinct(referenceOriginal_title)

### split
inhouseDbReference_publishingDetails <- inhouseDbSelected %>%
  filter(is.na(referenceOriginal_doi)) %>%
  filter(is.na(referenceOriginal_pubmed)) %>%
  filter(is.na(referenceOriginal_title)) %>%
  filter(!is.na(referenceOriginal_publishingDetails)) %>%
  distinct(referenceOriginal_publishingDetails)

### split
inhouseDbReference_split <- inhouseDbSelected %>%
  filter(is.na(referenceOriginal_doi)) %>%
  filter(is.na(referenceOriginal_pubmed)) %>%
  filter(is.na(referenceOriginal_title)) %>%
  filter(is.na(referenceOriginal_publishingDetails)) %>%
  filter(!is.na(referenceOriginal_split)) %>%
  distinct(referenceOriginal_split)

### original
inhouseDbReference_original <- inhouseDbSelected %>%
  filter(is.na(referenceOriginal_doi)) %>%
  filter(is.na(referenceOriginal_pubmed)) %>%
  filter(is.na(referenceOriginal_title)) %>%
  filter(is.na(referenceOriginal_publishingDetails)) %>%
  filter(is.na(referenceOriginal_split)) %>%
  filter(!is.na(referenceOriginal_original)) %>%
  distinct(referenceOriginal_original)

### full
inhouseDbReference_full <- inhouseDbSelected %>%
  distinct(
    organismOriginal,
    referenceOriginal_authors,
    referenceOriginal_doi,
    referenceOriginal_external,
    referenceOriginal_isbn,
    referenceOriginal_journal,
    referenceOriginal_original,
    referenceOriginal_pubmed,
    referenceOriginal_title,
    referenceOriginal_split
  ) %>%
  mutate_all(as.character) %>%
  pivot_longer(
    cols = 2:ncol(.),
    names_to = c("drop", "referenceType"),
    names_sep = "_",
    values_to = "referenceValue",
    values_drop_na = TRUE
  ) %>%
  select(organismOriginal,
         referenceType,
         referenceValue)

# structures
## with InChI
inhouseDbStructure_inchi <- inhouseDbSelected %>%
  filter(grepl(pattern = "^InChI=.*",
               x = structureOriginal_inchi)) %>%
  distinct(structureOriginal_inchi)

### without InChI but SMILES
inhouseDbStructure_smiles <- inhouseDbSelected %>%
  filter(!grepl(pattern = "^InChI=.*",
                x = structureOriginal_inchi)) %>%
  filter(!is.na(structureOriginal_smiles)) %>%
  distinct(structureOriginal_smiles)

### without InChI nor SMILES but name
inhouseDbStructure_nominal <- inhouseDbSelected %>%
  filter(!grepl(pattern = "^InChI=.*",
                x = structureOriginal_inchi)) %>%
  filter(is.na(structureOriginal_smiles)) %>%
  filter(!is.na(structureOriginal_nominal)) %>%
  distinct(structureOriginal_nominal)

### full
inhouseDbStructure_full <- inhouseDbSelected %>%
  distinct(structureOriginal_inchi,
           structureOriginal_smiles,
           structureOriginal_nominal) %>%
  mutate_all(as.character) %>%
  pivot_longer(
    cols = 1:ncol(.),
    names_to = c("drop", "structureType"),
    names_sep = "_",
    values_to = "structureValue",
    values_drop_na = TRUE
  ) %>%
  select(structureType, structureValue)


## full table
originalTable <- inhouseDbSelected %>%
  select(database, organismOriginal, everything()) %>%
  pivot_longer(
    cols = 6:ncol(.),
    names_to = c("drop", "referenceType"),
    names_sep = "_",
    values_to = "referenceValue",
    values_drop_na = TRUE
  ) %>%
  pivot_longer(
    cols = 3:5,
    names_to = c("drop2", "structureType"),
    names_sep = "_",
    values_to = "structureValue",
    values_drop_na = TRUE
  ) %>%
  select(-drop, -drop2)

# exporting
## creating directories if they do not exist
### interim tables
ifelse(!dir.exists(pathDataInterimTables),
       dir.create(pathDataInterimTables),
       FALSE)

#### original
ifelse(
  !dir.exists(pathDataInterimTablesOriginal),
  dir.create(pathDataInterimTablesOriginal),
  FALSE
)

##### organism
ifelse(
  !dir.exists(pathDataInterimTablesOriginalOrganism),
  dir.create(pathDataInterimTablesOriginalOrganism),
  FALSE
)

##### reference
ifelse(
  !dir.exists(pathDataInterimTablesOriginalReference),
  dir.create(pathDataInterimTablesOriginalReference),
  FALSE
)

ifelse(
  !dir.exists(pathDataInterimTablesOriginalReferenceOriginalFolder),
  dir.create(pathDataInterimTablesOriginalReferenceOriginalFolder),
  FALSE
)

##### structure
ifelse(
  !dir.exists(pathDataInterimTablesOriginalStructure),
  dir.create(pathDataInterimTablesOriginalStructure),
  FALSE
)

## writing files
### organism
split_data_table(
  x = inhouseDbOrganism,
  no_rows_per_frame = 10000,
  text = "",
  path_to_store = pathDataInterimTablesOriginalOrganism
)

### reference
#### DOI
write.table(
  x = inhouseDbReference_doi,
  file = gzfile(
    description = pathDataInterimTablesOriginalReferenceDoi,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = TRUE,
  quote = TRUE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

#### pubmed
write.table(
  x = inhouseDbReference_pubmed,
  file = gzfile(
    description = pathDataInterimTablesOriginalReferencePubmed,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = TRUE,
  quote = TRUE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

#### title
write.table(
  x = inhouseDbReference_title,
  file = gzfile(
    description = pathDataInterimTablesOriginalReferenceTitle,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = TRUE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

#### title
write.table(
  x = inhouseDbReference_publishingDetails,
  file = gzfile(
    description = pathDataInterimTablesOriginalReferencePublishingDetails,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = TRUE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

#### split
write.table(
  x = inhouseDbReference_split,
  file = gzfile(
    description = pathDataInterimTablesOriginalReferenceSplit,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = TRUE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

#### original
split_data_table(
  x = inhouseDbReference_original,
  no_rows_per_frame = 1000,
  text = "",
  path_to_store = pathDataInterimTablesOriginalReferenceOriginalFolder
)

#### full
write.table(
  x = inhouseDbReference_full,
  file = gzfile(
    description = pathDataInterimTablesOriginalReferenceFull,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = TRUE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

### structure
#### inchi
write.table(
  x = inhouseDbStructure_inchi,
  file = gzfile(
    description = pathDataInterimTablesOriginalStructureInchi,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

#### nominal
write.table(
  x = inhouseDbStructure_nominal,
  file = gzfile(
    description = pathDataInterimTablesOriginalStructureNominal,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

#### smiles
write.table(
  x = inhouseDbStructure_smiles,
  file = gzfile(
    description = pathDataInterimTablesOriginalStructureSmiles,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

#### full
write.table(
  x = inhouseDbStructure_full,
  file = gzfile(
    description = pathDataInterimTablesOriginalStructureFull,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = TRUE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

## table
write.table(
  x = originalTable,
  file = gzfile(
    description = pathDataInterimTablesOriginalTable,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
