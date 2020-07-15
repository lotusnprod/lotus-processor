# title: "Open NP DB (original) compileR"

# loading paths
source("paths.R")

source("functions/helpers.R")

library(data.table)
library(dplyr)
library(stringr)

# loading files
dbs <- lapply(pathDataInterimDbDir, function(x) {
  out <- db_loader(x)
  return(out)
})

inhouseDb <- rbindlist(l = dbs, fill = TRUE)

# selecting
inhouseDbSelected <- inhouseDb %>%
  mutate(structureOriginalNominal = name) %>%
  select(
    database,
    name,
    organismOriginal = biologicalsource,
    structureOriginalInchi = inchi,
    structureOriginalNominal = name,
    structureOriginalSmiles = smiles,
    referenceOriginalAuthors = reference_authors,
    referenceOriginalDoi = reference_doi,
    referenceOriginalExternal = reference_external,
    referenceOriginalIsbn = reference_isbn,
    referenceOriginalJournal = reference_journal,
    referenceOriginalPubmed = reference_pubmed,
    referenceOriginalTitle = reference_title,
    referenceOriginalUnsplit = reference_unsplittable,
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

inhouseDbSelected$name <- y_as_na(x = inhouseDbSelected$name,
                                  y = "n.a.")

## organism
inhouseDbOrganism <- inhouseDbSelected %>%
  filter(!is.na(organismOriginal)) %>%
  distinct(organismOriginal)

## reference
### authors
inhouseDbReferenceAuthors <- inhouseDbSelected %>%
  filter(!is.na(referenceOriginalAuthors)) %>%
  distinct(referenceOriginalAuthors)

### DOI
inhouseDbReferenceDoi <- inhouseDbSelected %>%
  filter(!is.na(referenceOriginalDoi)) %>%
  distinct(referenceOriginalDoi)

row.names(inhouseDbReferenceDoi) <-
  inhouseDbReferenceDoi$referenceOriginalDoi

# ### external
# inhouseDbReferenceExternal <- inhouseDbSelected %>%
#   filter(!is.na(referenceOriginalExternal)) %>%
#   distinct(referenceOriginalExternal)

# ### ISBN
# inhouseDbReferenceIsbn <- inhouseDbSelected %>%
#   filter(!is.na(referenceOriginalIsbn)) %>%
#   distinct(referenceOriginalIsbn)

# ### journal
# inhouseDbReferenceJournal <- inhouseDbSelected %>%
#   filter(!is.na(referenceOriginalJournal)) %>%
#   distinct(referenceOriginalJournal)

### pubmed
inhouseDbReferencePubmed <- inhouseDbSelected %>%
  filter(is.na(referenceOriginalDoi)) %>%
  filter(!is.na(referenceOriginalPubmed)) %>%
  distinct(referenceOriginalPubmed)

row.names(inhouseDbReferencePubmed) <-
  inhouseDbReferencePubmed$referenceOriginalPubmed

### title
inhouseDbReferenceTitle <- inhouseDbSelected %>%
  filter(is.na(referenceOriginalDoi)) %>%
  filter(is.na(referenceOriginalPubmed)) %>%
  filter(!is.na(referenceOriginalTitle)) %>%
  distinct(referenceOriginalTitle)

### unsplit
inhouseDbReferenceUnsplit <- inhouseDbSelected %>%
  filter(is.na(referenceOriginalDoi)) %>%
  filter(is.na(referenceOriginalPubmed)) %>%
  filter(is.na(referenceOriginalTitle)) %>%
  filter(!is.na(referenceOriginalUnsplit)) %>%
  distinct(referenceOriginalUnsplit)

### unsplit
inhouseDbReferenceFull <- inhouseDbSelected %>%
  distinct(
    referenceOriginalAuthors,
    referenceOriginalDoi,
    referenceOriginalExternal,
    referenceOriginalIsbn,
    referenceOriginalJournal,
    referenceOriginalPubmed,
    referenceOriginalTitle,
    referenceOriginalUnsplit
  ) %>%
  mutate_all(as.character)

# structures
## with InChI
inhouseDbStructureInchi <- inhouseDbSelected %>%
  filter(grepl(pattern = "^InChI=.*",
               x = structureOriginalInchi)) %>%
  distinct(structureOriginalInchi)

### without InChI nor SMILES but name
inhouseDbStructureNominal <- inhouseDbSelected %>%
  filter(!grepl(pattern = "^InChI=.*",
                x = structureOriginalInchi)) %>%
  filter(is.na(structureOriginalSmiles)) %>%
  filter(!is.na(structureOriginalNominal)) %>%
  distinct(structureOriginalNominal)

### without InChI but SMILES
inhouseDbStructureSmiles <- inhouseDbSelected %>%
  filter(!grepl(pattern = "^InChI=.*",
                x = structureOriginalInchi)) %>%
  filter(!is.na(structureOriginalSmiles)) %>%
  distinct(structureOriginalSmiles)

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

# #### authors
# write.table(
#   x = inhouseDbReferenceDoi,
#   file = gzfile(
#     description = pathDataInterimTablesOriginalReferenceAuthors,
#     compression = 9,
#     encoding = "UTF-8"
#   ),
#   row.names = FALSE,
#   quote = FALSE,
#   sep = "\t",
#   fileEncoding = "UTF-8"
# )

#### DOI
write.table(
  x = inhouseDbReferenceDoi,
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

# #### external
# write.table(
#   x = inhouseDbReferenceExternal,
#   file = gzfile(
#     description = pathDataInterimTablesOriginalReferenceExternal,
#     compression = 9,
#     encoding = "UTF-8"
#   ),
#   row.names = FALSE,
#   quote = FALSE,
#   sep = "\t",
#   fileEncoding = "UTF-8"
# )

# #### ISBN
# write.table(
#   x = inhouseDbReferenceIsbn,
#   file = gzfile(
#     description = pathDataInterimTablesOriginalReferenceIsbn,
#     compression = 9,
#     encoding = "UTF-8"
#   ),
#   row.names = FALSE,
#   quote = FALSE,
#   sep = "\t",
#   fileEncoding = "UTF-8"
# )

# #### journal
# write.table(
#   x = inhouseDbReferenceJournal,
#   file = gzfile(
#     description = pathDataInterimTablesOriginalReferenceJournal,
#     compression = 9,
#     encoding = "UTF-8"
#   ),
#   row.names = FALSE,
#   quote = FALSE,
#   sep = "\t",
#   fileEncoding = "UTF-8"
# )

#### pubmed
write.table(
  x = inhouseDbReferencePubmed,
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
  x = inhouseDbReferenceTitle,
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

#### unsplit
write.table(
  x = inhouseDbReferenceUnsplit,
  file = gzfile(
    description = pathDataInterimTablesOriginalReferenceUnsplit,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = TRUE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

#### full
write.table(
  x = inhouseDbReferenceFull,
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
  x = inhouseDbStructureInchi,
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
  x = inhouseDbStructureNominal,
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
  x = inhouseDbStructureSmiles,
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

## table
write.table(
  x = inhouseDbSelected,
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
