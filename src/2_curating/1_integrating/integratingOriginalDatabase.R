# title: "Open NP DB (original) compileR"
#setwd("/home/EPGL.UNIGE.LOCAL/allardp/opennaturalproductsdb/src")
# loading paths
source("paths.R")

source("functions/helpers.R")

library(data.table)
library(dplyr)

# loading files
dbs <- lapply(pathDataInterimDbDir, function(x) {
  out <- db_loader(x)
  return(out)
})
inhouseDb <- rbindlist(l = dbs, fill = TRUE)

#head(inhouseDbSelected, 3)

# selecting
set.seed(1234)
inhouseDbSelected <- inhouseDb %>%
mutate(structureOriginalNominal = name) %>%
  select(
    database,
    name,
    organismOriginal = biologicalsource,
    structureOriginalInchi = inchi,
    structureOriginalSmiles = smiles,
    # structureOriginalNumericalCas = cas, # cas numbers problematic (see 4-Guanidinobutanoate vs 4-Guanidinobutanoic acid)
    # structureOriginalNumericalPubchem = pubchem, # pubchem problematic (see vanillic acid / vanillate)
    structureOriginalNominal = name,
    referenceOriginal = reference
  )  %>% sample_n(10000)

inhouseDbSelected$name <- y_as_na(x = inhouseDbSelected$name,
                                  y = "n.a.")

# chemical
## entries with InChI
inhouseDbStructureInchi <- inhouseDbSelected %>%
  filter(grepl(pattern = "^InChI=.*",
               x = structureOriginalInchi)) %>%
  distinct(structureOriginalInchi)

### entries without InChI but SMILES
inhouseDbStructureSmiles <- inhouseDbSelected %>%
  filter(!grepl(pattern = "^InChI=.*",
                x = structureOriginalInchi)) %>%
  filter(!is.na(structureOriginalSmiles)) %>%
  distinct(structureOriginalSmiles)

### entries without InChI nor SMILES but name
inhouseDbStructureNominal <- inhouseDbSelected %>%
  filter(!grepl(pattern = "^InChI=.*",
                x = structureOriginalInchi)) %>%
  filter(is.na(structureOriginalSmiles)) %>%
  filter(!is.na(structureOriginalNominal)) %>%
  distinct(structureOriginalNominal)

## organism
inhouseDbOrganism <- inhouseDbSelected %>%
  filter(!is.na(organismOriginal)) %>%
  distinct(organismOriginal)

## reference
inhouseDbReference <- inhouseDbSelected %>%
  filter(!is.na(referenceOriginal)) %>%
  distinct(referenceOriginal)

# exporting
## inchi
write.table(
  x = inhouseDbStructureInchi,
  file = gzfile(
    description = pathOriginalStructureInchi,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

## smiles
write.table(
  x = inhouseDbStructureSmiles,
  file = gzfile(
    description = pathOriginalStructureSmiles,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

## nominal
write.table(
  x = inhouseDbStructureNominal,
  file = gzfile(
    description = pathOriginalStructureNominal,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

## organism
## organisms for gnfinder
split_data_table(
  x = inhouseDbOrganism,
  no_rows_per_frame = 10000,
  text = "originalOrganismGnfinderUntil_",
  path_to_store = pathOriginalOrganismDistinct
)

## ref
write.table(
  x = inhouseDbReference,
  file = gzfile(
    description = pathOriginalRef,
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
    description = pathOriginalTable,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
