# title: "Open NP DB (sanitized) compileR"

# loading paths
source("paths.R")

# loading functions
source("functions.R")

### STOPED HERE ### (PATH UPDATE)

#writing paths
##inputs
###original
inpathTranslated <-
  "outputs/tables/1_translated/table/translatedTable.tsv.zip"

###organism
inpathOrganism <-
  "outputs/tables/2_sanitized/sanitizedOrganism.tsv.zip"

###structure
inpathStructure <-
  # "outputs/tables/2_sanitized/sanitizedStructure.tsv.zip"
  "outputs/tables/2_sanitized/sanitizedStructure.tsv"

###reference
inpathReference <-
  "outputs/tables/2_sanitized/sanitizedReference.tsv.zip"

##outputs
outpath <- "outputs/tables/2_sanitized/sanitizedTable.tsv.zip"

# loading files
## translated
dataTranslated <- read_delim(
  file = gzfile(inpathTranslated),
  col_types = cols(.default = "c"),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  data.frame()

## organism
dataSanitizedOrganism <- read_delim(
  file = gzfile(inpathOrganism),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  select(organismOriginal,
         organismTranslated,
         organismSanitized) %>%
  data.frame()

## structure
dataSanitizedStructure <- read_delim(
  file = inpathStructure,
  # file = gzfile(inpathStructure),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  select(structureTranslated,
         inchi = inchi_sanitized,) %>%
  data.frame()

## sanitized reference
dataSanitizedReference <- read_delim(
  file = gzfile(inpathReference),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  data.frame() %>%
  mutate_all(as.character)

# joining
## structure
dataSanitized <- left_join(dataTranslated, dataSanitizedStructure)

## organism
dataSanitized <- left_join(dataSanitized, dataSanitizedOrganism)

## reference
dataSanitized <- left_join(dataSanitized, dataSanitizedReference)

# exporting
write.table(
  x = dataSanitized,
  file = gzfile(
    description = outpath,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
