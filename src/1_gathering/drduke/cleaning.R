# title: "DrDuke cleaneR"

# setting working directory
setwd("~/GitLab/opennaturalproductsdb/src/")

# loading paths
source("paths.R")

# loading functions
source("functions.R")

## files
data_common <- read_delim(
  file = pathDataExternalDbSourceDrdukeCommonNames,
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character) %>%
  select(FNFNUM, CNNAM)

data_farmacy <- read_delim(
  file = pathDataExternalDbSourceDrdukeFarmacy,
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

data_fntax <- read_delim(
  file = pathDataExternalDbSourceDrdukeTaxa,
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character) %>%
  select(FNFNUM, TAXON)

data_reference <- read_delim(
  file = pathDataExternalDbSourceDrdukeReferences,
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character) %>%
  select(REFERENCE, LONGREF)

# joining
data_joined <- left_join(data_farmacy, data_fntax)

data_joined <- left_join(data_joined, data_reference)

# selecting
data_selected <- data_joined %>%
  select(name = CHEM,
         biologicalsource = TAXON,
         reference = LONGREF) %>%
  mutate(reference = ifelse(
    test = !is.na(reference),
    yes = reference,
    no = "DrDuke"
  ))

# standardizing
data_standard <-
  standardizing_original(data_selected = data_selected,
                         db = "duk_1",
                         structure_field = "name")

# exporting
write.table(
  x = data_standard,
  file = gzfile(
    description = pathDataInterimDbDrduke,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
