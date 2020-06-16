# title: "TPPT cleaneR"

# loading paths
source("paths.R")

# loading functions
source("functions.R")

## files
data_original_1 <- read_excel(pathDataExternalDbSourceTpptOriginal,
                              sheet = 1) %>%
  mutate_all(as.character)

data_original_2 <- read_excel(pathDataExternalDbSourceTpptOriginal,
                              sheet = 3) %>%
  mutate_all(as.character)

data_filled <- data_original_1 %>%
  mutate(smiles = ifelse(
    Stereo_SMILES == "NI",
    Canonical_SMILES,
    ifelse(
      Stereo_SMILES == "NS",
      Canonical_SMILES,
      ifelse(Stereo_SMILES == "racemat",
             Canonical_SMILES,
             Stereo_SMILES)
    )
  ))

# joining
data_original <- left_join(data_filled, data_original_2)

# selecting
data_selected <- data_original %>%
  select(
    Phytotoxin_number,
    name = Phytotoxin_name,
    CASRN,
    smiles,
    PubChem_CID,
    biologicalsource = Latin_plant_name,
    reference = References
  ) %>%
  mutate(reference = gsub(",", "|", reference))

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "tpp_1",
    structure_field = c("name", "smiles")
  )

# exporting
write.table(
  x = data_standard,
  file = gzfile(
    description = pathDataInterimDbTmmc,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
