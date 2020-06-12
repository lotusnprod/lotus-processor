#title: "TPPT cleaneR"

#loading
##functions
source("../../functions.R")

##db
db <- "TPPT"
originalfile <- "0_initial_files/TPPT_database.xlsx"

##paths
outpath <- paste(db,
                 "_std.tsv.zip",
                 sep = "")

##files
data_original_1 <- read_excel(originalfile,
                              sheet = 1) %>%
  mutate_all(as.character)

data_original_2 <- read_excel(originalfile,
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

#joining
data_original <- left_join(data_filled, data_original_2)

#selecting
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

#standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "tpp_1",
    structure_field = c("name", "smiles")
  )

#exporting
write.table(
  x = data_standard,
  file = gzfile(description = outpath,
                compression = 9,
                encoding = "UTF-8"),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
