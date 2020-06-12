#title: "NPASS cleaneR"

#loading
##functions
source("../../functions.R")

##db
db <- "NPASS"
originalfile_1 <-
  "0_initial_files/NPASSv1.0_download_naturalProducts_generalInfo.txt"
originalfile_2 <-
  "0_initial_files/NPASSv1.0_download_naturalProducts_properties.txt"
originalfile_3 <-
  "0_initial_files/NPASSv1.0_download_naturalProducts_speciesInfo.txt"
originalfile_4 <-
  "0_initial_files/NPASSv1.0_download_naturalProducts_species_pair.txt"

##paths
outpath <- paste(db,
                 "_std.tsv.zip",
                 sep = "")

##files
data_original_1 <- read_delim(
  file = originalfile_1,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

data_original_2 <- read_delim(
  file = originalfile_2,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

data_original_3 <- read_delim(
  file = originalfile_3,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

data_original_4 <- read_delim(
  file = originalfile_4,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

#joining
data_original <- left_join(data_original_1, data_original_2)

data_original <- left_join(data_original, data_original_4)

data_original <- left_join(data_original, data_original_3)

#selecting
data_selected <- data_original %>%
  mutate(reference = paste(ref_id, ref_id_type, sep = "§")) %>%
  select(
    pubchem = pubchem_cid,
    np_id,
    name = pref_name,
    inchi = standard_inchi,
    standard_inchi_key,
    smiles = canonical_smiles,
    biologicalsource = org_name,
    reference
  )

#standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "npa_1",
    structure_field = c("name", "inchi", "smiles")
  )

data_standard$name <- y_as_na(data_standard$name, "n.a.")

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
