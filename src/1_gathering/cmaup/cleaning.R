# title: "CMAUP cleaneR"

# setting working directory
setwd("~/GitLab/opennaturalproductsdb/src/")

# loading paths
source("paths.R")

# loading functions
source("functions.R")

## files
data_original_1 <- read_delim(
  file = pathDataExternalDbSourceCmaupIngredients,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE,
  col_names = FALSE
) %>%
  mutate_all(as.character)

data_original_2 <- read_delim(
  file = pathDataExternalDbSourceCmaupPlants,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

data_original_3 <- read_delim(
  file = pathDataExternalDbSourceCmaupAssociations,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE,
  col_names = FALSE
) %>%
  mutate_all(as.character)

## naming
colnames(data_original_1)[1] <- "Ingredient_ID"
colnames(data_original_1)[2] <- "pref_name"
colnames(data_original_1)[3] <- "iupac_name"
colnames(data_original_1)[4] <- "chembl_id"
colnames(data_original_1)[5] <- "pubchem_id"
colnames(data_original_1)[6] <- "zinc_id"
colnames(data_original_1)[7] <- "unpd_id"
colnames(data_original_1)[8] <- "tcmsp_id"
colnames(data_original_1)[9] <- "tcmid_id"
colnames(data_original_1)[10] <- "npass_id"
colnames(data_original_1)[11] <- "formula"
colnames(data_original_1)[12] <- "mw"
colnames(data_original_1)[13] <- "alogp"
colnames(data_original_1)[14] <- "mlogp"
colnames(data_original_1)[15] <- "xlogp"
colnames(data_original_1)[16] <- "hba"
colnames(data_original_1)[17] <- "hbd"
colnames(data_original_1)[18] <- "psa"
colnames(data_original_1)[19] <- "rotatable_bond"
colnames(data_original_1)[20] <- "rings"
colnames(data_original_1)[21] <- "heavy_atom"
colnames(data_original_1)[22] <- "lipinski_failure"
colnames(data_original_1)[23] <- "standard_inchi"
colnames(data_original_1)[24] <- "standard_inchi_key"
colnames(data_original_1)[25] <- "canonical_smiles"

colnames(data_original_3)[1] <- "Plant_ID"
colnames(data_original_3)[2] <- "Ingredient_ID"


# joining
data_original <- left_join(data_original_1, data_original_3)

data_original <- left_join(data_original, data_original_2)


# selecting
data_selected <- data_original %>%
  select(
    uniqueid = Ingredient_ID,
    name = pref_name,
    # pubchem = pubchem_id,
    smiles = canonical_smiles,
    inchi = standard_inchi,
    biologicalsource = Plant_Name
  )


# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "cma_1",
    structure_field = c("name", "inchi", "smiles")
  )

data_standard$name <- y_as_na(data_standard$name, "n.a.")

# exporting
write.table(
  x = data_standard,
  file = gzfile(
    description = pathDataInterimDbCmaup,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
