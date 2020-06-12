#title: "SYMMAP cleaneR"

#loading
##functions
source("../../functions.R")

##db
db <- "SYMMAP"

originalfile_bio <- "0_initial_files/SymMap v1.0, SMHB file.xlsx"

originalfile_chemo <- "0_initial_files/SymMap v1.0, SMIT file.xlsx"

##paths
filenames <- list.files(path = "0_initial_files/data/",
                        pattern = "*.csv",
                        full.names = TRUE)

outpath <- paste(db,
                 "_std.tsv.zip",
                 sep = "")

##files
data_original <- do.call("rbind",
                         lapply(filenames,
                                function(x) {
                                  dat <- read.csv(x, header = TRUE, sep = ",")
                                  dat$fileName <-
                                    tools::file_path_sans_ext(basename(x))
                                  dat
                                })) %>%
  mutate_all(as.character)

data_bio <- read_excel(originalfile_bio) %>%
  mutate_all(as.character)

data_chemo <- read_excel(originalfile_chemo) %>%
  mutate_all(as.character)

#cleaning
data_original$fileName <- gsub('data-', '', data_original$fileName)

data_original$Ingredient.id <-
  gsub('SMIT', '', data_original$Ingredient.id)

data_original <- data_original %>%
  select(
    MOL_id = Ingredient.id,
    Molecule_name = Molecule.name,
    Molecule_formula = Molecule.formula,
    Molecule_weight = Molecule.weight,
    OB_score = OB.score,
    Pubchem_id = PubChem.id,
    CAS_id = CAS.id,
    Herb_id = fileName
  )

data_chemo$MOL_id <-
  str_pad (data_chemo$MOL_id, width = 5, pad = "0")

colnames(data_bio)[12] <- "TCMID_id_bio"
colnames(data_bio)[13] <- "TCM-ID_id_bio"
colnames(data_bio)[14] <- "TCMSP_id_bio"

data_full <- full_join(data_bio, data_original)

data_full <- full_join(data_chemo, data_full)

#selecting
data_selected <- data_full %>%
  select(
    uniqueid = MOL_id,
    name = Molecule_name,
    pubchem = PubChem_id,
    cas = CAS_id,
    TCMID_id,
    `TCM-ID_id`,
    TCMSP_id,
    TCMID_id_bio,
    `TCM-ID_id_bio`,
    TCMSP_id_bio,
    Latin_name,
    English_name
  ) %>%
  mutate(biologicalsource = paste(Latin_name, English_name, sep = " "))


#standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "sym_1",
    structure_field = c("name")
  )

data_standard$biologicalsource <-
  y_as_na(data_standard$biologicalsource, "NA NA")


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
