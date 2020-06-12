#title: "MITISHAMBA cleaneR"

#loading
##functions
source("../../functions.R")

##db
db <- "MITISHAMBA"
originalfile <- "0_initial_files/Mitishamba_db_scraped.tsv.zip"

##paths
outpath <- paste(db,
                 "_std.tsv.zip",
                 sep = "")

##files
data_original <- read_delim(
  file = gzfile(originalfile),
  delim = "\t",
  trim_ws = TRUE,
  escape_backslash = TRUE
) %>%
  mutate_all(as.character)

#selecting
data_selected <- data_original %>%
  select(
    smiles,
    biologicalsource = plant_species,
    name = common_name,
    reference = authors
  )

#standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "mit_1",
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
