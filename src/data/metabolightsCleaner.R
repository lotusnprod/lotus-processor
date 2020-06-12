#title: "Metabolights cleaneR"

#loading
##functions
source("../../functions.R")

inpathData <- "0_initial_files/metabolights_data_std.tsv.zip"

inpathStudies <- "0_initial_files/metabolights_studies_std.tsv.zip"

outpath <- "METABOLIGHTS_std.tsv.zip"

data_clean_final <- read_delim(
  file = gzfile(inpathData),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

species_studies <- read_delim(
  file = gzfile(inpathStudies),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  distinct(species)

species <- capitalize(species_studies$species)

data_selected <- data_clean_final %>%
  filter(!biologicalsource %in% species)

data_selected$inchi <- y_as_na(data_selected$inchi, "NULL")
data_selected$name <- y_as_na(data_selected$name, "NULL")
data_selected$biologicalsource <-
  y_as_na(data_selected$biologicalsource, "NULL")
data_selected$biologicalsource <-
  y_as_na(data_selected$biologicalsource, "reference compound")

data_selected <- data_selected %>%
  filter(!is.na(biologicalsource))

#standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "met_1",
    structure_field = c("name", "inchi")
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
