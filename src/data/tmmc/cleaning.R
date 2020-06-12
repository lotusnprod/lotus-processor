#title: "TMMC cleaneR"

#loading
##functions
source("../../functions.R")

##db
db <- "TMMC"
originalfile <- "0_initial_files/compound.xlsx"

##paths
outpath <- paste(db,
                 "_std.tsv.zip",
                 sep = "")

##files
data_original <- read_excel(originalfile,
                            sheet = 1) %>%
  mutate_all(as.character)

data_original_long <- data_original %>%
  cSplit("CSID", "|") %>%
  pivot_longer(17:ncol(.)) %>%
  filter(!is.na(value)) %>%
  distinct(SCIENCE, COMPOUND, .keep_all = TRUE) %>%
  mutate(
    name = COMPOUND,
    biologicalsource = str_extract(SCIENCE, "(?<=\\[).+?(?=\\])"),
    reference = LINK
  ) %>%
  distinct(name, .keep_all = TRUE) %>%
  mutate(
    biologicalsource = gsub("<i>", "", biologicalsource),
    biologicalsource = gsub("</i>", "", biologicalsource)
  )

#standardizing
data_standard <-
  standardizing_original(
    data_selected = data_original_long,
    db = "tmm_1",
    structure_field = c("name")
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
