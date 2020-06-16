# title: "TMMC cleaneR"

# loading paths
source("paths.R")

# loading functions
source("functions.R")

## files
data_original <- read_excel(pathDataExternalDbSourceTmmcOriginal,
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

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_original_long,
    db = "tmm_1",
    structure_field = c("name")
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
