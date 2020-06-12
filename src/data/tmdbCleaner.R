#title: "TMDB cleaneR"

#loading
##functions

source("../../functions.R")

##db
db <- "TMDB"

##paths
originalfile <- "0_initial_files/TMDB_scraped.tsv.zip"

outpath <- paste(db,
                 "_std.tsv.zip",
                 sep = "")


##files
data_original <- read_delim(
  file = gzfile(originalfile),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)


#pivoting
data_pivoted <- data_original %>%
  mutate(level = as.numeric(gl(nrow(.) / 28, 28))) %>%
  group_by(level) %>%
  pivot_wider(names_from = 1, values_from = 2) %>%
  unnest() %>%
  ungroup()


#selecting
data_selected <- data_pivoted %>%
  select(name = `Entry name`,
         biologicalsource = `Latin name`,
         reference = References)


#standardizing
data_standard <-
  standardizing_original(data_selected = data_selected,
                         db = "tmd_1",
                         structure_field = "name")

data_standard[] <-
  lapply(data_standard, function(x)
    gsub("Not Available", NA, x))

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
