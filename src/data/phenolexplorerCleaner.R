#title: "PhenolExplorer cleaneR"

#loading
##functions
source("../../functions.R")

##db
db <- "PHENOLEXPLORER"

##paths
outpath <- paste(db,
                 "_std.tsv.zip",
                 sep = "")

#Loading all files
compounds_classification <- read_delim(
  "0_initial_files/compounds-classification.csv",
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
)

compounds_structures <- read_delim(
  "0_initial_files/compounds-structures.csv",
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
)

compounds <- read_delim(
  "0_initial_files/compounds.csv",
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
)

foods_classification <- read_delim(
  "0_initial_files/foods-classification.csv",
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
)

foods <- read_delim(
  "0_initial_files/foods.csv",
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
)

metabolites_structures <- read_delim(
  "0_initial_files/metabolites-structures.csv",
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
)

metabolites <- read_delim(
  "0_initial_files/metabolites.csv",
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
)

publications <- read_delim(
  "0_initial_files/publications.csv",
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
)

composition <- read_excel("0_initial_files/composition-data.xlsx",
                          sheet = 1)

###joining
a <- full_join(compounds, compounds_structures)
b <- full_join(a, composition, by = c("name" = "compound"))
c <- full_join(b, foods, by = c("food" = "name"))


#pivoting to join right references
colnames(c)[29] <- "publicationids"

data_pivoted <- c %>%
  cSplit("publicationids", ";") %>%
  group_by(id.x) %>%
  pivot_longer(
    38:ncol(.),
    names_to = c(".value", "level"),
    names_sep = "_",
    values_to = "taxonomy",
    values_drop_na = TRUE
  )

#adding references
data_referenced <-
  left_join(data_pivoted, publications, by = c("publicationids" = "id")) %>%
  select(
    uniqueid = id.x,
    name,
    cas = cas_number,
    pubchem = pubchem_compound_id,
    smiles,
    standardcontent = mean,
    biologicalsource = food_source_scientific_name,
    level,
    reference = title
  ) %>%
  pivot_wider(names_from = level,
              values_from = reference) %>%
  mutate_all(as.character)

data_referenced[] <-
  lapply(data_referenced, function(x)
    gsub("NULL", NA, x))

#tailing and selecting
data_selected <- data_referenced

data_selected$reference <-
  apply(data_selected[, 8:ncol(data_selected)] , 1 , paste , collapse = "|")

data_selected$reference <- gsub("|NA", "", data_selected$reference)

data_selected <- data_selected %>%
  select(uniqueid,
         name,
         cas,
         pubchem,
         smiles,
         standardcontent,
         biologicalsource,
         reference)

#standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "phe_1",
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
