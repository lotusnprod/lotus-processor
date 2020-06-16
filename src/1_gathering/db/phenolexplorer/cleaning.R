# title: "PhenolExplorer cleaneR"

# setting working directory
setwd("~/GitLab/opennaturalproductsdb/src/")

# loading paths
source("paths.R")

# loading functions
source("functions.R")

# loading all files
compounds_classification <- read_delim(
  file = pathDataExternalDbSourcePhenolexplorerCompoundsClassification,
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
)

compounds_structures <- read_delim(
  file = pathDataExternalDbSourcePhenolexplorerCompoundsStructures,
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
)

compounds <- read_delim(
  file = pathDataExternalDbSourcePhenolexplorerCompounds,
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
)

foods_classification <- read_delim(
  file = pathDataExternalDbSourcePhenolexplorerFoodsClassification,
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
)

foods <- read_delim(
  file = pathDataExternalDbSourcePhenolexplorerFoods,
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
)

metabolites_structures <- read_delim(
  file = pathDataExternalDbSourcePhenolexplorerMetabolitesStructures,
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
)

metabolites <- read_delim(
  file = pathDataExternalDbSourcePhenolexplorerMetabolites,
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
)

publications <- read_delim(
  file = pathDataExternalDbSourcePhenolexplorerPublications,
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
)

composition <-
  read_excel(path = pathDataExternalDbSourcePhenolexplorerComposition,
             sheet = 1)

### joining
a <- full_join(compounds, compounds_structures)
b <- full_join(a, composition, by = c("name" = "compound"))
c <- full_join(b, foods, by = c("food" = "name"))


# pivoting to join right references
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

# adding references
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

# tailing and selecting
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

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "phe_1",
    structure_field = c("name", "smiles")
  )

# exporting
write.table(
  x = data_standard,
  file = gzfile(
    description = pathDataInterimDbPhenolexplorer,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
