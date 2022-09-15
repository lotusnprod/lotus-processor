# title: "wikidata cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

library(dplyr)
library(readr)
library(splitstackshape)
library(WikidataQueryServiceR)

# get paths
database <- databases$get("wikidata")

## files
data_organism <-
  read_delim(
    file = wikidataLotusExporterDataOutputTaxaPath
  ) %>%
  cSplit("names_pipe_separated",
    sep = "|",
    direction = "long"
  ) %>%
  distinct()

data_structures <-
  read_delim(file = wikidataLotusExporterDataOutputStructuresPath) 

data_structures_1 <- data_structures %>%
  filter(!grepl(pattern = "\\|", x = inchiKey)) 

data_structures_2 <- data_structures %>%
  filter(grepl(pattern = "\\|", x = inchiKey)) %>%
  cSplit(c("canonicalSmiles", "isomericSmiles", "inchi", "inchiKey"),
         sep = "|") %>%
  pivot_longer(cols = contains("_")) %>%
  filter(!is.na(value)) %>%
  cSplit("name",
         sep = "_") %>%
  pivot_wider(names_from = "name_1") %>%
  select(-name_2) %>%
  distinct()

data_structures <- data_structures_1 %>%
  bind_rows(data_structures_2)

data_references <-
  read_delim(
    file = wikidataLotusExporterDataOutputReferencesPath
  ) %>%
  cSplit("dois_pipe_separated",
    sep = "|",
    direction = "long"
  ) %>%
  distinct()

data_triples <-
  read_delim(
    file = wikidataLotusExporterDataOutputTriplesPath
  )

# step to discard ambiguous uninomials from wikidata
names <- WikidataQueryServiceR::query_wikidata(
  sparql_query = "
SELECT DISTINCT ?name (COUNT(DISTINCT ?id) AS ?count) WHERE {
  ?id wdt:P31 wd:Q16521;
    wdt:P105 wd:Q34740;
    wdt:P225 ?name.
}
GROUP BY ?name
HAVING ((COUNT(?id)) > 1 )
ORDER BY DESC (?count)
"
)

data_organism <- data_organism %>%
  filter(!names_pipe_separated %in% names$name)

# manipulating
data_manipulated <- data_triples %>%
  inner_join(data_organism, by = c("taxon" = "wikidataId")) %>%
  inner_join(data_structures, by = c("compound" = "wikidataId")) %>%
  inner_join(data_references, by = c("reference" = "wikidataId")) %>%
  mutate(structure_smiles = if_else(
    condition = !is.na(isomericSmiles),
    true = isomericSmiles,
    false = canonicalSmiles
  )) %>%
  select(
    structure_smiles,
    structure_inchi = inchi,
    organism_clean = names_pipe_separated,
    reference_doi = dois_pipe_separated,
    reference_title = title
  ) %>%
  distinct() %>%
  data.frame()

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_manipulated,
    db = "wikidata",
    structure_field = "structure_inchi",
    organism_field = "organism_clean",
    reference_field = "reference_doi"
  )

# exporting
database$writeInterim(data_standard)
