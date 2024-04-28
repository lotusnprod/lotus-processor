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
  readr::read_delim(file = wikidataLotusExporterDataOutputTaxaPath) |>
  splitstackshape::cSplit("names_pipe_separated",
    sep = "|",
    direction = "long"
  ) |>
  dplyr::distinct()

data_structures <-
  readr::read_delim(file = wikidataLotusExporterDataOutputStructuresPath)

data_structures_1 <- data_structures |>
  dplyr::filter(!grepl(pattern = "\\|", x = inchiKey))

data_structures_2 <- data_structures |>
  dplyr::filter(grepl(pattern = "\\|", x = inchiKey)) |>
  splitstackshape::cSplit(c("canonicalSmiles", "isomericSmiles", "inchi", "inchiKey"),
    sep = "|"
  ) |>
  tidyr::pivot_longer(cols = contains("_")) |>
  dplyr::filter(!is.na(value)) |>
  splitstackshape::cSplit("name",
    sep = "_"
  ) |>
  tidyr::pivot_wider(names_from = "name_1") |>
  dplyr::select(-name_2) |>
  dplyr::distinct()

data_structures <- data_structures_1 |>
  dplyr::bind_rows(data_structures_2)

data_references <-
  readr::read_delim(file = wikidataLotusExporterDataOutputReferencesPath) |>
  splitstackshape::cSplit("dois_pipe_separated",
    sep = "|",
    direction = "long"
  ) |>
  dplyr::distinct()

data_triples <-
  readr::read_delim(file = wikidataLotusExporterDataOutputTriplesPath)

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

data_organism <- data_organism |>
  dplyr::filter(!names_pipe_separated %in% names$name)

# manipulating
data_manipulated <- data_triples |>
  dplyr::inner_join(data_organism, by = c("taxon" = "wikidataId")) |>
  dplyr::inner_join(data_structures, by = c("compound" = "wikidataId")) |>
  dplyr::inner_join(data_references, by = c("reference" = "wikidataId")) |>
  dplyr::mutate(structure_smiles = if_else(
    condition = !is.na(isomericSmiles),
    true = isomericSmiles,
    false = canonicalSmiles
  )) |>
  dplyr::select(
    structure_smiles,
    # structure_inchi = inchi,
    organism_clean = names_pipe_separated,
    reference_doi = dois_pipe_separated
    # reference_title = title
  ) |>
  distinct() |>
  data.frame()

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_manipulated,
    db = "wikidata",
    structure_field = "structure_smiles",
    organism_field = "organism_clean",
    reference_field = "reference_doi"
  )

# exporting
database$writeInterim(data_standard)
