# title: "Nubbe cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

library(dplyr)
library(splitstackshape) # provides cSplit
library(XML)

# get paths
database <- databases$get("nubbe")

## files
X <- database$sourceFiles$tsv

list <- list()

for (i in seq_along(X)) {
  data <- xmlParse(file = X[i])

  xml_data <- xmlToList(node = data)

  id <- (xml_data$id)
  name <- (xml_data$nome)
  inchi <- (xml_data$inchi)
  inchikey <- (xml_data$inchikey)
  smiles <- (xml_data$smiles)
  referenceYear <-
    ifelse(
      test = is.null(xml_data$publicacoes$artigo),
      yes = xml_data$publicacoes$capitulo_libro$ano,
      no = xml_data$publicacoes$artigo$ano
    )
  referenceAuthors <-
    ifelse(
      test = is.null(xml_data$publicacoes$artigo),
      yes = xml_data$publicacoes$capitulo_libro$autores,
      no = xml_data$publicacoes$artigo$autores
    )
  referenceTitle <-
    ifelse(
      test = is.null(xml_data$publicacoes$artigo),
      yes = xml_data$publicacoes$capitulo_libro$titulo,
      no = xml_data$publicacoes$artigo$titulo
    )
  referenceFull <-
    ifelse(
      test = is.null(xml_data$publicacoes$artigo),
      yes = xml_data$publicacoes$capitulo_libro$compilado,
      no = xml_data$publicacoes$artigo$compilado
    )
  organismFamily <-
    ifelse(
      test = xml_data$especies == "\n    ",
      yes = NA,
      no = (ifelse(
        test = is.null(xml_data$especies$origem$familia),
        yes = NA,
        no = (xml_data$especies$origem$familia)
      ))
    )
  organismGenus <-
    ifelse(
      test = xml_data$especies == "\n    ",
      yes = NA,
      no = (ifelse(
        test = is.null(xml_data$especies$origem$genero),
        yes = NA,
        no = (xml_data$especies$origem$genero)
      ))
    )
  organismSpecies <- ifelse(
    test = xml_data$especies == "\n    ",
    yes = NA,
    no = (ifelse(
      test = is.null(xml_data$especies$origem$especie),
      yes = NA,
      no = (xml_data$especies$origem$especie)
    ))
  )

  df <- data.frame(
    id,
    name,
    inchi,
    inchikey,
    smiles,
    referenceYear,
    referenceAuthors,
    referenceTitle,
    referenceFull,
    organismFamily,
    organismGenus,
    organismSpecies
  )

  list[[i]] <- df
}

# manipulating
data_original <- bind_rows(list)

data_manipulated <- data_original %>%
  mutate(
    biologicalsource = paste(organismGenus,
      organismSpecies,
      organismFamily,
      sep = " "
    ),
    inchi = paste0("InChI=", inchi)
  ) %>%
  cSplit("referenceFull", sep = ";") %>%
  mutate_all(as.character) %>%
  mutate(
    reference_doi = gsub(
      pattern = "DOI: ",
      replacement = "",
      x = referenceFull_01
    ),
    biologicalsource = gsub(
      pattern = " NA",
      replacement = "",
      x = biologicalsource,
      fixed = TRUE
    )
  ) %>%
  select(
    id,
    structure_name = name,
    structure_inchi = inchi,
    structure_inchikey = inchikey,
    structure_smiles = smiles,
    organism_clean = biologicalsource,
    reference_authors = referenceAuthors,
    reference_title = referenceTitle,
    reference_doi
  ) %>%
  data.frame()

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_manipulated,
    db = "nubbe",
    structure_field = "structure_smiles",
    organism_field = "organism_clean",
    reference_field = c("reference_authors", "reference_title", "reference_doi")
  )

# exporting
database$writeInterim(data_standard)
