# title: "Nubbe cleaneR"

# loading paths
source("paths.R")
source("functions/helpers.R")
source("functions/parallel.R")
source("functions/standardizing.R")

library(dplyr)
library(pbmcapply)
library(parallel)
library(data.table)
library(splitstackshape) # provides cSplit
library(rvest)  # provides read_html
library(tidyr) #provides pivot_wider
library(XML)

# get paths
database <- databases$get("nubbe")

## files
ids <- read_delim(
  file = database$sourceFiles$tsv,
  delim = "\t",
  escape_double = FALSE,
  col_names = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

X <- ids$X1

list <- list()

for (i in 1:length(X)) {
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
      test = xml_data$especies == "\n\t",
      yes = NA,
      no = (ifelse(
        test = is.null(xml_data$especies$origem$familia),
        yes = NA,
        no = (xml_data$especies$origem$familia)
      ))
    )
  organismGenus <-
    ifelse(
      test = xml_data$especies == "\n\t",
      yes = NA,
      no = (ifelse(
        test = is.null(xml_data$especies$origem$genero),
        yes = NA,
        no = (xml_data$especies$origem$genero)
      ))
    )
  organismSpecies <- ifelse(
    test = xml_data$especies == "\n\t",
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
    biologicalsource = paste(organismFamily,
                             organismGenus,
                             organismSpecies,
                             sep = " "),
    inchi = paste("InChI=", inchi, sep = "")
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
      pattern = "NA ",
      replacement = "",
      x = biologicalsource,
      fixed = TRUE
    )
  ) %>%
  select(
    id,
    name,
    inchi,
    inchikey,
    smiles,
    #temporary
    biologicalsource,
    reference_authors = referenceAuthors,
    reference_title = referenceTitle,
    reference_doi
  ) %>%
  data.frame()

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_manipulated,
    db = "nub_1",
    structure_field = c("inchi", "name", "smiles"),
    reference_field = c("reference_authors", "reference_title", "reference_doi")
  )

# exporting
database$writeInterim(data_standard)
