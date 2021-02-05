# title: "Metabolights cleaneR"

# loading paths
source("paths.R")
source("r/parallel.R")
source("r/standardizing_original.R")

library(dplyr)
library(parallel)
library(Hmisc)
library(pbmcapply)
library(splitstackshape)
library(tidyverse)
library(vroom)
library(XML)
library(xml2)

database <- databases$get("metabolights")
xml <- xmlParse(database$sourceFiles$xmlComplete)

ref <- xpathSApply(
  doc = xml,
  path = "//ref",
  fun = xmlAttrs,
  simplify = TRUE
)

cross_references <- data.frame(t(ref))

field <- xpathSApply(
  doc = xml,
  path = "//field",
  fun = xmlValue,
  simplify = FALSE
)

field_name <- xpathSApply(
  doc = xml,
  path = "//field",
  fun = xmlAttrs,
  simplify = FALSE
)

additional_fields <- rbind(field)
additional_fields <- data.frame(t(additional_fields))

field_names <- rbind(field_name)
field_names <- data.frame(t(field_names))

entry <- xpathSApply(
  doc = xml,
  path = "//entry",
  fun = xmlAttrs,
  simplify = FALSE
)

entries <- rbind(entry)
entries <- data.frame(t(entries)) %>%
  mutate_all(as.character)

cross_references_l <-
  xml_find_all(
    x = read_xml(database$sourceFiles$xmlComplete),
    xpath = ".//entry/cross_references"
  )

additional_fields_l <-
  xml_find_all(
    x = read_xml(database$sourceFiles$xmlComplete),
    xpath = ".//entry/additional_fields"
  )

X <- seq_len(nrow(entries))

get_length_cross_references <- function(X) {
  length <- length(cross_references_l[[X]] %>% as_list())

  return(length)
}

get_length_additional_fields <- function(X) {
  length <- length(additional_fields_l[[X]] %>% as_list())

  return(length)
}

get_length <- function(.function) {
  df1 <- invisible(
    pbmclapply(
      FUN = .function,
      X = X,
      mc.preschedule = TRUE,
      mc.set.seed = TRUE,
      mc.silent = TRUE,
      mc.cores = (parallel::detectCores() - 2),
      mc.cleanup = TRUE,
      mc.allow.recursive = TRUE,
      ignore.interactive = TRUE
    )
  )

  df2 <- rbind(df1)
  df3 <- data.frame(t(df2))
  df3$id <- rownames(df3)

  df4 <- df3[rep(seq(nrow(df3)), df3$df1), ] %>%
    select(id)
  return(df4)
}

entries$id <- row.names(entries)

id_cross_references <-
  get_length(.function = get_length_cross_references)

id_additional_fields <-
  get_length(.function = get_length_additional_fields)

cross_references_full <-
  cbind(id_cross_references, cross_references) %>%
  mutate_all(as.character)

additional_fields_organism <- cbind(
  id_additional_fields,
  additional_fields,
  field_names
) %>%
  filter(field_name == "Organism") %>%
  select(id,
    organism = field
  ) %>%
  mutate_all(as.character)

additional_fields_structure <- cbind(
  id_additional_fields,
  additional_fields,
  field_names
) %>%
  filter(field_name == "inchi") %>%
  select(id,
    inchi = field
  ) %>%
  mutate_all(as.character)

global_df <- left_join(entries, cross_references_full) %>%
  left_join(., entries, by = c("dbkey" = "entry")) %>%
  left_join(., additional_fields_organism, by = c("id.x" = "id")) %>%
  left_join(., additional_fields_structure, by = c("id.y" = "id")) %>%
  select(
    id_study = entry,
    pubmed_ref = dbkey,
    organism,
    inchi
  ) %>%
  filter((!is.na(organism) &
    !is.na(inchi)) | !is.na(as.integer(pubmed_ref)))

ref_df <- global_df %>%
  filter(!is.na(as.integer(pubmed_ref))) %>%
  distinct(id_study, pubmed_ref)

organism_structure_df <- global_df %>%
  filter(!is.na(organism) & !is.na(inchi)) %>%
  distinct(id_study, organism, inchi)

final_df <- left_join(ref_df, organism_structure_df) %>%
  filter(!is.na(organism) &
    !is.na(inchi) &
    !is.na(pubmed_ref)) %>%
  mutate(organism = gsub(
    pattern = "NCBITAXON:",
    replacement = "",
    x = organism,
    fixed = TRUE
  )) %>%
  mutate(organism = gsub(
    pattern = "EFO:",
    replacement = "",
    x = organism,
    fixed = TRUE
  )) %>%
  distinct() %>%
  cSplit("organism", sep = ";", direction = "long") %>%
  filter(!grepl(
    pattern = "blank",
    x = organism,
    ignore.case = TRUE
  )) %>%
  filter(!grepl(
    pattern = "reference",
    x = organism,
    ignore.case = TRUE
  )) %>%
  filter(!grepl(
    pattern = "standard",
    x = organism,
    ignore.case = TRUE
  )) %>%
  filter(!grepl(
    pattern = "solvent",
    x = organism,
    ignore.case = TRUE
  )) %>%
  filter(!grepl(
    pattern = "culture",
    x = organism,
    ignore.case = TRUE
  )) %>%
  filter(!grepl(
    pattern = "buffer",
    x = organism,
    ignore.case = TRUE
  )) %>%
  filter(!grepl(
    pattern = " feed",
    x = organism,
    ignore.case = TRUE
  )) %>%
  filter(!grepl(
    pattern = "soil",
    ## indicates two organisms
    x = organism,
    ignore.case = TRUE
  )) %>%
  filter(!grepl(
    pattern = "matrix",
    ## indicates two organisms
    x = organism,
    ignore.case = TRUE
  )) %>%
  filter(!grepl(
    pattern = "control",
    ## indicates two organisms
    x = organism,
    ignore.case = TRUE
  )) %>%
  filter(!grepl(
    pattern = "/",
    ## indicates two organisms
    x = organism,
    ignore.case = TRUE
  )) %>%
  mutate(organism = gsub(
    pattern = "[",
    replacement = "",
    x = organism,
    fixed = TRUE
  )) %>%
  mutate(organism = gsub(
    pattern = "]",
    replacement = "",
    x = organism,
    fixed = TRUE
  )) %>%
  mutate(organism = tolower(organism)) %>%
  mutate(organism = capitalize(organism)) %>%
  filter(pubmed_ref != "32331455") %>%
  distinct()

data_selected <- final_df %>%
  select(
    structure_inchi = inchi,
    organism_clean = organism,
    reference_pubmed = pubmed_ref
  ) %>%
  distinct() %>%
  data.frame()

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "metabolights",
    structure_field = "structure_inchi",
    organism_field = "organism_clean",
    reference_field = "reference_pubmed"
  )

# exporting
database$writeInterim(data_standard)
