# title: "Metabolights cleaneR"

# loading paths
source("paths.R")
source("r/capitalize.R")
source("r/standardizing_original.R")

library(dplyr)
library(splitstackshape)
library(XML)
library(xml2)

database <- databases$get("metabolights")
xml <- XML::xmlParse(file = database$sourceFiles$xmlComplete)

ref <- XML::xpathSApply(
  doc = xml,
  path = "//ref",
  fun = xmlAttrs,
  simplify = TRUE
)

cross_references <- data.frame(t(ref))

field <- XML::xpathSApply(
  doc = xml,
  path = "//field",
  fun = xmlValue,
  simplify = FALSE
)

field_name <- XML::xpathSApply(
  doc = xml,
  path = "//field",
  fun = xmlAttrs,
  simplify = FALSE
)

additional_fields <- rbind(field)
additional_fields <- data.frame(t(additional_fields))

field_names <- rbind(field_name)
field_names <- data.frame(t(field_names))

entry <- XML::xpathSApply(
  doc = xml,
  path = "//entry",
  fun = xmlAttrs,
  simplify = FALSE
)

entries <- rbind(entry)
entries <- data.frame(t(entries)) |>
  dplyr::mutate_all(as.character)

cross_references_l <-
  xml2::xml_find_all(
    x = xml2::read_xml(x = database$sourceFiles$xmlComplete),
    xpath = ".//entry/cross_references"
  )

additional_fields_l <-
  xml2::xml_find_all(
    x = xml2::read_xml(x = database$sourceFiles$xmlComplete),
    xpath = ".//entry/additional_fields"
  )

X <- seq_len(nrow(entries))

get_length_cross_references <- function(X) {
  length <- length(cross_references_l[[X]] |>
    as_list())

  return(length)
}

get_length_additional_fields <- function(X) {
  length <- length(additional_fields_l[[X]] |>
    as_list())

  return(length)
}

get_length <- function(.function) {
  df1 <- invisible(lapply(
    FUN = .function,
    X = X
  ))

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
  cbind(id_cross_references, cross_references) |>
  dplyr::mutate_all(as.character)

additional_fields_organism <- cbind(
  id_additional_fields,
  additional_fields,
  field_names
) |>
  dplyr::filter(field_name == "Organism") |>
  dplyr::select(id,
    organism = field
  ) |>
  dplyr::mutate_all(as.character)

additional_fields_structure <- cbind(
  id_additional_fields,
  additional_fields,
  field_names
) |>
  dplyr::filter(field_name == "inchi") |>
  dplyr::select(id,
    inchi = field
  ) |>
  dplyr::mutate_all(as.character)

global_df <- dplyr::left_join(entries, cross_references_full) |>
  dplyr::left_join(entries, by = c("dbkey" = "entry")) |>
  dplyr::left_join(additional_fields_organism, by = c("id.x" = "id")) |>
  dplyr::left_join(additional_fields_structure, by = c("id.y" = "id")) |>
  dplyr::select(
    id_study = entry,
    pubmed_ref = dbkey,
    organism,
    inchi
  ) |>
  dplyr::filter((!is.na(organism) &
    !is.na(inchi)) | !is.na(as.integer(pubmed_ref)))

ref_df <- global_df |>
  dplyr::filter(!is.na(as.integer(pubmed_ref))) |>
  dplyr::distinct(id_study, pubmed_ref)

organism_structure_df <- global_df |>
  dplyr::filter(!is.na(organism) & !is.na(inchi)) |>
  dplyr::distinct(id_study, organism, inchi)

final_df <- dplyr::left_join(ref_df, organism_structure_df) |>
  dplyr::filter(!is.na(organism) &
    !is.na(inchi) &
    !is.na(pubmed_ref)) |>
  dplyr::mutate(organism = gsub(
    pattern = "NCBITAXON:",
    replacement = "",
    x = organism,
    fixed = TRUE
  )) |>
  dplyr::mutate(organism = gsub(
    pattern = "EFO:",
    replacement = "",
    x = organism,
    fixed = TRUE
  )) |>
  dplyr::distinct() |>
  splitstackshape::cSplit("organism", sep = ";", direction = "long") |>
  dplyr::filter(!grepl(
    pattern = "blank",
    x = organism,
    ignore.case = TRUE
  )) |>
  dplyr::filter(!grepl(
    pattern = "reference",
    x = organism,
    ignore.case = TRUE
  )) |>
  dplyr::filter(!grepl(
    pattern = "standard",
    x = organism,
    ignore.case = TRUE
  )) |>
  dplyr::filter(!grepl(
    pattern = "solvent",
    x = organism,
    ignore.case = TRUE
  )) |>
  dplyr::filter(!grepl(
    pattern = "culture",
    x = organism,
    ignore.case = TRUE
  )) |>
  dplyr::filter(!grepl(
    pattern = "buffer",
    x = organism,
    ignore.case = TRUE
  )) |>
  dplyr::filter(!grepl(
    pattern = " feed",
    x = organism,
    ignore.case = TRUE
  )) |>
  dplyr::filter(!grepl(
    pattern = "soil",
    ## indicates two organisms
    x = organism,
    ignore.case = TRUE
  )) |>
  dplyr::filter(!grepl(
    pattern = "matrix",
    ## indicates two organisms
    x = organism,
    ignore.case = TRUE
  )) |>
  dplyr::filter(!grepl(
    pattern = "control",
    ## indicates two organisms
    x = organism,
    ignore.case = TRUE
  )) |>
  dplyr::filter(!grepl(
    pattern = "/",
    ## indicates two organisms
    x = organism,
    ignore.case = TRUE
  )) |>
  dplyr::mutate(organism = gsub(
    pattern = "[",
    replacement = "",
    x = organism,
    fixed = TRUE
  )) |>
  dplyr::mutate(organism = gsub(
    pattern = "]",
    replacement = "",
    x = organism,
    fixed = TRUE
  )) |>
  dplyr::mutate(organism = tolower(organism)) |>
  dplyr::mutate(organism = capitalize(organism)) |>
  dplyr::filter(pubmed_ref != "32331455") |>
  dplyr::distinct()

data_selected <- final_df |>
  dplyr::select(
    structure_inchi = inchi,
    organism_clean = organism,
    reference_pubmed = pubmed_ref
  ) |>
  distinct() |>
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
