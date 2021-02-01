source("paths.R")
library(tidyverse)

#' Title
#'
#' @param data_selected
#' @param db
#' @param organism_field
#' @param structure_field
#' @param reference_field
#'
#' @return
#' @export
#'
#' @examples
standardizing_original <- function(data_selected,
                                   db,
                                   organism_field,
                                   # possibilities: c("organism_clean","organism_dirty")
                                   structure_field,
                                   # possibilities: c("structure_inchi","structure_smiles","structure_name")
                                   reference_field)
# possibilities: c("reference_authors","reference_doi","reference_external","reference_isbn","reference_journal", "reference_original","reference_publishingDetails","reference_pubmed","reference_split","reference_title")
{
  data_selected[setdiff(
    c(
      "structure_name",
      "structure_inchi",
      "structure_smiles",
      "organism_clean",
      "organism_dirty",
      "reference_authors",
      "reference_doi",
      "reference_external",
      "reference_isbn",
      "reference_journal",
      "reference_original",
      "reference_publishingDetails",
      "reference_pubmed",
      "reference_split",
      "reference_title"
    ),
    names(data_selected)
  )] <- NA

  data_standard <- data.frame(data_selected) %>%
    mutate(database = db) %>%
    select(
      database,
      all_of(structure_field),
      all_of(organism_field),
      all_of(reference_field)
    ) %>%
    filter_at(vars(all_of(structure_field)), any_vars(!is.na(.))) %>%
    filter_at(vars(all_of(structure_field)), any_vars(grepl(pattern = "[[:alnum:]]", x = .))) %>%
    filter_at(vars(all_of(organism_field)), any_vars(!is.na(.))) %>%
    filter_at(vars(all_of(organism_field)), any_vars(grepl(pattern = "[[:alpha:]]", x = .))) %>%
    filter_at(vars(all_of(reference_field)), any_vars(!is.na(.))) %>%
    filter_at(vars(all_of(reference_field)), any_vars(grepl(pattern = "[[:alnum:]]", x = .))) %>%
    distinct_at(vars(
      all_of(structure_field),
      all_of(organism_field),
      all_of(reference_field)
    ),
    .keep_all = TRUE
    )

  data_standard[] <-
    lapply(data_standard, function(x) {
      gsub("\r\n", " ", x)
    })
  data_standard[] <-
    lapply(data_standard, function(x) {
      gsub("\r", " ", x)
    })
  data_standard[] <-
    lapply(data_standard, function(x) {
      gsub("\n", " ", x)
    })
  data_standard[] <-
    lapply(data_standard, function(x) {
      gsub("\t", " ", x)
    })

  data_standard <- data_standard %>%
    mutate_all(~ iconv(x = ., from = "utf-8", to = "utf-8//ignore")) %>%
    mutate_all(
      .tbl = .,
      .funs = trimws
    )

  return(data_standard)
}
