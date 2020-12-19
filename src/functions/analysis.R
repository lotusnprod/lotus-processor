####################   Functions   ####################

library(tidyverse)

#######################################################

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
distinct_biosources <- function(x) {
  newdf <- x %>%
    filter(!is.na(organismLowestTaxon)) %>%
    distinct(organismLowestTaxon,
      .keep_all = TRUE
    ) %>%
    group_by(organism_7_species) %>%
    add_count() %>%
    ungroup() %>%
    filter(!is.na(organism_7_species) |
      !n > 1) %>%
    select(-n) %>%
    group_by(organism_6_genus) %>%
    add_count() %>%
    ungroup() %>%
    filter(!is.na(organism_6_genus) |
      !n > 1) %>%
    select(-n) %>%
    group_by(organism_5_family) %>%
    add_count() %>%
    ungroup() %>%
    filter(!is.na(organism_5_family) |
      !n > 1) %>%
    select(-n) %>%
    group_by(organism_4_order) %>%
    add_count() %>%
    ungroup() %>%
    filter(!is.na(organism_4_order) |
      !n > 1) %>%
    select(-n) %>%
    group_by(organism_3_class) %>%
    add_count() %>%
    ungroup() %>%
    filter(!is.na(organism_3_class) |
      !n > 1) %>%
    select(-n) %>%
    group_by(organism_2_phylum) %>%
    add_count() %>%
    ungroup() %>%
    filter(!is.na(organism_2_phylum) |
      !n > 1) %>%
    select(-n) %>%
    group_by(organism_1_kingdom) %>%
    add_count() %>%
    ungroup() %>%
    filter(!is.na(organism_1_kingdom) |
      !n > 1) %>%
    select(-n)

  return(newdf)
}

#######################################################

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
distinct_pairs <- function(x) {
  newdf <- x %>%
    filter(
      !is.na(structureCurated) &
        !is.na(organismLowestTaxon) &
        !is.na(referenceOriginal) |
        # this will have to be adapted later on
        database == "dnp_1"
    ) %>% # this will have to be adapted later on
    distinct(structureCurated,
      organismLowestTaxon,
      .keep_all = TRUE
    ) %>%
    group_by(
      structureCurated,
      organism_7_species
    ) %>%
    add_count() %>%
    ungroup() %>%
    filter(!is.na(organismLowestTaxon) |
      !n > 1) %>%
    select(-n) %>%
    group_by(
      structureCurated,
      organism_6_genus
    ) %>%
    add_count() %>%
    ungroup() %>%
    filter(!is.na(organismLowestTaxon) |
      !n > 1) %>%
    select(-n) %>%
    group_by(
      structureCurated,
      organism_5_family
    ) %>%
    add_count() %>%
    ungroup() %>%
    filter(!is.na(organismLowestTaxon) |
      !n > 1) %>%
    select(-n) %>%
    group_by(
      structureCurated,
      organism_4_order
    ) %>%
    add_count() %>%
    ungroup() %>%
    filter(!is.na(organismLowestTaxon) |
      !n > 1) %>%
    select(-n) %>%
    group_by(
      structureCurated,
      organism_3_class
    ) %>%
    add_count() %>%
    ungroup() %>%
    filter(!is.na(organismLowestTaxon) |
      !n > 1) %>%
    select(-n) %>%
    group_by(
      structureCurated,
      organism_2_phylum
    ) %>%
    add_count() %>%
    ungroup() %>%
    filter(!is.na(organismLowestTaxon) |
      !n > 1) %>%
    select(-n) %>%
    group_by(
      structureCurated,
      organism_1_kingdom
    ) %>%
    add_count() %>%
    ungroup() %>%
    filter(!is.na(organismLowestTaxon) |
      !n > 1) %>%
    select(-n)

  return(newdf)
}

#######################################################
