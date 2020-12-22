library(splitstackshape)
library(stringi)
source("r/y_as_na.R")

#' Title
#'
#' @param dfsel
#' @param dic
#'
#' @return
#' @export
#'
#' @examples
manipulating_taxo <- function(dfsel, dic) {
  # creating variables for replacement by dictionary
  a <- paste("\\b", dic$taxaRank, "\\b", sep = "")
  b <- dic$taxaRankStandard

  dfsel$rank <- gsub(
    pattern = "[.]",
    replacement = "",
    x = dfsel$rank
  )
  dfsel$rank <- stri_replace_all_regex(
    str = dfsel$rank,
    pattern = a,
    replacement = b,
    case_insensitive = FALSE,
    vectorize_all = FALSE
  )

  # removing false non-empty cells
  dfsel$rank <- y_as_na(
    x = dfsel$rank,
    y = ""
  )

  dfsel$taxonomy <- y_as_na(
    x = dfsel$taxonomy,
    y = ""
  )

  dfsel$rank <- y_as_na(
    x = dfsel$rank,
    y = ""
  )

  # selecting and splitting taxonomy and ranks
  df1 <- dfsel %>%
    select(
      organismOriginal,
      organismCleaned,
      organismDbTaxo,
      taxonId,
      # because some organisms can have multiple ids
      dbQuality,
      name = taxonomy,
      rank,
      id = ids
    ) %>%
    distinct(organismCleaned,
      organismDbTaxo,
      taxonId, # because some organisms can have multiple ids
      .keep_all = TRUE
    ) %>%
    cSplit(
      splitCols = "name",
      sep = "|"
    ) %>%
    cSplit(
      splitCols = "rank",
      sep = "|"
    ) %>%
    cSplit(
      splitCols = "id",
      sep = "|"
    ) %>%
    mutate_all(as.character) %>%
    tibble()

  # manipulating taxa
  df2 <- df1 %>%
    pivot_longer(
      cols = 6:ncol(.),
      names_to = c(".value", "level"),
      names_sep = "_",
      values_to = "taxonomy",
      values_drop_na = TRUE
    ) %>%
    distinct(organismOriginal,
      organismCleaned,
      organismDbTaxo,
      taxonId,
      # because some organisms can have multiple ids
      level,
      .keep_all = TRUE
    )

  df2$rank <- ifelse(test = is.na(df2$rank),
    yes = "NA",
    no = df2$rank
  )

  # manipulating taxa
  df3 <- df2 %>%
    filter(
      rank == "kingdom" |
        rank == "phylum" |
        rank == "class" |
        rank == "order" |
        rank == "family" |
        rank == "genus" |
        rank == "species"
    ) %>%
    pivot_wider(
      names_from = rank,
      values_from = c(name, id)
    ) %>%
    select_if(
      names(.) %in%
        c(
          "organismOriginal",
          "organismCleaned",
          "organismDbTaxo",
          "taxonId",
          "dbQuality",
          "name_kingdom",
          "name_phylum",
          "name_class",
          "name_order",
          "name_family",
          "name_genus",
          "name_species",
          "id_kingdom",
          "id_phylum",
          "id_class",
          "id_order",
          "id_family",
          "id_genus",
          "id_species"
        )
    )

  # pivoting (long)
  if (nrow(df3) != 0) {
    df4 <- df3 %>%
      pivot_longer(
        cols = 6:ncol(.),
        names_to = c(".value", "rank"),
        names_sep = "_",
        values_to = "taxonomy",
        values_drop_na = TRUE
      )

    df4$id <- sub(
      pattern = "urn:lsid:marinespecies.org:taxname:",
      replacement = "",
      x = df4$id,
      fixed = TRUE
    )
  }

  # pivoting (wide)
  if (nrow(df3) != 0) {
    df5 <- df4 %>%
      group_by(organismCleaned, organismDbTaxo, taxonId) %>%
      distinct(rank,
        name,
        id,
        .keep_all = TRUE
      ) %>%
      pivot_wider(
        names_from = rank,
        values_from = c(name, id)
      ) %>%
      ungroup() %>%
      select_if(
        names(.) %in%
          c(
            "organismCleaned",
            "organismDbTaxo",
            "taxonId",
            "dbQuality",
            "name_kingdom",
            "name_phylum",
            "name_class",
            "name_order",
            "name_family",
            "name_genus",
            "name_species",
            "id_kingdom",
            "id_phylum",
            "id_class",
            "id_order",
            "id_family",
            "id_genus",
            "id_species"
          )
      )
  }

  if (nrow(df3) != 0) {
    df5[setdiff(
      x = c(
        "organismCleaned",
        "organismDbTaxo",
        "taxonId",
        "dbQuality",
        "name_kingdom",
        "name_phylum",
        "name_class",
        "name_order",
        "name_family",
        "name_genus",
        "name_species",
        "id_kingdom",
        "id_phylum",
        "id_class",
        "id_order",
        "id_family",
        "id_genus",
        "id_species"
      ),
      y = names(df5)
    )] <- NA

    df5 <- df5 %>%
      mutate(organismCleanedNew = apply(.[, 5:11], 1, function(x) {
        tail(na.omit(x), 1)
      }))

    df5 <- df5 %>%
      mutate(taxonIdNew = apply(.[, 12:18], 1, function(x) {
        tail(na.omit(x), 1)
      }))
  }

  # adding taxa to initial df
  if (nrow(df3) != 0) {
    df6 <- left_join(dfsel, df5) %>%
      select(
        organismOriginal,
        organismDetected = organismCleaned,
        organismCleaned = organismCleanedNew,
        organismCleanedId = taxonIdNew,
        organismDbTaxo,
        organismDbTaxoQuality = dbQuality,
        organismTaxonIds = ids,
        organismTaxonRanks = rank,
        organismTaxonomy = taxonomy,
        organism_1_kingdom = name_kingdom,
        organism_2_phylum = name_phylum,
        organism_3_class = name_class,
        organism_4_order = name_order,
        organism_5_family = name_family,
        organism_6_genus = name_genus,
        organism_7_species = name_species,
        organism_1_kingdom_id = id_kingdom,
        organism_2_phylum_id = id_phylum,
        organism_3_class_id = id_class,
        organism_4_order_id = id_order,
        organism_5_family_id = id_family,
        organism_6_genus_id = id_genus,
        organism_7_species_id = id_species
      ) %>%
      filter(!is.na(organismCleaned)) %>%
      distinct()
  }

  if (nrow(df3) == 0) {
    df6 <- data.frame() %>%
      mutate(
        organismOriginal = NA,
        organismDetected = NA,
        organismCleaned = NA,
        organismCleanedId = NA,
        organismDbTaxo = NA,
        organismDbTaxoQuality = NA,
        organismTaxonIds = NA,
        organismTaxonRanks = NA,
        organismTaxonomy = NA,
        organism_1_kingdom = NA,
        organism_2_phylum = NA,
        organism_3_class = NA,
        organism_4_order = NA,
        organism_5_family = NA,
        organism_6_genus = NA,
        organism_7_species = NA,
        organism_1_kingdom_id = NA,
        organism_2_phylum_id = NA,
        organism_3_class_id = NA,
        organism_4_order_id = NA,
        organism_5_family_id = NA,
        organism_6_genus_id = NA,
        organism_7_species_id = NA
      )
  }

  return(df6)
}
