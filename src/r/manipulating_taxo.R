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

  # dfsel$ids <- y_as_na(
  #   x = dfsel$ids,
  #   y = ""
  # )

  # selecting and splitting taxonomy and ranks
  df1 <- dfsel %>%
    select(
      organismCleaned,
      organismDbTaxo,
      taxonId,
      # because some organisms can have multiple ids
      # dbQuality,
      name = taxonomy,
      rank,
      # id = ids
    ) %>%
    distinct(
      organismCleaned,
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
    # cSplit(
    #   splitCols = "id",
    #   sep = "|"
    # ) %>%
    mutate_all(as.character) %>%
    tibble()

  # manipulating taxa
  df2 <- df1 %>%
    pivot_longer(
      cols = 4:ncol(.),
      names_to = c(".value", "level"),
      names_sep = "_",
      values_to = "taxonomy",
      values_drop_na = TRUE
    ) %>%
    distinct(
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

  # df2$id <- sub(
  #   pattern = "urn:lsid:marinespecies.org:taxname:",
  #   replacement = "",
  #   x = df2$id,
  #   fixed = TRUE
  # )

  df2$taxonId <- sub(
    pattern = "urn:lsid:marinespecies.org:taxname:",
    replacement = "",
    x = df2$taxonId,
    fixed = TRUE
  )

  df2_a <- df2 %>%
    filter(!is.na(rank)) %>%
    filter(rank != "") %>%
    arrange(desc(level)) %>%
    distinct(
      organismCleaned,
      organismDbTaxo,
      taxonId,
      .keep_all = TRUE
    ) %>%
    select(
      taxon_rank = rank,
      everything(),
      -level
    ) %>%
    distinct()

  # manipulating taxa
  df3 <- df2 %>%
    filter(
      rank == "kingdom" |
        rank == "phylum" |
        rank == "class" |
        rank == "order" |
        # rank == "infraorder" |
        rank == "family" |
        # rank == "subfamily" |
        # rank == "tribe" |
        # rank == "subtribe" |
        rank == "genus" |
        rank == "subgenus" |
        rank == "species" |
        rank == "subspecies" |
        rank == "variety"
    ) %>%
    pivot_wider(
      names_from = rank,
      values_from = name
    ) %>%
    select(-level)

  # pivoting (long)
  if (nrow(df3) != 0) {
    df4 <- df3 %>%
      pivot_longer(
        cols = 4:ncol(.),
        names_to = "rank",
        values_to = "name",
        values_drop_na = TRUE
      )
  }

  # pivoting (wide)
  if (nrow(df3) != 0) {
    df5 <- df4 %>%
      group_by(organismCleaned, organismDbTaxo, taxonId) %>%
      distinct(rank,
        name,
        .keep_all = TRUE
      ) %>%
      pivot_wider(
        names_from = rank,
        values_from = name
      ) %>%
      ungroup() %>%
      select_if(
        names(.) %in%
          c(
            "organismCleaned",
            "organismDbTaxo",
            "taxonId",
            "dbQuality",
            "kingdom",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "subgnus",
            "species",
            "subspecies",
            "variety"
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
        "kingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "subgenus",
        "species",
        "subspecies",
        "variety"
        # "id_kingdom",
        # "id_phylum",
        # "id_class",
        # "id_order",
        # "id_family",
        # "id_genus",
        # "id_species",
        # "id_variety"
      ),
      y = names(df5)
    )] <- NA
  }

  # adding taxa to initial df
  if (nrow(df3) != 0) {
    df6 <- left_join(dfsel, df5) %>%
      left_join(., df2_a) %>%
      select(
        organismOriginal,
        organismDetected = organismCleaned,
        organismCleaned = currentCanonicalFull,
        organismCleanedId = taxonId,
        organismCleanedRank = taxon_rank,
        organismDbTaxo,
        organismDbTaxoQuality = dbQuality,
        name,
        # organismTaxonIds = ids,
        organismTaxonRanks = rank,
        organismTaxonomy = taxonomy,
        organism_1_kingdom = kingdom,
        organism_2_phylum = phylum,
        organism_3_class = class,
        organism_4_order = order,
        organism_5_family = family,
        organism_6_genus = genus,
        organism_6_1_subgenus = subgenus,
        organism_7_species = species,
        organism_7_1_subspecies = subspecies,
        organism_8_variety = variety,
        # organism_1_kingdom_id = id_kingdom,
        # organism_2_phylum_id = id_phylum,
        # organism_3_class_id = id_class,
        # organism_4_order_id = id_order,
        # organism_5_family_id = id_family,
        # organism_6_genus_id = id_genus,
        # organism_7_species_id = id_species,
        # organism_8_variety_id = id_variety
      ) %>%
      filter(!is.na(organismCleaned)) %>%
      filter(!is.na(organismCleanedRank)) %>%
      filter(organismCleanedRank != "NA") %>%
      distinct()
  }

  if (nrow(df3) == 0) {
    df6 <- data.frame() %>%
      mutate(
        organismOriginal = NA,
        organismDetected = NA,
        organismCleaned = NA,
        organismCleanedId = NA,
        organismCleanedRank = NA,
        organismDbTaxo = NA,
        organismDbTaxoQuality = NA,
        # organismTaxonIds = NA,
        organismTaxonRanks = NA,
        organismTaxonomy = NA,
        organism_1_kingdom = NA,
        organism_2_phylum = NA,
        organism_3_class = NA,
        organism_4_order = NA,
        organism_5_family = NA,
        organism_6_genus = NA,
        organism_6_1_subgenus = NA,
        organism_7_species = NA,
        organism_7_1_subspecies = NA,
        organism_8_variety = NA,
        # organism_1_kingdom_id = NA,
        # organism_2_phylum_id = NA,
        # organism_3_class_id = NA,
        # organism_4_order_id = NA,
        # organism_5_family_id = NA,
        # organism_6_genus_id = NA,
        # organism_7_species_id = NA,
        # organism_8_variety_id = NA
      )
  }

  return(df6)
}
