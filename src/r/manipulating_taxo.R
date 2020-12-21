library(splitstackshape)
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
      taxonomy,
      rank,
      ids
    ) %>%
    distinct(organismCleaned,
      organismDbTaxo,
      taxonId, # because some organisms can have multiple ids
      .keep_all = TRUE
    ) %>%
    cSplit(
      splitCols = "taxonomy",
      sep = "|"
    ) %>%
    cSplit(
      splitCols = "rank",
      sep = "|"
    ) %>%
    mutate_all(as.character) %>%
    tibble()

  # manipulating taxa
  df2 <- df1 %>%
    pivot_longer(
      cols = 7:ncol(.),
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
        rank == "species" |
        rank == "variety"
    ) %>%
    pivot_wider(
      names_from = rank,
      values_from = taxonomy
    ) %>%
    select_if(
      names(.) %in%
        c(
          "organismOriginal",
          "organismCleaned",
          "organismDbTaxo",
          "taxonId",
          "ids",
          "dbQuality",
          "kingdom",
          "phylum",
          "class",
          "order",
          "family",
          "genus",
          "species",
          "variety"
        )
    )

  # pasting suffix to colnames to pivot then (the double pivot allows to tidy the data)
  colnames(df3)[7:ncol(df3)] <-
    paste("bio_", colnames(df3)[7:ncol(df3)], sep = "")

  # pivoting (long)
  if (nrow(df3) != 0) {
    df4 <- df3 %>%
      pivot_longer(
        cols = 7:ncol(.),
        names_to = c(".value", "level"),
        names_sep = "_",
        values_to = "taxonomy",
        values_drop_na = TRUE
      )
  }

  # pivoting (wide)
  if (nrow(df3) != 0) {
    df5 <- df4 %>%
      group_by(organismCleaned, organismDbTaxo, taxonId) %>%
      distinct(ids,
        level,
        .keep_all = TRUE
      ) %>%
      pivot_wider(
        names_from = level,
        values_from = bio
      ) %>%
      select_if(
        names(.) %in%
          c(
            "organismCleaned",
            "organismDbTaxo",
            "ids",
            "dbQuality",
            "kingdom",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species",
            "variety"
          )
      )
  }

  if (nrow(df3) != 0) {
    df5[setdiff(
      x = c(
        "organismCleaned",
        "organismDbTaxo",
        "ids",
        "dbQuality",
        "kingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
        "variety"
      ),
      y = names(df5)
    )] <- NA
  }

  if (nrow(df3) != 0) {
    df5 <- df5 %>%
      select(
        organismCleaned,
        organismDbTaxo,
        organismDbTaxoQuality = dbQuality,
        organismTaxonId = ids,
        organism_1_kingdom = kingdom,
        organism_2_phylum = phylum,
        organism_3_class = class,
        organism_4_order = order,
        organism_5_family = family,
        organism_6_genus = genus,
        organism_7_species = species,
        organism_8_variety = variety
      )
  }

  # adding taxa to initial df
  if (nrow(df3) != 0) {
    df6 <- left_join(dfsel, df5) %>%
      select(
        organismOriginal,
        organismCleaned,
        organismDbTaxo,
        organismDbTaxoQuality = dbQuality,
        organismTaxonIds = ids,
        organismTaxonRanks = rank,
        organismTaxonomy = taxonomy,
        organism_1_kingdom,
        organism_2_phylum,
        organism_3_class,
        organism_4_order,
        organism_5_family,
        organism_6_genus,
        organism_7_species,
        organism_8_variety
      )
  }

  if (nrow(df3) == 0) {
    df6 <- data.frame() %>%
      mutate(
        organismOriginal = NA,
        organismCleaned = NA,
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
        organism_8_variety = NA
      )
  }

  return(df6)
}
