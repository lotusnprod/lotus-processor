library(dplyr)
library(stringr)

#' Title
#'
#' @param dfsel
#'
#' @return
#' @export
#'
#' @examples
taxo_cleaning_auto <- function(dfsel) {
  cat("cleaning inconsistent taxa \n")
  df1_a <- dfsel %>%
    distinct(
      organismDbTaxo,
      organismValue,
      organismDetected,
      organism_5_family,
      organism_6_genus,
      organism_7_species,
      organism_8_variety,
      .keep_all = TRUE
    ) %>%
    group_by(organismValue, organismDetected) %>% # intra detection consistency
    add_count() %>% # see Delea versicolor example
    group_by(organismValue) %>%
    add_count(name = "m") %>%
    mutate(ratio = n / m) %>%
    ungroup() %>%
    filter(ratio > 0.34 | organismValue == organismCleaned) %>%
    select(-n, -m, -ratio) %>%
    group_by(organismDbTaxo, organism_5_family, organism_6_genus) %>% # intraDB consistency
    add_count() %>%
    group_by(organismDbTaxo, organism_6_genus) %>%
    add_count(name = "m") %>%
    mutate(ratio = n / m) %>%
    ungroup() %>%
    filter(ratio > 0.1 | is.na(organism_6_genus)) %>%
    select(-n, -m, -ratio) %>%
    group_by(organism_5_family, organism_6_genus) %>% # interDB consistency
    add_count() %>%
    group_by(organism_6_genus) %>%
    add_count(name = "m") %>%
    mutate(ratio = n / m) %>%
    ungroup() %>%
    filter(ratio > 0.1 | is.na(organism_6_genus)) %>%
    select(-n, -m, -ratio) %>%
    mutate(organismDetected_1 = word(
      string = organismDetected,
      start = 1,
      end = 1
    )) %>%
    group_by(organismDbTaxo, organismDetected_1, organism_1_kingdom) %>% # intraDB consistency complex (Echinacea example)
    add_count() %>%
    group_by(organismDbTaxo, organismDetected_1) %>%
    add_count(name = "m") %>%
    mutate(ratio = n / m) %>%
    ungroup() %>%
    filter((ratio > 0.5 |
      is.na(organism_1_kingdom)) |
      str_count(string = organismDetected, pattern = "\\b") > 2) %>%
    select(-n, -m, -ratio, -organismDetected_1) %>%
    group_by(organismDetected, organism_5_family, organism_6_genus) %>% # interDB consistency
    add_count() %>%
    group_by(organismDetected, organism_5_family) %>%
    add_count(name = "m") %>%
    mutate(ratio = n / m) %>%
    ungroup() %>%
    filter(ratio >= 0.5 | is.na(organism_6_genus)) %>%
    select(-n, -m, -ratio)

  df1_b <- dfsel %>%
    distinct(
      organismType,
      organismValue,
      organismDbTaxo,
      organismDetected
    )

  df1_c <- left_join(df1_b, df1_a)

  df2_a <- semi_join(dfsel, df1_c)

  df2_b <- anti_join(dfsel, df1_c) %>%
    distinct(
      organismType,
      organismValue,
      organismDbTaxo,
      organismDetected
    )

  df2_c <- df2_b %>%
    distinct(organismDbTaxo, organismDetected)

  df3_a <- left_join(df2_c, df1_c) %>%
    filter(!is.na(organismCleaned)) %>%
    select(-organismType, -organismValue) %>%
    distinct()

  df3_b <- left_join(df2_b, df3_a) %>%
    filter(!is.na(organismCleaned))

  df3_c <- bind_rows(df2_a, df3_b) %>%
    distinct()

  # in case you want to see the horror
  horror <- anti_join(dfsel, df3_c)

  cat("cleaning duplicate upstream taxa \n")
  # df4 <- df_sel %>%
  df4 <- df3_c %>%
    group_by(organismType, organismValue, organism_7_1_subspecies) %>%
    fill(organism_8_variety, .direction = "downup") %>%
    group_by(organismType, organismValue, organism_7_species) %>%
    fill(organism_7_1_subspecies, .direction = "downup") %>%
    group_by(organismType, organismValue, organism_6_1_subgenus) %>%
    fill(organism_7_species, .direction = "downup") %>%
    group_by(organismType, organismValue, organism_6_genus) %>%
    fill(organism_6_1_subgenus, .direction = "downup") %>%
    group_by(organismType, organismValue, organism_5_family) %>%
    fill(organism_6_genus, .direction = "downup") %>%
    group_by(organismType, organismValue, organism_4_order) %>%
    fill(organism_5_family, .direction = "downup") %>%
    group_by(organismType, organismValue, organism_3_class) %>%
    fill(organism_4_order, .direction = "downup") %>%
    group_by(organismType, organismValue, organism_2_phylum) %>%
    fill(organism_3_class, .direction = "downup") %>%
    group_by(organismType, organismValue, organism_1_kingdom) %>%
    fill(organism_2_phylum, .direction = "downup") %>%
    group_by(organismType, organismValue) %>%
    fill(organism_1_kingdom, .direction = "downup") %>%
    ungroup() %>%
    mutate(organismCleanedBis = apply(.[, grepl(
      pattern = "organism_",
      x = colnames(.),
      fixed = TRUE
    )], 1, function(x) {
      tail(na.omit(x), 1)
    }))

  if (nrow(df4) != 0) {
    df5 <- df4 %>%
      filter(name == organismCleanedBis) %>%
      distinct(
        organismType,
        organismValue,
        organismDetected,
        organismDbTaxo,
        organismDbTaxoQuality,
        # organismTaxonIds,
        organismTaxonRanks,
        organismTaxonomy,
        organismCleaned
      )
  } else {
    df5 <- tibble(
      organismType = NA,
      organismValue = NA,
      organismDetected = NA,
      organismDbTaxo = NA,
      organismDbTaxoQuality = NA,
      # organismTaxonIds=NA,
      organismTaxonRanks = NA,
      organismTaxonomy = NA,
      organismCleaned = NA,
      name = NA,
    )
  }

  # df6 <- left_join(df5, df_sel)
  df6 <- left_join(df5, df3_c) %>%
    select(-name)

  return(df6)
}
