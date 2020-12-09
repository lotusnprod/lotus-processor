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
      taxonId, # because some organisms can have multiple ids
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
      taxonId, # because some organisms can have multiple ids
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

taxo_cleaning_manual <- function(dfsel) {
  inhouse_db <- dfsel

  inhouse_db <- inhouse_db %>%
    mutate_at(
      .vars = vars(
        organism_1_kingdom,
        organism_2_phylum,
        organism_3_class,
        organism_4_order,
        organism_5_family,
        organism_6_genus,
        organism_7_species,
        organism_8_variety
      ),
      .funs = function(x) {
        gsub(
          pattern = "^c\\(|\\)$",
          replacement = "",
          x = x
        )
      }
    ) %>%
    mutate_at(
      .vars = vars(
        organism_1_kingdom,
        organism_2_phylum,
        organism_3_class,
        organism_4_order,
        organism_5_family,
        organism_6_genus,
        organism_7_species,
        organism_8_variety
      ),
      .funs = function(x) {
        gsub(
          pattern = "\"",
          replacement = "",
          x = x
        )
      }
    )

  inhouse_db$organism_1_kingdom <- gsub(
    pattern = "Viridiplantae",
    replacement = "Plantae",
    x = inhouse_db$organism_1_kingdom
  )

  inhouse_db$organism_1_kingdom <- gsub(
    pattern = "Metazoa",
    replacement = "Animalia",
    x = inhouse_db$organism_1_kingdom
  )

  inhouse_db$organism_2_phylum <- gsub(
    pattern = "Streptophyta",
    replacement = "Tracheophyta",
    x = inhouse_db$organism_2_phylum
  )

  inhouse_db$organism_5_family <- gsub(
    pattern = "Trichocomaceae",
    replacement = "Aspergillaceae",
    x = inhouse_db$organism_5_family
  )

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Allomyrina dichotoma"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Allomyrina dichotoma"] <-
    "Animalia"
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Allomyrina dichotoma"] <-
    "Arthropoda"
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Allomyrina dichotoma"] <-
    "Insecta"
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Allomyrina dichotoma"] <-
    "Coleoptera"
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Allomyrina dichotoma"] <-
    "Scarabaeidae"
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Allomyrina dichotoma"] <-
    "Trypoxylus"
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Allomyrina dichotoma"] <-
    "Trypoxylus dichotomus"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Agaricus pattersonae"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Agaricus pattersonae"] <-
    "Fungi"
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Agaricus pattersonae"] <-
    "Basidiomycota"
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Agaricus pattersonae"] <-
    "Agaricomycetes"
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Agaricus pattersonae"] <-
    "Agaricales"
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Agaricus pattersonae"] <-
    "Agaricaceae"
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Agaricus pattersonae"] <-
    "Agaricus"
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Agaricus pattersonae"] <-
    "Agaricus pattersoniae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Melaphis chinensis"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Melaphis chinensis"] <-
    "Animalia"
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Melaphis chinensis"] <-
    "Arthropoda"
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Melaphis chinensis"] <-
    "Insecta"
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Melaphis chinensis"] <-
    "Hemiptera"
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Melaphis chinensis"] <-
    "Aphididae"
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Melaphis chinensis"] <-
    "Schlechtendalia"
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Melaphis chinensis"] <-
    "Schlechtendalia chinensis"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Chloroclysta truncata"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Chloroclysta truncata"] <-
    "Animalia"
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Chloroclysta truncata"] <-
    "Arthropoda"
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Chloroclysta truncata"] <-
    "Insecta"
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Chloroclysta truncata"] <-
    "Lepidoptera"
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Chloroclysta truncata"] <-
    "Geometridae"
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Chloroclysta truncata"] <-
    "Dysstroma"
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Chloroclysta truncata"] <-
    "Dysstroma truncata"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Lindenbergia urticaefolia"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Lindenbergia urticaefolia"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Lindenbergia urticaefolia"] <-
    "Tracheophyta"
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Lindenbergia urticaefolia"] <-
    "Magnoliopsida"
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Lindenbergia urticaefolia"] <-
    "Lamiales"
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Lindenbergia urticaefolia"] <-
    "Orobanchaceae"
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Lindenbergia urticaefolia"] <-
    "Lindenbergia"
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Lindenbergia urticaefolia"] <-
    "Lindenbergia urticifolia"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Tetraselmis chui"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Tetraselmis chui"] <-
    "Bacteria"
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Tetraselmis chui"] <-
    "Chlorophyta"
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Tetraselmis chui"] <-
    "Chlorodendrophyceae"
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Tetraselmis chui"] <-
    "Chlorodendrales"
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Tetraselmis chui"] <-
    "Chlorodendraceae"
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Tetraselmis chui"] <-
    "Tetraselmis"
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Tetraselmis chui"] <-
    "Tetraselmis chuii"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Cyanospira rippkae"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Cyanospira rippkae"] <-
    "Bacteria"
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Cyanospira rippkae"] <-
    "Cyanobacteria"
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Cyanospira rippkae"] <-
    "Cyanophyceae"
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Cyanospira rippkae"] <-
    "Nostocales"
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Cyanospira rippkae"] <-
    "Aphanizomenonaceae"
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Cyanospira rippkae"] <-
    "Cyanospira"
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Cyanospira rippkae"] <-
    "Cyanospira rippkae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Nicandra physaloides"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Nicandra physaloides"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Nicandra physaloides"] <-
    "Tracheophyta"
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Nicandra physaloides"] <-
    "Magnoliopsida"
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Nicandra physaloides"] <-
    "Solanales"
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Nicandra physaloides"] <-
    "Solanaceae"
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Nicandra physaloides"] <-
    "Nicandra"
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Nicandra physaloides"] <-
    "Nicandra physalodes"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Salvia shannoni"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Salvia shannoni"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Salvia shannoni"] <-
    "Tracheophyta"
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Salvia shannoni"] <-
    "Magnoliopsida"
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Salvia shannoni"] <-
    "Lamiales"
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Salvia shannoni"] <-
    "Lamiaceae"
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Salvia shannoni"] <-
    "Salvia"
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Salvia shannoni"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Streptomyces tsukubaensis"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Streptomyces tsukubaensis"] <-
    "Bacteria"
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Streptomyces tsukubaensis"] <-
    "Actinobacteria"
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Streptomyces tsukubaensis"] <-
    "Actinobacteria"
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Streptomyces tsukubaensis"] <-
    "Actinomycetales"
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Streptomyces tsukubaensis"] <-
    "Streptomycetaceae"
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Streptomyces tsukubaensis"] <-
    "Streptomyces"
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Streptomyces tsukubaensis"] <-
    "Streptomyces tsukubensis"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Evea brasiliensis"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Evea brasiliensis"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Evea brasiliensis"] <-
    "Tracheophyta"
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Evea brasiliensis"] <-
    "Magnoliopsida"
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Evea brasiliensis"] <-
    "Malpighiales"
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Evea brasiliensis"] <-
    "Euphorbiaceae"
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Evea brasiliensis"] <-
    "Hevea"
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Evea brasiliensis"] <-
    "Hevea brasiliensis"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismTranslated == "japanese yew"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismTranslated == "japanese yew"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organismTranslated == "japanese yew"] <-
    "Tracheophyta"
  inhouse_db$organism_3_class[inhouse_db$organismTranslated == "japanese yew"] <-
    "Pinopsida"
  inhouse_db$organism_4_order[inhouse_db$organismTranslated == "japanese yew"] <-
    "Pinales"
  inhouse_db$organism_5_family[inhouse_db$organismTranslated == "japanese yew"] <-
    "Taxaceae"
  inhouse_db$organism_6_genus[inhouse_db$organismTranslated == "japanese yew"] <-
    "Taxus"
  inhouse_db$organism_7_species[inhouse_db$organismTranslated == "japanese yew"] <-
    "Taxus cuspidata"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismTranslated == "Galla Chinensis Rhus"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismTranslated == "Galla Chinensis Rhus"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organismTranslated == "Galla Chinensis Rhus"] <-
    "Tracheophyta"
  inhouse_db$organism_3_class[inhouse_db$organismTranslated == "Galla Chinensis Rhus"] <-
    "Magnoliopsida"
  inhouse_db$organism_4_order[inhouse_db$organismTranslated == "Galla Chinensis Rhus"] <-
    "Sapindales"
  inhouse_db$organism_5_family[inhouse_db$organismTranslated == "Galla Chinensis Rhus"] <-
    "Anacardiaceae"
  inhouse_db$organism_6_genus[inhouse_db$organismTranslated == "Galla Chinensis Rhus"] <-
    "Rhus"
  inhouse_db$organism_7_species[inhouse_db$organismTranslated == "Galla Chinensis Rhus"] <-
    "Rhus chinensis"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Chinensis"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Chinensis"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Chinensis"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Chinensis"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Chinensis"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Chinensis"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Chinensis"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Chinensis"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Sinensis"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Sinensis"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Sinensis"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Sinensis"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Sinensis"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Sinensis"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Sinensis"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Sinensis"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Ootheca"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Ootheca"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Ootheca"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Ootheca"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Ootheca"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Ootheca"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Ootheca"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Ootheca"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Uncultured"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Uncultured"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Uncultured"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Uncultured"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Uncultured"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Uncultured"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Uncultured"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Uncultured"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Stigma"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Stigma"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Stigma"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Stigma"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Stigma"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Stigma"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Stigma"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Stigma"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Spica"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Spica"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Spica"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Spica"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Spica"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Spica"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Spica"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Spica"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Semen"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Semen"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Semen"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Semen"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Semen"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Semen"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Semen"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Semen"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Rotundus"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Rotundus"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Rotundus"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Rotundus"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Rotundus"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Rotundus"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Rotundus"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Rotundus"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Rhizoma"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Rhizoma"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Rhizoma"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Rhizoma"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Rhizoma"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Rhizoma"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Rhizoma"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Rhizoma"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Ramulus"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Ramulus"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Ramulus"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Ramulus"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Ramulus"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Ramulus"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Ramulus"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Ramulus"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Radix"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Radix"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Radix"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Radix"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Radix"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Radix"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Radix"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Radix"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Pollen"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Pollen"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Pollen"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Pollen"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Pollen"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Pollen"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Pollen"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Pollen"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Lignum"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Lignum"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Lignum"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Lignum"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Lignum"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Lignum"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Lignum"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Lignum"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Fructus"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Fructus"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Fructus"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Fructus"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Fructus"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Fructus"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Fructus"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Fructus"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Flos"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Flos"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Flos"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Flos"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Flos"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Flos"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Flos"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Flos"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Corolla"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Corolla"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Corolla"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Corolla"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Corolla"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Corolla"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Corolla"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Corolla"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Cacumen"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Cacumen"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Cacumen"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Cacumen"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Cacumen"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Cacumen"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Cacumen"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Cacumen"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Bulbus"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Bulbus"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Bulbus"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Bulbus"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Bulbus"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Bulbus"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Bulbus"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Bulbus"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Megaleia"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Megaleia"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Megaleia"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Megaleia"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Megaleia"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Megaleia"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Megaleia"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Megaleia"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Candidatus"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Candidatus"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Candidatus"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Candidatus"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Candidatus"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Candidatus"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Candidatus"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Candidatus"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Tasmanian"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Tasmanian"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Tasmanian"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Tasmanian"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Tasmanian"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Tasmanian"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Tasmanian"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Tasmanian"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Asian"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Asian"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Asian"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Asian"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Asian"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Asian"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Asian"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Asian"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Mammalian"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Mammalian"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Mammalian"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Mammalian"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Mammalian"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Mammalian"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Mammalian"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Mammalian"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Red"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Red"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Red"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Red"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Red"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Red"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Red"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Red"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Turkey"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Turkey"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Turkey"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Turkey"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Turkey"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Turkey"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Turkey"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Turkey"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Synthetis"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Synthetis"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Synthetis"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Synthetis"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Synthetis"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Synthetis"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Synthetis"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Synthetis"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Pagellus erythrinus"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Pagellus erythrinus"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Pagellus erythrinus"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Pagellus erythrinus"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Pagellus erythrinus"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Pagellus erythrinus"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Pagellus erythrinus"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Pagellus erythrinus"] <-
    ""
  # comes from Becker translation

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Lagenorhynchus obliquidens"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Lagenorhynchus obliquidens"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Lagenorhynchus obliquidens"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Lagenorhynchus obliquidens"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Lagenorhynchus obliquidens"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Lagenorhynchus obliquidens"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Lagenorhynchus obliquidens"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Lagenorhynchus obliquidens"] <-
    ""
  # comes from Lag translation

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Actinomycetales bacterium	"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Actinomycetales bacterium	"] <-
    "Bacteria"
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Actinomycetales bacterium	"] <-
    "Actinobacteria"
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Actinomycetales bacterium	"] <-
    "Actinobacteria"
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Actinomycetales bacterium	"] <-
    "Actinomycetales"
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Actinomycetales bacterium	"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Actinomycetales bacterium	"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Actinomycetales bacterium	"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Galla"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Galla"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Galla"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Galla"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Galla"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Galla"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Galla"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Galla"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Red Sea bacterium KT-2K1"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Red Sea bacterium KT-2K1"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Red Sea bacterium KT-2K1"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Red Sea bacterium KT-2K1"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Red Sea bacterium KT-2K1"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Red Sea bacterium KT-2K1"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Red Sea bacterium KT-2K1"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Red Sea bacterium KT-2K1"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Green Pelican GFP transformation vector"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Green Pelican GFP transformation vector"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Green Pelican GFP transformation vector"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Green Pelican GFP transformation vector"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Green Pelican GFP transformation vector"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Green Pelican GFP transformation vector"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Green Pelican GFP transformation vector"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Green Pelican GFP transformation vector"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Peripatoides"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Peripatoides"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Peripatoides"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Peripatoides"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Peripatoides"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Peripatoides"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Peripatoides"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Peripatoides"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Bacterium MPBA1"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Bacterium MPBA1"] <-
    "Bacteria"
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Bacterium MPBA1"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Bacterium MPBA1"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Bacterium MPBA1"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Bacterium MPBA1"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Bacterium MPBA1"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Bacterium MPBA1"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Rhizobiaceae bacterium"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Rhizobiaceae bacterium"] <-
    "Bacteria"
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Rhizobiaceae bacterium"] <-
    "Proteobacteria"
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Rhizobiaceae bacterium"] <-
    "Alphaproteobacteria"
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Rhizobiaceae bacterium"] <-
    "Rhizobiales"
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Rhizobiaceae bacterium"] <-
    "Rhizobiaceae"
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Rhizobiaceae bacterium"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Rhizobiaceae bacterium"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Cyanophyta"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Cyanophyta"] <-
    "Bacteria"
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Cyanophyta"] <-
    "Cyanobacteria"
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Cyanophyta"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Cyanophyta"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Cyanophyta"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Cyanophyta"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Cyanophyta"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Candida"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Candida"] <-
    "Fungi"
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Candida"] <-
    "Ascomycota"
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Candida"] <-
    "Saccharomycetes"
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Candida"] <-
    "Saccharomycetales"
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Candida"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Candida"] <-
    "Candida"
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Candida"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Microsorium"] <-
    "y"
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Microsorium"] <-
    "Microsorum"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Pseudoeurotium"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Pseudoeurotium"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Pseudoeurotium"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Pseudoeurotium"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Pseudoeurotium"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Pseudoeurotium"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Pseudoeurotium"] <-
    "Pseudeurotium"
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Pseudoeurotium"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Lanea"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Lanea"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Lanea"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Lanea"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Lanea"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Lanea"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Lanea"] <-
    "Lannea"
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Lanea"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Aspidium"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Aspidium"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Aspidium"] <-
    "Tracheophyta"
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Aspidium"] <-
    "Polypodiopsida"
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Aspidium"] <-
    "Polypodiales"
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Aspidium"] <-
    "Dryopteridaceae"
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Aspidium"] <-
    "Polystichum"
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Aspidium"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Iris"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Iris"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Iris"] <-
    "Tracheophyta"
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Iris"] <-
    "Liliopsida"
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Iris"] <-
    "Asparagales"
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Iris"] <-
    "Iridaceae"
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Iris"] <-
    "Iris"
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Iris"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Plantae"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Plantae"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Plantae"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Plantae"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Plantae"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Plantae"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Plantae"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Plantae"] <-
    ""

  # sadly can be multiple kingdoms, therefore droped
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Algae"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Algae"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Algae"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Algae"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Algae"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Algae"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Algae"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Algae"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Fungi"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Fungi"] <-
    "Fungi"
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Fungi"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Fungi"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Fungi"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Fungi"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Fungi"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Fungi"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Anaerobic"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Anaerobic"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Anaerobic"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Anaerobic"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Anaerobic"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Anaerobic"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Anaerobic"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismCleaned == "Anaerobic"] <-
    ""

  # inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Myxococcus hansupus"] <-
  #   "y"
  # inhouse_db$organism_6_genus[inhouse_db$organism_7_species == "Myxococcus hansupus"] <-
  #   "Myxococcus"

  # inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Aubrietia deltoidea"] <-
  #   "y"
  # inhouse_db$organism_6_genus[inhouse_db$organism_7_species == "Aubrietia deltoidea"] <-
  #   "Aubrietia"

  # inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Cystoseira granulata"] <-
  #   "y"
  # inhouse_db$organism_6_genus[inhouse_db$organism_7_species == "Cystoseira granulata"] <-
  #   "Cystoseira"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Chromulina ochromonoides"] <-
    "y"
  inhouse_db$organism_5_family[inhouse_db$organism_7_species == "Chromulina ochromonoides"] <-
    "Chromulinaceae"
  inhouse_db$organism_6_genus[inhouse_db$organism_7_species == "Chromulina ochromonoides"] <-
    "Chromulina"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Microchaete loktakensis"] <-
    "y"
  inhouse_db$organism_4_order[inhouse_db$organism_7_species == "Microchaete loktakensis"] <-
    "Nostocales"
  inhouse_db$organism_5_family[inhouse_db$organism_7_species == "Microchaete loktakensis"] <-
    "Microchaetaceae"
  inhouse_db$organism_6_genus[inhouse_db$organism_7_species == "Microchaete loktakensis"] <-
    "Microchaete"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Umezakia natans"] <-
    "y"
  inhouse_db$organism_4_order[inhouse_db$organism_7_species == "Umezakia natans"] <-
    "Nostocales"
  inhouse_db$organism_5_family[inhouse_db$organism_7_species == "Umezakia natans"] <-
    "Aphanizomenonaceae"
  inhouse_db$organism_6_genus[inhouse_db$organism_7_species == "Umezakia natans"] <-
    "Umezakia"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Aphanothece sacrum"] <-
    "y"
  inhouse_db$organism_4_order[inhouse_db$organism_7_species == "Aphanothece sacrum"] <-
    "Chroococcales"
  inhouse_db$organism_5_family[inhouse_db$organism_7_species == "Aphanothece sacrum"] <-
    "Aphanothecaceae"
  inhouse_db$organism_6_genus[inhouse_db$organism_7_species == "Aphanothece sacrum"] <-
    "Aphanothece"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Aphanothece halophytica"] <-
    "y"
  inhouse_db$organism_4_order[inhouse_db$organism_7_species == "Aphanothece halophytica"] <-
    "Chroococcales"
  inhouse_db$organism_5_family[inhouse_db$organism_7_species == "Aphanothece halophytica"] <-
    "Aphanothecaceae"
  inhouse_db$organism_6_genus[inhouse_db$organism_7_species == "Aphanothece halophytica"] <-
    "Aphanothece"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Syringoderma phinneyi"] <-
    "y"
  inhouse_db$organism_4_order[inhouse_db$organism_7_species == "Syringoderma phinneyi"] <-
    "Syringodermatales"
  inhouse_db$organism_5_family[inhouse_db$organism_7_species == "Syringoderma phinneyi"] <-
    "Syringodermataceae"
  inhouse_db$organism_6_genus[inhouse_db$organism_7_species == "Syringoderma phinneyi"] <-
    "Syringoderma"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Scytonema varium"] <-
    "y"
  inhouse_db$organism_2_phylum[inhouse_db$organism_7_species == "Scytonema varium"] <-
    "Cyanobacteria"
  inhouse_db$organism_4_order[inhouse_db$organism_7_species == "Scytonema varium"] <-
    "Nostocales"
  inhouse_db$organism_5_family[inhouse_db$organism_7_species == "Scytonema varium"] <-
    "Scytonemataceae"
  inhouse_db$organism_6_genus[inhouse_db$organism_7_species == "Scytonema varium"] <-
    "Scytonema"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Scytonema hofmanni"] <-
    "y"
  inhouse_db$organism_2_phylum[inhouse_db$organism_7_species == "Scytonema hofmanni"] <-
    "Cyanobacteria"
  inhouse_db$organism_4_order[inhouse_db$organism_7_species == "Scytonema hofmanni"] <-
    "Nostocales"
  inhouse_db$organism_5_family[inhouse_db$organism_7_species == "Scytonema hofmanni"] <-
    "Scytonemataceae"
  inhouse_db$organism_6_genus[inhouse_db$organism_7_species == "Scytonema hofmanni"] <-
    "Scytonema"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Cylindrospermum muscicola"] <-
    "y"
  inhouse_db$organism_2_phylum[inhouse_db$organism_7_species == "Cylindrospermum muscicola"] <-
    "Cyanobacteria"
  inhouse_db$organism_4_order[inhouse_db$organism_7_species == "Cylindrospermum muscicola"] <-
    "Nostocales"
  inhouse_db$organism_5_family[inhouse_db$organism_7_species == "Cylindrospermum muscicola"] <-
    "Nostocaceae"
  inhouse_db$organism_6_genus[inhouse_db$organism_7_species == "Cylindrospermum muscicola"] <-
    "Cylindrospermum"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Tolypothrix nodosa"] <-
    "y"
  inhouse_db$organism_2_phylum[inhouse_db$organism_7_species == "Tolypothrix nodosa"] <-
    "Cyanobacteria"
  inhouse_db$organism_4_order[inhouse_db$organism_7_species == "Tolypothrix nodosa"] <-
    "Nostocales"
  inhouse_db$organism_5_family[inhouse_db$organism_7_species == "Tolypothrix nodosa"] <-
    "Scytonemataceae"
  inhouse_db$organism_6_genus[inhouse_db$organism_7_species == "Tolypothrix nodosa"] <-
    "Tolypothrix"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Tolypothrix tjipanasensis"] <-
    "y"
  inhouse_db$organism_2_phylum[inhouse_db$organism_7_species == "Tolypothrix tjipanasensis"] <-
    "Cyanobacteria"
  inhouse_db$organism_4_order[inhouse_db$organism_7_species == "Tolypothrix tjipanasensis"] <-
    "Nostocales"
  inhouse_db$organism_5_family[inhouse_db$organism_7_species == "Tolypothrix tjipanasensis"] <-
    "Scytonemataceae"
  inhouse_db$organism_6_genus[inhouse_db$organism_7_species == "Tolypothrix tjipanasensis"] <-
    "Tolypothrix"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Dolichospermum flos-aquae"] <-
    "y"
  inhouse_db$organism_6_genus[inhouse_db$organism_7_species == "Dolichospermum flos-aquae"] <-
    "Dolichospermum"

  # inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Cystoseira barbata"] <-
  #   "y"
  # inhouse_db$organism_6_genus[inhouse_db$organism_7_species == "Cystoseira barbata"] <-
  #   "Cystoseira"

  # inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Cystoseira barbatula"] <-
  #   "y"
  # inhouse_db$organism_6_genus[inhouse_db$organism_7_species == "Cystoseira barbatula"] <-
  #   "Cystoseira"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Scytosiphon lomentaria"] <-
    "y"
  inhouse_db$organism_3_class[inhouse_db$organism_7_species == "Scytosiphon lomentaria"] <-
    "Phaeophyceae"
  inhouse_db$organism_4_order[inhouse_db$organism_7_species == "Scytosiphon lomentaria"] <-
    "Ectocarpales"
  inhouse_db$organism_5_family[inhouse_db$organism_7_species == "Scytosiphon lomentaria"] <-
    "Scytosiphonaceae"
  inhouse_db$organism_6_genus[inhouse_db$organism_7_species == "Scytosiphon lomentaria"] <-
    "Scytosiphon"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Homoeostrichus sinclairii"] <-
    "y"
  inhouse_db$organism_3_class[inhouse_db$organism_7_species == "Homoeostrichus sinclairii"] <-
    "Phaeophyceae"
  inhouse_db$organism_4_order[inhouse_db$organism_7_species == "Homoeostrichus sinclairii"] <-
    "Dictyotales"
  inhouse_db$organism_5_family[inhouse_db$organism_7_species == "Homoeostrichus sinclairii"] <-
    "Dictyotaceae"
  inhouse_db$organism_6_genus[inhouse_db$organism_7_species == "Homoeostrichus sinclairii"] <-
    "Homoeostrichus"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Dictyopteris prolifera"] <-
    "y"
  inhouse_db$organism_3_class[inhouse_db$organism_7_species == "Dictyopteris prolifera"] <-
    "Phaeophyceae"
  inhouse_db$organism_4_order[inhouse_db$organism_7_species == "Dictyopteris prolifera"] <-
    "Dictyotales"
  inhouse_db$organism_5_family[inhouse_db$organism_7_species == "Dictyopteris prolifera"] <-
    "Dictyotaceae"
  inhouse_db$organism_6_genus[inhouse_db$organism_7_species == "Dictyopteris prolifera"] <-
    "Dictyopteris"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Xiphophora gladiata"] <-
    "y"
  inhouse_db$organism_6_genus[inhouse_db$organism_7_species == "Xiphophora gladiata"] <-
    "Xiphophora"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Iyengaria stellata"] <-
    "y"
  inhouse_db$organism_3_class[inhouse_db$organism_7_species == "Iyengaria stellata"] <-
    "Phaeophyceae"
  inhouse_db$organism_4_order[inhouse_db$organism_7_species == "Iyengaria stellata"] <-
    "Ectocarpales"
  inhouse_db$organism_5_family[inhouse_db$organism_7_species == "Iyengaria stellata"] <-
    "Scytosiphonaceae"
  inhouse_db$organism_6_genus[inhouse_db$organism_7_species == "Iyengaria stellata"] <-
    "Iyengaria"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Astrosporangium hypotensionis"] <-
    "y"
  inhouse_db$organism_2_phylum[inhouse_db$organism_7_species == "Astrosporangium hypotensionis"] <-
    "Actinobacteria"
  inhouse_db$organism_3_class[inhouse_db$organism_7_species == "Astrosporangium hypotensionis"] <-
    "Actinobacteria"
  inhouse_db$organism_4_order[inhouse_db$organism_7_species == "Astrosporangium hypotensionis"] <-
    "Actinomycetales"
  inhouse_db$organism_5_family[inhouse_db$organism_7_species == "Astrosporangium hypotensionis"] <-
    "Streptosporangiaceae"
  inhouse_db$organism_6_genus[inhouse_db$organism_7_species == "Astrosporangium hypotensionis"] <-
    "Astrosporangium"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Caulerpa scalpelliformis"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_7_species == "Caulerpa scalpelliformis"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organism_7_species == "Caulerpa scalpelliformis"] <-
    "Chlorophyta"
  inhouse_db$organism_3_class[inhouse_db$organism_7_species == "Caulerpa scalpelliformis"] <-
    "Ulvophyceae"
  inhouse_db$organism_4_order[inhouse_db$organism_7_species == "Caulerpa scalpelliformis"] <-
    "Bryopsidales"
  inhouse_db$organism_5_family[inhouse_db$organism_7_species == "Caulerpa scalpelliformis"] <-
    "Caulerpaceae"
  inhouse_db$organism_6_genus[inhouse_db$organism_7_species == "Caulerpa scalpelliformis"] <-
    "Caulerpa"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Gram-negative bacterium"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_7_species == "Gram-negative bacterium"] <-
    "Bacteria"
  inhouse_db$organism_2_phylum[inhouse_db$organism_7_species == "Gram-negative bacterium"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organism_7_species == "Gram-negative bacterium"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organism_7_species == "Gram-negative bacterium"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organism_7_species == "Gram-negative bacterium"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organism_7_species == "Gram-negative bacterium"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organism_7_species == "Gram-negative bacterium"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Laurencia viridis"] <-
    "y"
  inhouse_db$organism_4_order[inhouse_db$organism_7_species == "Laurencia viridis"] <-
    "Ceramiales"
  inhouse_db$organism_5_family[inhouse_db$organism_7_species == "Laurencia viridis"] <-
    "Rhodomelaceae"
  inhouse_db$organism_6_genus[inhouse_db$organism_7_species == "Laurencia viridis"] <-
    "Laurencia"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Laurencia complanata"] <-
    "y"
  inhouse_db$organism_4_order[inhouse_db$organism_7_species == "Laurencia complanata"] <-
    "Ceramiales"
  inhouse_db$organism_5_family[inhouse_db$organism_7_species == "Laurencia complanata"] <-
    "Rhodomelaceae"
  inhouse_db$organism_6_genus[inhouse_db$organism_7_species == "Laurencia complanata"] <-
    "Laurencia"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Chondracanthus harveyanus"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_7_species == "Chondracanthus harveyanus"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organism_7_species == "Chondracanthus harveyanus"] <-
    "Rhodophyta"
  inhouse_db$organism_3_class[inhouse_db$organism_7_species == "Chondracanthus harveyanus"] <-
    "Florideophyceae"
  inhouse_db$organism_4_order[inhouse_db$organism_7_species == "Chondracanthus harveyanus"] <-
    "Gigartinales"
  inhouse_db$organism_5_family[inhouse_db$organism_7_species == "Chondracanthus harveyanus"] <-
    "Gigartinaceae"
  inhouse_db$organism_6_genus[inhouse_db$organism_7_species == "Chondracanthus harveyanus"] <-
    "Chondracanthus"

  # examples mismatched genera
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Streptomyces varius"] <-
    "y"
  inhouse_db$organism_6_genus[inhouse_db$organism_7_species == "Streptomyces varius"] <-
    "Streptomyces"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Streptomyces gangtokensis"] <-
    "y"
  inhouse_db$organism_6_genus[inhouse_db$organism_7_species == "Streptomyces gangtokensis"] <-
    "Streptomyces"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Streptomyces kaniharaensis"] <-
    "y"
  inhouse_db$organism_6_genus[inhouse_db$organism_7_species == "Streptomyces kaniharaensis"] <-
    "Streptomyces"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Mus striatus"] <-
    "y"
  inhouse_db$organism_6_genus[inhouse_db$organism_7_species == "Mus striatus"] <-
    "Mus"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Peridinium foliaceum"] <-
    "y"
  inhouse_db$organism_6_genus[inhouse_db$organism_7_species == "Peridinium foliaceum"] <-
    "Peridinium"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Agrobacterium aurantiacum"] <-
    "y"
  inhouse_db$organism_6_genus[inhouse_db$organism_7_species == "Agrobacterium aurantiacum"] <-
    "Agrobacterium"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Aucklandia lappa"] <-
    "y"
  inhouse_db$organism_6_genus[inhouse_db$organism_7_species == "Aucklandia lappa"] <-
    "Aucklandia"

  # example
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Solanum etuberosum"] <-
    "y"
  inhouse_db$organism_7_species[inhouse_db$organism_7_species == "Solanum etuberosum"] <-
    "Solanum tuberosum"

  # example
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Toricellia angulata"] <-
    "y"
  inhouse_db$organism_7_species[inhouse_db$organism_7_species == "Toricellia angulata"] <-
    "Torricellia angulata"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Boerhaavia diffusa"] <-
    "y"
  inhouse_db$organism_7_species[inhouse_db$organism_7_species == "Boerhaavia diffusa"] <-
    "Boerhavia diffusa"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Boerhaavia coccinea"] <-
    "y"
  inhouse_db$organism_7_species[inhouse_db$organism_7_species == "Boerhaavia coccinea"] <-
    "Boerhavia coccinea"

  # mismatched genus
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Achromobacter cycloclastes"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_7_species == "Achromobacter cycloclastes"] <-
    "Bacteria"
  inhouse_db$organism_2_phylum[inhouse_db$organism_7_species == "Achromobacter cycloclastes"] <-
    "Proteobacteria"
  inhouse_db$organism_3_class[inhouse_db$organism_7_species == "Achromobacter cycloclastes"] <-
    "Betaproteobacteria"
  inhouse_db$organism_4_order[inhouse_db$organism_7_species == "Achromobacter cycloclastes"] <-
    "Burkholderiales"
  inhouse_db$organism_5_family[inhouse_db$organism_7_species == "Achromobacter cycloclastes"] <-
    "Alcaligenaceae"
  inhouse_db$organism_6_genus[inhouse_db$organism_7_species == "Achromobacter cycloclastes"] <-
    "Achromobacter"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Pseudomonas reactans"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_7_species == "Pseudomonas reactans"] <-
    "Bacteria"
  inhouse_db$organism_2_phylum[inhouse_db$organism_7_species == "Pseudomonas reactans"] <-
    "Proteobacteria"
  inhouse_db$organism_3_class[inhouse_db$organism_7_species == "Pseudomonas reactans"] <-
    "Gammaproteobacteria"
  inhouse_db$organism_4_order[inhouse_db$organism_7_species == "Pseudomonas reactans"] <-
    "Pseudomonadales"
  inhouse_db$organism_5_family[inhouse_db$organism_7_species == "Pseudomonas reactans"] <-
    "Pseudomonadaceae"
  inhouse_db$organism_6_genus[inhouse_db$organism_7_species == "Pseudomonas reactans"] <-
    "Pseudomonas"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Micromonospora megalomicea"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_7_species == "Micromonospora megalomicea"] <-
    "Bacteria"
  inhouse_db$organism_2_phylum[inhouse_db$organism_7_species == "Micromonospora megalomicea"] <-
    "Actinobacteria"
  inhouse_db$organism_3_class[inhouse_db$organism_7_species == "Micromonospora megalomicea"] <-
    "Actinobacteria"
  inhouse_db$organism_4_order[inhouse_db$organism_7_species == "Micromonospora megalomicea"] <-
    "Actinomycetales"
  inhouse_db$organism_5_family[inhouse_db$organism_7_species == "Micromonospora megalomicea"] <-
    "Micromonosporaceae"
  inhouse_db$organism_6_genus[inhouse_db$organism_7_species == "Micromonospora megalomicea"] <-
    "Micromonospora"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Hypnea musciformis"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_7_species == "Hypnea musciformis"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organism_7_species == "Hypnea musciformis"] <-
    "Rhodophyta"
  inhouse_db$organism_3_class[inhouse_db$organism_7_species == "Hypnea musciformis"] <-
    "Florideophyceae"
  inhouse_db$organism_4_order[inhouse_db$organism_7_species == "Hypnea musciformis"] <-
    "Gigartinales"
  inhouse_db$organism_5_family[inhouse_db$organism_7_species == "Hypnea musciformis"] <-
    "Cystocloniaceae"
  inhouse_db$organism_6_genus[inhouse_db$organism_7_species == "Hypnea musciformis"] <-
    "Hypnea"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Hypnea valentiae"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_7_species == "Hypnea valentiae"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organism_7_species == "Hypnea valentiae"] <-
    "Rhodophyta"
  inhouse_db$organism_3_class[inhouse_db$organism_7_species == "Hypnea valentiae"] <-
    "Florideophyceae"
  inhouse_db$organism_4_order[inhouse_db$organism_7_species == "Hypnea valentiae"] <-
    "Gigartinales"
  inhouse_db$organism_5_family[inhouse_db$organism_7_species == "Hypnea valentiae"] <-
    "Cystocloniaceae"
  inhouse_db$organism_6_genus[inhouse_db$organism_7_species == "Hypnea valentiae"] <-
    "Hypnea"

  # example
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Spirodela polyrrhiza"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_7_species == "Spirodela polyrrhiza"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organism_7_species == "Spirodela polyrrhiza"] <-
    "Tracheophyta"
  inhouse_db$organism_3_class[inhouse_db$organism_7_species == "Spirodela polyrrhiza"] <-
    "Magnoliopsida"
  inhouse_db$organism_4_order[inhouse_db$organism_7_species == "Spirodela polyrrhiza"] <-
    "Alismatales"
  inhouse_db$organism_5_family[inhouse_db$organism_7_species == "Spirodela polyrrhiza"] <-
    "Araceae"
  inhouse_db$organism_6_genus[inhouse_db$organism_7_species == "Spirodela polyrrhiza"] <-
    "Spirodela"
  inhouse_db$organism_7_species[inhouse_db$organism_7_species == "Spirodela polyrrhiza"] <-
    "Spirodela polyrhiza"

  # example_2
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Hyacinthoides nonscripta"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_7_species == "Hyacinthoides nonscripta"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organism_7_species == "Hyacinthoides nonscripta"] <-
    "Tracheophyta"
  inhouse_db$organism_3_class[inhouse_db$organism_7_species == "Hyacinthoides nonscripta"] <-
    "Magnoliopsida"
  inhouse_db$organism_4_order[inhouse_db$organism_7_species == "Hyacinthoides nonscripta"] <-
    "Asparagales"
  inhouse_db$organism_5_family[inhouse_db$organism_7_species == "Hyacinthoides nonscripta"] <-
    "Asparagaceae"
  inhouse_db$organism_6_genus[inhouse_db$organism_7_species == "Hyacinthoides nonscripta"] <-
    "Hyacinthoides"
  inhouse_db$organism_7_species[inhouse_db$organism_7_species == "Hyacinthoides nonscripta"] <-
    "Hyacinthoides non-scripta"

  # double taxonomies -> no "y", just choose one
  # catalogue of Life as reference, then NCBI

  inhouse_db$organism_1_kingdom[inhouse_db$organism_7_species == "Turbinaria ornata"] <-
    "Chromista"
  inhouse_db$organism_2_phylum[inhouse_db$organism_7_species == "Turbinaria ornata"] <-
    "Ochrophyta"
  inhouse_db$organism_3_class[inhouse_db$organism_7_species == "Turbinaria ornata"] <-
    "Phaeophyceae"
  inhouse_db$organism_4_order[inhouse_db$organism_7_species == "Turbinaria ornata"] <-
    "Fucales"
  inhouse_db$organism_5_family[inhouse_db$organism_7_species == "Turbinaria ornata"] <-
    "Sargassaceae"
  inhouse_db$organism_6_genus[inhouse_db$organism_7_species == "Turbinaria ornata"] <-
    "Turbinaria"

  inhouse_db$organism_2_phylum[inhouse_db$organism_6_genus == "Dolichospermum"] <-
    "Cyanobacteria"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Dolichospermum"] <-
    "Nostocales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Dolichospermum"] <-
    "Aphanizomenonaceae"

  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Aphanothece"] <-
    "Chroococcaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Eschscholtzia"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_6_genus == "Eschscholtzia"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organism_6_genus == "Eschscholtzia"] <-
    "Tracheophyta"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Eschscholtzia"] <-
    "Magnoliopsida"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Eschscholtzia"] <-
    "Ranunculales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Eschscholtzia"] <-
    "Papaveraceae"
  inhouse_db$organism_6_genus[inhouse_db$organism_6_genus == "Eschscholtzia"] <-
    "Eschscholzia"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Ishige"] <-
    "y"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Ishige"] <-
    "Ishigeaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Scytothamnus"] <-
    "y"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Scytothamnus"] <-
    "Phaeophyceae"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Scytothamnus"] <-
    "Scytothamnales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Scytothamnus"] <-
    "Splachnidiaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Cystoseira"] <-
    "y"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Cystoseira"] <-
    "Phaeophyceae"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Cystoseira"] <-
    "Fucales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Cystoseira"] <-
    "Sargassaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Nigritella"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_6_genus == "Nigritella"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organism_6_genus == "Nigritella"] <-
    "Tracheophyta"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Nigritella"] <-
    "Magnoliopsida"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Nigritella"] <-
    "Asparagales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Nigritella"] <-
    "Orchidaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Lannea"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_6_genus == "Lannea"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organism_6_genus == "Lannea"] <-
    "Tracheophyta"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Lannea"] <-
    "Magnoliopsida"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Lannea"] <-
    "Sapindales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Lannea"] <-
    "Anacardiaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Dichotomomyces"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_6_genus == "Dichotomomyces"] <-
    "Fungi"
  inhouse_db$organism_2_phylum[inhouse_db$organism_6_genus == "Dichotomomyces"] <-
    "Ascomycota"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Dichotomomyces"] <-
    "Eurotiomycetes"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Dichotomomyces"] <-
    "Eurotiales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Dichotomomyces"] <-
    "Aspergillaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Dichotomomyces"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_6_genus == "Dichotomomyces"] <-
    "Fungi"
  inhouse_db$organism_2_phylum[inhouse_db$organism_6_genus == "Dichotomomyces"] <-
    "Ascomycota"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Dichotomomyces"] <-
    "Eurotiomycetes"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Dichotomomyces"] <-
    "Eurotiales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Dichotomomyces"] <-
    "Aspergillaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Microsorum"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_6_genus == "Microsorum"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organism_6_genus == "Microsorum"] <-
    "Tracheophyta"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Microsorum"] <-
    "Polypodiopsida"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Microsorum"] <-
    "Polypodiales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Microsorum"] <-
    "Polypodiaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Leishmania"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_6_genus == "Leishmania"] <-
    "Protozoa"
  inhouse_db$organism_2_phylum[inhouse_db$organism_6_genus == "Leishmania"] <-
    "Euglenozoa"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Leishmania"] <-
    "Kinetoplastea"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Leishmania"] <-
    "Trypanosomatida"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Leishmania"] <-
    "Trypanosomatidae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Pseudogymnoascus"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_6_genus == "Pseudogymnoascus"] <-
    "Fungi"
  inhouse_db$organism_2_phylum[inhouse_db$organism_6_genus == "Pseudogymnoascus"] <-
    "Ascomycota"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Pseudogymnoascus"] <-
    "Dothideomycetes"
  # no organism_4_order
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Pseudogymnoascus"] <-
    "Pseudeurotiaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Methanothrix"] <-
    "y"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Methanothrix"] <-
    "Methanotrichaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Daphnia"] <-
    "y"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Daphnia"] <-
    "Daphniidae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Ophryoscolex"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_6_genus == "Ophryoscolex"] <-
    "Protozoa"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Ophryoscolex"] <-
    "Entodiniomorphida"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Ophryoscolex"] <-
    "Ophryoscolecidae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Plasmodium"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_6_genus == "Plasmodium"] <-
    "Protozoa"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Plasmodium"] <-
    "Haemospororida"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Plasmodium"] <-
    "Plasmodiidae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Megalobulimus"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_6_genus == "Megalobulimus"] <-
    "Animalia"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Megalobulimus"] <-
    "Gastropoda"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Megalobulimus"] <-
    "Stylommatophora"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Megalobulimus"] <-
    "Megalobulimidae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Euglena"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_6_genus == "Euglena"] <-
    "Protozoa"
  inhouse_db$organism_2_phylum[inhouse_db$organism_6_genus == "Euglena"] <-
    "Euglenozoa"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Euglena"] <-
    "Euglenida"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Euglena"] <-
    "Euglenales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Euglena"] <-
    "Euglenaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Zyzza"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_6_genus == "Zyzza"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organism_6_genus == "Zyzza"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Zyzza"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Zyzza"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Zyzza"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Prochloron"] <-
    "y"
  inhouse_db$organism_2_phylum[inhouse_db$organism_6_genus == "Prochloron"] <-
    "Cyanobacteria"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Prochloron"] <-
    "Chroococcales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Prochloron"] <-
    "Prochloraceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Glomospora"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_6_genus == "Glomospora"] <-
    "Fungi"
  inhouse_db$organism_2_phylum[inhouse_db$organism_6_genus == "Glomospora"] <-
    "Basidiomycota"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Glomospora"] <-
    "Pucciniomycetes"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Glomospora"] <-
    "Platygloeales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Glomospora"] <-
    "Platygloeaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Flacourtia"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_6_genus == "Flacourtia"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organism_6_genus == "Flacourtia"] <-
    "Tracheophyta"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Flacourtia"] <-
    "Magnoliopsida"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Flacourtia"] <-
    "Malpighiales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Flacourtia"] <-
    "Salicaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Geomyces"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_6_genus == "Geomyces"] <-
    "Fungi"
  inhouse_db$organism_2_phylum[inhouse_db$organism_6_genus == "Geomyces"] <-
    "Ascomycota"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Geomyces"] <-
    "Leotiomycetes"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Geomyces"] <-
    "Helotiales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Geomyces"] <-
    "Myxotrichaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Robbsia"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_6_genus == "Robbsia"] <-
    "Bacteria"
  inhouse_db$organism_2_phylum[inhouse_db$organism_6_genus == "Robbsia"] <-
    "Proteobacteria"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Robbsia"] <-
    "Betaproteobacteria"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Robbsia"] <-
    "Burkholderiales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Robbsia"] <-
    "Burkholderiaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Metanarthecium"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_6_genus == "Metanarthecium"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organism_6_genus == "Metanarthecium"] <-
    "Tracheophyta"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Metanarthecium"] <-
    "Magnoliopsida"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Metanarthecium"] <-
    "Dioscoreales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Metanarthecium"] <-
    "Nartheciaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Tichocarpus"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_6_genus == "Tichocarpus"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organism_6_genus == "Tichocarpus"] <-
    "Rhodophyta"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Tichocarpus"] <-
    "Florideophyceae"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Tichocarpus"] <-
    "Gigartinales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Tichocarpus"] <-
    "Tichocarpaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Flabellina"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_6_genus == "Flabellina"] <-
    "Animalia"
  inhouse_db$organism_2_phylum[inhouse_db$organism_6_genus == "Flabellina"] <-
    "Mollusca"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Flabellina"] <-
    "Gastropoda"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Flabellina"] <-
    "Nudibranchia"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Flabellina"] <-
    "Flabellinidae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Calcarisporiella"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_6_genus == "Calcarisporiella"] <-
    "Fungi"
  inhouse_db$organism_2_phylum[inhouse_db$organism_6_genus == "Calcarisporiella"] <-
    "Mucoromycota"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Calcarisporiella"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Calcarisporiella"] <-
    "Calcarisporiellales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Calcarisporiella"] <-
    "Calcarisporiellaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Trypanosoma"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_6_genus == "Trypanosoma"] <-
    "Protozoa"
  inhouse_db$organism_2_phylum[inhouse_db$organism_6_genus == "Trypanosoma"] <-
    "Euglenozoa"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Trypanosoma"] <-
    "Kinetoplastea"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Trypanosoma"] <-
    "Trypanosomatida"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Trypanosoma"] <-
    "Trypanosomatidae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Anatheca"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_6_genus == "Anatheca"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organism_6_genus == "Anatheca"] <-
    "Rhodophyta"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Anatheca"] <-
    "Florideophyceae"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Anatheca"] <-
    "Gigartinales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Anatheca"] <-
    "Solieriaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Echinocystis"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_6_genus == "Echinocystis"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organism_6_genus == "Echinocystis"] <-
    "Tracheophyta"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Echinocystis"] <-
    "Magnoliopsida"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Echinocystis"] <-
    "Cucurbitales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Echinocystis"] <-
    "Cucurbitaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Rhabdonia"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_6_genus == "Rhabdonia"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organism_6_genus == "Rhabdonia"] <-
    "Rhodophyta"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Rhabdonia"] <-
    "Florideophyceae"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Rhabdonia"] <-
    "Gigartinales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Rhabdonia"] <-
    "Areschougiaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Hypsizigus"] <-
    "y"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Hypsizigus"] <-
    "Agaricomycetes"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Hypsizigus"] <-
    "Agaricales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Hypsizigus"] <-
    "Lyophyllaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Solenopora"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_6_genus == "Solenopora"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organism_6_genus == "Solenopora"] <-
    "Rhodophyta"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Solenopora"] <-
    "Florideophyceae"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Solenopora"] <-
    "Halymeniales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Solenopora"] <-
    "Solenoporaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismCleaned == "Mayodendron"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismCleaned == "Mayodendron"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismCleaned == "Mayodendron"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismCleaned == "Mayodendron"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismCleaned == "Mayodendron"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismCleaned == "Mayodendron"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismCleaned == "Mayodendron"] <-
    ""

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Macfadyena"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_6_genus == "Macfadyena"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organism_6_genus == "Macfadyena"] <-
    "Tracheophyta"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Macfadyena"] <-
    "Magnoliopsida"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Macfadyena"] <-
    "Lamiales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Macfadyena"] <-
    "Bignoniaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Pseudeurotium"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_6_genus == "Pseudeurotium"] <-
    "Fungi"
  inhouse_db$organism_2_phylum[inhouse_db$organism_6_genus == "Pseudeurotium"] <-
    "Ascomycota"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Pseudeurotium"] <-
    "Leotiomycetes"
  # no order
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Pseudeurotium"] <-
    "Pseudeurotiaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Ohtaekwangia"] <-
    "y"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Ohtaekwangia"] <-
    "Cytophagia"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Ohtaekwangia"] <-
    "Cytophagales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Ohtaekwangia"] <-
    "Cytophagaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Nigrospora"] <-
    "y"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Nigrospora"] <-
    "Trichosphaeriales"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Glyphium"] <-
    "y"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Glyphium"] <-
    "Patellariales"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Coleophoma"] <-
    "y"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Coleophoma"] <-
    "Helotiales"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Papulaspora"] <-
    "y"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Papulaspora"] <-
    "Sordariales"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Passeriniella"] <-
    "y"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Passeriniella"] <-
    "Pleosporales"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Phialemoniopsis"] <-
    "y"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Phialemoniopsis"] <-
    "Xylariales"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Resinicium"] <-
    "y"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Resinicium"] <-
    "Hymenochaetales"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Stilbella"] <-
    "y"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Stilbella"] <-
    "Hypocreales"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Pseudohyphozyma"] <-
    "y"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Pseudohyphozyma"] <-
    "Chrysozymaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Thyronectria"] <-
    "y"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Thyronectria"] <-
    "Hypocreales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Thyronectria"] <-
    "Nectriaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Racomitrium"] <-
    "y"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Racomitrium"] <-
    "Grimmiales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Racomitrium"] <-
    "Grimmiaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Pseudobotrytis"] <-
    "y"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Pseudobotrytis"] <-
    "Coniochaetales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Pseudobotrytis"] <-
    "Coniochaetaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Prototheca"] <-
    "y"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Prototheca"] <-
    "Chlorellales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Prototheca"] <-
    "Chlorellaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Oxyporus"] <-
    "y"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Oxyporus"] <-
    "Hymenochaetales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Oxyporus"] <-
    "Schizoporaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Leprocaulon"] <-
    "y"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Leprocaulon"] <-
    "Leprocaulales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Leprocaulon"] <-
    "Leprocaulaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Haliphthoros"] <-
    "y"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Haliphthoros"] <-
    "Lagenidiales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Haliphthoros"] <-
    "Haliphthoraceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Cora"] <-
    "y"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Cora"] <-
    "Agaricales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Cora"] <-
    "Hygrophoraceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Acaromyces"] <-
    "y"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Acaromyces"] <-
    "Exobasidiales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Acaromyces"] <-
    "Cryptobasidiaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Dichotomophthora"] <-
    "y"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Dichotomophthora"] <-
    "Dothideomycetes"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Dichotomophthora"] <-
    "Pleosporales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Dichotomophthora"] <-
    "Pleosporaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Lauriomyces"] <-
    "y"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Lauriomyces"] <-
    "Leotiomycetes"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Microcyclospora"] <-
    "y"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Microcyclospora"] <-
    "Dothideomycetes"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Microcyclospora"] <-
    "Capnodiales"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Zygosporium"] <-
    "y"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Zygosporium"] <-
    "Sordariomycetes"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Zygosporium"] <-
    "Xylariales"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Hansfordia"] <-
    "y"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Hansfordia"] <-
    "Sordariomycetes"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Hansfordia"] <-
    "Xylariales"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Veronaea"] <-
    "y"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Veronaea"] <-
    "Eurotiomycetes"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Veronaea"] <-
    "Chaetothyriales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Veronaea"] <-
    "Herpotrichiellaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Stachylidium"] <-
    "y"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Stachylidium"] <-
    "Sordariomycetes"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Stachylidium"] <-
    "Glomerellales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Stachylidium"] <-
    "Plectosphaerellaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Culicinomyces"] <-
    "y"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Culicinomyces"] <-
    "Sordariomycetes"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Culicinomyces"] <-
    "Hypocreales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Culicinomyces"] <-
    "Clavicipitaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Plectophomella"] <-
    "y"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Plectophomella"] <-
    "Dothideomycetes"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Plectophomella"] <-
    "Pleosporales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Plectophomella"] <-
    "Leptosphaeriaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Chaetosphaeronema"] <-
    "y"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Chaetosphaeronema"] <-
    "Dothideomycetes"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Chaetosphaeronema"] <-
    "Pleosporales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Chaetosphaeronema"] <-
    "Phaeosphaeriaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Phialomyces"] <-
    "y"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Phialomyces"] <-
    "Eurotiomycetes"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Phialomyces"] <-
    "Eurotiales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Phialomyces"] <-
    "Aspergillaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Ovadendron"] <-
    "y"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Ovadendron"] <-
    "Eurotiomycetes"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Ovadendron"] <-
    "Onygenales"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Pleiochaeta"] <-
    "y"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Pleiochaeta"] <-
    "Dothideomycetes"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Pleiochaeta"] <-
    "Pleosporales"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Monodictys"] <-
    "y"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Monodictys"] <-
    "Dothideomycetes"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Monodictys"] <-
    "Pleosporales"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Albophoma"] <-
    "y"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Albophoma"] <-
    "Sordariomycetes"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Acrophialophora"] <-
    "y"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Acrophialophora"] <-
    "Sordariomycetes"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Acrophialophora"] <-
    "Sordariales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Acrophialophora"] <-
    "Chaetomiaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Acrodontium"] <-
    "y"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Acrodontium"] <-
    "Dothideomycetes"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Acrodontium"] <-
    "Capnodiales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Acrodontium"] <-
    "Teratosphaeriaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Prasinococcus"] <-
    "y"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Prasinococcus"] <-
    "Prasinococcales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Prasinococcus"] <-
    "Prasinococcaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Gyrodinium"] <-
    "y"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Gyrodinium"] <-
    "Gymnodiniales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Gyrodinium"] <-
    "Gymnodiniaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Pantoneura"] <-
    "y"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Pantoneura"] <-
    "Ceramiales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Pantoneura"] <-
    "Delesseriaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Prochlorococcus"] <-
    "y"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Prochlorococcus"] <-
    "Synechococcales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Prochlorococcus"] <-
    "Prochloraceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Westiella"] <-
    "y"
  inhouse_db$organism_2_phylum[inhouse_db$organism_6_genus == "Westiella"] <-
    "Cyanobacteria"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Westiella"] <-
    "Nostocales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Westiella"] <-
    "Hapalosiphonaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Fritschiella"] <-
    "y"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Fritschiella"] <-
    "Fritschiellaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Sporobolomyces"] <-
    "y"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Sporobolomyces"] <-
    "Sporidiobolaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Rhodotorula"] <-
    "y"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Rhodotorula"] <-
    "Sporidiobolaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Pseudomuriella"] <-
    "y"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Pseudomuriella"] <-
    "Pseudomuriellaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Melanocarpus"] <-
    "y"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Melanocarpus"] <-
    "Chaetomiaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Nakazawaea"] <-
    "y"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Nakazawaea"] <-
    "Pichiaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Cyberlindnera"] <-
    "y"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Cyberlindnera"] <-
    "Phaffomycetaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Candida"] <-
    "y"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Candida"] <-
    "Debaryomycetaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Prasinoderma"] <-
    "y"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Prasinoderma"] <-
    "Prasinococcaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Plenodomus"] <-
    "y"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Plenodomus"] <-
    "Pleosporineae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Periconia"] <-
    "y"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Periconia"] <-
    "Periconiaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Mycocentrospora"] <-
    "y"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Mycocentrospora"] <-
    "Mycosphaerellaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Paraphoma"] <-
    "y"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Paraphoma"] <-
    "Pleosporineae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Clavariopsis"] <-
    "y"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Clavariopsis"] <-
    "Halosphaeriaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Berkleasmium"] <-
    "y"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Berkleasmium"] <-
    "Tubeufiaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Astrosphaeriella"] <-
    "y"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Astrosphaeriella"] <-
    "Astrosphaeriellaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Oltmannsiellopsis"] <-
    "y"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Oltmannsiellopsis"] <-
    "Oltmannsiellopsidaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Myrmecridium"] <-
    "y"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Myrmecridium"] <-
    "Myrmecridiaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Endoconidiophora"] <-
    "y"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Endoconidiophora"] <-
    "Ceratocystidaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Myrothecium"] <-
    "y"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Myrothecium"] <-
    "Stachybotryaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Sarocladium"] <-
    "y"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Sarocladium"] <-
    "Sarocladiaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Memnoniella"] <-
    "y"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Memnoniella"] <-
    "Stachybotryaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Ilyonectria"] <-
    "y"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Ilyonectria"] <-
    "Nectriaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Gliomastix"] <-
    "y"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Gliomastix"] <-
    "Bionectriaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Geosmithia"] <-
    "y"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Geosmithia"] <-
    "Bionectriaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Circinotrichum"] <-
    "y"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Circinotrichum"] <-
    "Xylariaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Trichosporiella"] <-
    "y"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Trichosporiella"] <-
    "Dermateaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Tapesia"] <-
    "y"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Tapesia"] <-
    "Dermateaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Glarea"] <-
    "y"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Glarea"] <-
    "Helotiaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Dactylaria"] <-
    "y"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Dactylaria"] <-
    "Orbiliaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Thermomyces"] <-
    "y"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Thermomyces"] <-
    "Trichocomaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Tubakia"] <-
    "y"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Tubakia"] <-
    "Melanconiellaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Stenocarpella"] <-
    "y"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Stenocarpella"] <-
    "Diaporthaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Poterioochromonas"] <-
    "y"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Poterioochromonas"] <-
    "Ochromonadaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Sirodesmium"] <-
    "y"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Sirodesmium"] <-
    "Paradictyoarthriniaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Solibacillus"] <-
    "y"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Solibacillus"] <-
    "Planococcaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Elmerina"] <-
    "y"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Elmerina"] <-
    "Aporpiaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_6_genus == "Leucocybe"] <-
    "y"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Leucocybe"] <-
    "Tricholomataceae"

  # double taxonomies -> no "y", just choose one
  # catalogue of Life as reference, then NCBI

  inhouse_db$organism_1_kingdom[inhouse_db$organism_6_genus == "Calophyllum"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organism_6_genus == "Calophyllum"] <-
    "Tracheophyta"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Calophyllum"] <-
    "Magnoliopsida"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Calophyllum"] <-
    "Malpighiales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Calophyllum"] <-
    "Calophyllaceae"

  inhouse_db$organism_1_kingdom[inhouse_db$organism_6_genus == "Heliotropium"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organism_6_genus == "Heliotropium"] <-
    "Tracheophyta"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Heliotropium"] <-
    "Magnoliopsida"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Heliotropium"] <-
    "Boraginales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Heliotropium"] <-
    "Heliotropiaceae"

  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Xiphophora"] <-
    "Phaeophyceae"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Xiphophora"] <-
    "Fucales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Xiphophora"] <-
    "Fucaceae"

  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Hymenomonas"] <-
    "Coccolithophyceae"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Hymenomonas"] <-
    "Coccosphaerales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Hymenomonas"] <-
    "Hymenomonadaceae"

  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Chondrosia"] <-
    "Chondrosiida"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Chondrosia"] <-
    "Chondrosiidae"

  inhouse_db$organism_1_kingdom[inhouse_db$organism_6_genus == "Hemerocallis"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organism_6_genus == "Hemerocallis"] <-
    "Tracheophyta"
  inhouse_db$organism_3_class[inhouse_db$organism_6_genus == "Hemerocallis"] <-
    "Magnoliopsida"
  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Hemerocallis"] <-
    "Asparagales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Hemerocallis"] <-
    "Asphodelaceae"

  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Fischerella"] <-
    "Stigonematales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Fischerella"] <-
    "Stigonemataceae"

  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Grateloupia"] <-
    "Halymeniales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Grateloupia"] <-
    "Halymeniaceae"

  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Pellia"] <-
    "Pelliales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Pellia"] <-
    "Pelliaceae"

  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Phormidium"] <-
    "Nostocales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Phormidium"] <-
    "Oscillatoriaceae"

  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Pleopsidium"] <-
    "Acarosporales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Pleopsidium"] <-
    "Acarosporaceae"

  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Pyrenochaeta"] <-
    "Pleosporales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Pyrenochaeta"] <-
    "Cucurbitariaceae"

  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Pyrenula"] <-
    "Pyrenulales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Pyrenula"] <-
    "Pyrenulaceae"

  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Synechococcus"] <-
    "Chroococcales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Synechococcus"] <-
    "Chroococcaceae"

  inhouse_db$organism_4_order[inhouse_db$organism_6_genus == "Xanthoceras"] <-
    "Sapindales"
  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Xanthoceras"] <-
    "Sapindaceae"

  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Yucca"] <-
    "Asparagaceae"

  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Viburnum"] <-
    "Adoxaceae"

  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Viscum"] <-
    "Santalaceae"

  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Glenodinium"] <-
    "Peridiniaceae"

  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Sirococcus"] <-
    "Gnomoniaceae"

  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Monostroma"] <-
    "Gomontiaceae"

  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Panaeolus"] <-
    "Bolbitiaceae"

  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Patrinia"] <-
    "Caprifoliaceae"

  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Pedicularis"] <-
    "Orobanchaceae"

  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Saccharothrix"] <-
    "Pseudonocardiaceae"

  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Spinacia"] <-
    "Amaranthaceae"

  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Brachychiton"] <-
    "Malvaceae"

  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Campylopus"] <-
    "Dicranaceae"

  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Caraipa"] <-
    "Calophyllaceae"

  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Caraipa"] <-
    "Calophyllaceae"

  inhouse_db$organism_5_family[inhouse_db$organism_6_genus == "Cochlospermum"] <-
    "Bixaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_5_family == "Scalibregmatidae"] <-
    "y"
  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Scalibregmatidae"] <-
    "Capitellida"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_5_family == "Ropalosporaceae"] <-
    "y"
  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Ropalosporaceae"] <-
    "Umbilicariales"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_5_family == "Holopodidae"] <-
    "y"
  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Holopodidae"] <-
    "Cyrtocrinida"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_5_family == "Fuscideaceae"] <-
    "y"
  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Fuscideaceae"] <-
    "Umbilicariales"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_5_family == "Eremomycetaceae"] <-
    "y"
  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Eremomycetaceae"] <-
    "Eremomycetales"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_5_family == "Chaetopteridae"] <-
    "y"
  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Chaetopteridae"] <-
    "Spionida"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_5_family == "Cephalothecaceae"] <-
    "y"
  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Cephalothecaceae"] <-
    "Sordariales"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_5_family == "Cephalodiscidae"] <-
    "y"
  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Cephalodiscidae"] <-
    "Cephalodiscida"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_5_family == "Capitellidae"] <-
    "y"
  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Capitellidae"] <-
    "Capitellida"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_5_family == "Astasiaceae"] <-
    "y"
  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Astasiaceae"] <-
    "Rhabdomonadales"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_5_family == "Arenicolidae"] <-
    "y"
  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Arenicolidae"] <-
    "Capitellida"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_5_family == "Apiosporaceae"] <-
    "y"
  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Apiosporaceae"] <-
    "Xylariales"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_5_family == "Acinetosporaceae"] <-
    "y"
  inhouse_db$organism_3_class[inhouse_db$organism_5_family == "Acinetosporaceae"] <-
    "Phaeophyceae"
  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Acinetosporaceae"] <-
    "Ectocarpales"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_5_family == "Leguminosae"] <-
    "y"
  inhouse_db$organism_5_family[inhouse_db$organism_5_family == "Leguminosae"] <-
    "Fabaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_5_family == "Cervantesiaceae"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_5_family == "Cervantesiaceae"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organism_5_family == "Cervantesiaceae"] <-
    "Tracheophyta"
  inhouse_db$organism_3_class[inhouse_db$organism_5_family == "Cervantesiaceae"] <-
    "Magnoliopsida"
  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Cervantesiaceae"] <-
    "Santalales"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_5_family == "Globulariaceae"] <-
    "y"
  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Globulariaceae"] <-
    "Lamiales"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_5_family == "Myoporaceae"] <-
    "y"
  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Myoporaceae"] <-
    "Lamiales"
  inhouse_db$organism_5_family[inhouse_db$organism_5_family == "Myoporaceae"] <-
    "Scrophulariaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_5_family == "Flacourtiaceae"] <-
    "y"
  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Flacourtiaceae"] <-
    "Malpighiales"
  inhouse_db$organism_5_family[inhouse_db$organism_5_family == "Flacourtiaceae"] <-
    "Salicaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_5_family == "Dracaenaceae"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_5_family == "Dracaenaceae"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organism_5_family == "Dracaenaceae"] <-
    "Tracheophyta"
  inhouse_db$organism_3_class[inhouse_db$organism_5_family == "Dracaenaceae"] <-
    "Magnoliopsida"
  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Dracaenaceae"] <-
    "Asparagales"
  inhouse_db$organism_5_family[inhouse_db$organism_5_family == "Dracaenaceae"] <-
    "Ruscaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_5_family == "Agavaceae"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_5_family == "Agavaceae"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organism_5_family == "Agavaceae"] <-
    "Tracheophyta"
  inhouse_db$organism_3_class[inhouse_db$organism_5_family == "Agavaceae"] <-
    "Magnoliopsida"
  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Agavaceae"] <-
    "Asparagales"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_5_family == "Ternstroemiaceae"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_5_family == "Ternstroemiaceae"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organism_5_family == "Ternstroemiaceae"] <-
    "Tracheophyta"
  inhouse_db$organism_3_class[inhouse_db$organism_5_family == "Ternstroemiaceae"] <-
    "Magnoliopsida"
  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Ternstroemiaceae"] <-
    "Ericales"
  inhouse_db$organism_5_family[inhouse_db$organism_5_family == "Ternstroemiaceae"] <-
    "Pentaphylacaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_5_family == "Aroideae"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_5_family == "Aroideae"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organism_5_family == "Aroideae"] <-
    "Tracheophyta"
  inhouse_db$organism_3_class[inhouse_db$organism_5_family == "Aroideae"] <-
    "Magnoliopsida"
  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Aroideae"] <-
    "Alismatales"
  inhouse_db$organism_5_family[inhouse_db$organism_5_family == "Aroideae"] <-
    "Araceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_5_family == "Epacridaceae"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_5_family == "Epacridaceae"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organism_5_family == "Epacridaceae"] <-
    "Tracheophyta"
  inhouse_db$organism_3_class[inhouse_db$organism_5_family == "Epacridaceae"] <-
    "Magnoliopsida"
  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Epacridaceae"] <-
    "Ericales"
  inhouse_db$organism_5_family[inhouse_db$organism_5_family == "Epacridaceae"] <-
    "Ericaceae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_5_family == "Takakiaceae"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_5_family == "Takakiaceae"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organism_5_family == "Takakiaceae"] <-
    "Tracheophyta"
  inhouse_db$organism_3_class[inhouse_db$organism_5_family == "Takakiaceae"] <-
    "Takakiopsida"
  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Takakiaceae"] <-
    "Takakiales"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_5_family == "Cyanophoraceae"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_5_family == "Cyanophoraceae"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organism_5_family == "Cyanophoraceae"] <-
    "Glaucophyta"
  inhouse_db$organism_3_class[inhouse_db$organism_5_family == "Cyanophoraceae"] <-
    "Glaucophyceae"
  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Cyanophoraceae"] <-
    "Glaucocystales"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_5_family == "Phormidiaceae"] <-
    "y"
  inhouse_db$organism_2_phylum[inhouse_db$organism_5_family == "Phormidiaceae"] <-
    "Cyanobacteria"
  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Phormidiaceae"] <-
    "Oscillatoriales"

  # double taxonomies -> no "y", just choose one
  # catalogue of Life as reference, then NCBI

  inhouse_db$organism_1_kingdom[inhouse_db$organism_5_family == "Fabaceae"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organism_5_family == "Fabaceae"] <-
    "Tracheophyta"
  inhouse_db$organism_3_class[inhouse_db$organism_5_family == "Fabaceae"] <-
    "Magnoliopsida"
  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Fabaceae"] <-
    "Fabales"

  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Daphniidae"] <-
    "Diplostraca"

  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Streptomycetaceae"] <-
    "Actinomycetales"

  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Myxotrichaceae"] <-
    "Helotiales"

  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Solanaceae"] <-
    "Solanales"

  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Botrydiaceae"] <-
    "Botrydiales"

  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Bonnemaisoniaceae"] <-
    "Nemaliales"

  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Chattonellaceae"] <-
    "Chattonellales"

  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Chlamydomonadaceae"] <-
    "Chlamydomonadales"

  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Chordariaceae"] <-
    "Ectocarpales"

  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Chromulinaceae"] <-
    "Chromulinales"

  inhouse_db$organism_2_phylum[inhouse_db$organism_5_family == "Chroomonadaceae"] <-
    "Cryptophyta"

  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Corallinaceae"] <-
    "Corallinales"

  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Cupressaceae"] <-
    "Pinales"

  inhouse_db$organism_3_class[inhouse_db$organism_5_family == "Euglenaceae"] <-
    "Euglenophyceae"
  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Euglenaceae"] <-
    "Euglenales"

  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Gracilariaceae"] <-
    "Gigartinales"

  inhouse_db$organism_3_class[inhouse_db$organism_5_family == "Halosphaeriaceae"] <-
    "Sordariomycetes"
  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Halosphaeriaceae"] <-
    "Microascales"

  inhouse_db$organism_3_class[inhouse_db$organism_5_family == "Isochrysidaceae"] <-
    "Coccolithophyceae"

  # debate
  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Kallymeniaceae"] <-
    "Gigartinales"

  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Micrococcaceae"] <-
    "Actinomycetales"

  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Mycosphaerellaceae"] <-
    "Mycosphaerellales"

  inhouse_db$organism_3_class[inhouse_db$organism_5_family == "Orbiliaceae"] <-
    "Orbiliomycetes"
  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Orbiliaceae"] <-
    "Orbiliales"

  inhouse_db$organism_3_class[inhouse_db$organism_5_family == "Paradictyoarthriniaceae"] <-
    "Dothideomycetes"
  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Paradictyoarthriniaceae"] <-
    "Pleosporales"

  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Prochloraceae"] <-
    "Synechococcales"

  inhouse_db$organism_3_class[inhouse_db$organism_5_family == "Prymnesiaceae"] <-
    "Coccolithophyceae"

  inhouse_db$organism_3_class[inhouse_db$organism_5_family == "Pseudeurotiaceae"] <-
    "Leotiomycetes"

  inhouse_db$organism_3_class[inhouse_db$organism_5_family == "Pycnococcaceae"] <-
    "Pyramimonadophyceae"

  # not even in the initial two
  inhouse_db$organism_3_class[inhouse_db$organism_5_family == "Sarcinochrysidaceae"] <-
    "Pelagophyceae"

  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Scytosiphonaceae"] <-
    "Ectocarpales"

  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Tubeufiaceae"] <-
    "Tubeufiales"

  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Tetractinellida"] <-
    "Demospongiae"

  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Dunaliellaceae"] <-
    "Chlamydomonadales"

  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Gelidiaceae"] <-
    "Nemaliales"

  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Geodiidae"] <-
    "Tetractinellida"

  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Haematococcaceae"] <-
    "Chlamydomonadales"

  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Halichondriidae"] <-
    "Suberitida"

  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Holothuriidae"] <-
    "Holothuriida"

  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Lecideaceae"] <-
    "Lecanorales"

  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Ommastrephidae"] <-
    "Oegopsida"

  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Ostreopsidaceae"] <-
    "Peridiniales"

  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Peyssonneliaceae"] <-
    "Halymeniales"

  inhouse_db$organism_3_class[inhouse_db$organism_5_family == "Phaeocystaceae"] <-
    "Coccolithophyceae"
  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Phaeocystaceae"] <-
    "Prymnesiales"

  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Plocamiaceae"] <-
    "Gigartinales"

  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Podocarpaceae"] <-
    "Pinales"

  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Rhizocarpaceae"] <-
    "Lecanorales"

  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Scenedesmaceae"] <-
    "Sphaeropleales"

  inhouse_db$organism_3_class[inhouse_db$organism_5_family == "Sparidae"] <-
    "Actinopterygii"
  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Sparidae"] <-
    "Perciformes"

  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Synechococcaceae"] <-
    "Synechococcales"

  inhouse_db$organism_3_class[inhouse_db$organism_5_family == "Ishigeaceae"] <-
    "Phaeophyceae"
  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Ishigeaceae"] <-
    "Ectocarpales"

  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Ochromonadaceae"] <-
    "Ochromonadales"

  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Xylariaceae"] <-
    "Xylariales"

  inhouse_db$organism_4_order[inhouse_db$organism_5_family == "Ancorinidae"] <-
    "Tetractinellida"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_4_order == "Thermoleophilales"] <-
    "y"
  inhouse_db$organism_3_class[inhouse_db$organism_4_order == "Thermoleophilales"] <-
    "Thermoleophilia"

  inhouse_db$organism_3_class[inhouse_db$organism_4_order == "Xylariales"] <-
    "Sordariomycetes"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_4_order == "Rubrobacterales"] <-
    "y"
  inhouse_db$organism_3_class[inhouse_db$organism_4_order == "Rubrobacterales"] <-
    "Rubrobacteria"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_4_order == "Umbelopsidales"] <-
    "y"
  inhouse_db$organism_3_class[inhouse_db$organism_4_order == "Umbelopsidales"] <-
    "Umbelopsidomycetes"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_4_order == "Tricladida"] <-
    "y"
  inhouse_db$organism_3_class[inhouse_db$organism_4_order == "Tricladida"] <-
    "Rhabditophora"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_4_order == "Polycladida"] <-
    "y"
  inhouse_db$organism_3_class[inhouse_db$organism_4_order == "Polycladida"] <-
    "Rhabditophora"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_4_order == "Mycoplasmatales"] <-
    "y"
  inhouse_db$organism_3_class[inhouse_db$organism_4_order == "Mycoplasmatales"] <-
    "Mollicutes"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_4_order == "Bifidobacteriales"] <-
    "y"
  inhouse_db$organism_3_class[inhouse_db$organism_4_order == "Bifidobacteriales"] <-
    "Actinobacteria"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_4_order == "Coriobacteriales"] <-
    "y"
  inhouse_db$organism_3_class[inhouse_db$organism_4_order == "Coriobacteriales"] <-
    "Coriobacteriia"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_4_order == "Cryptonemiales"] <-
    "y"
  inhouse_db$organism_3_class[inhouse_db$organism_4_order == "Cryptonemiales"] <-
    "Florideophyceae"
  inhouse_db$organism_4_order[inhouse_db$organism_4_order == "Cryptonemiales"] <-
    "Halymeniales"

  # double taxonomies -> no "y", just choose one
  # catalogue of Life as reference, then NCBI

  inhouse_db$organism_3_class[inhouse_db$organism_4_order == "Cladophorales"] <-
    "Ulvophyceae"

  inhouse_db$organism_3_class[inhouse_db$organism_4_order == "Isochrysidales"] <-
    "Coccolithophyceae"

  inhouse_db$organism_1_kingdom[inhouse_db$organism_4_order == "Alcyonacea"] <-
    "Animalia"
  inhouse_db$organism_2_phylum[inhouse_db$organism_4_order == "Alcyonacea"] <-
    "Cnidaria"
  inhouse_db$organism_3_class[inhouse_db$organism_4_order == "Alcyonacea"] <-
    "Anthozoa"

  inhouse_db$organism_3_class[inhouse_db$organism_4_order == "Actinomycetales"] <-
    "Actinobacteria"

  inhouse_db$organism_3_class[inhouse_db$organism_4_order == "Bryopsidales"] <-
    "Ulvophyceae"

  inhouse_db$organism_3_class[inhouse_db$organism_4_order == "Chlorodendrales"] <-
    "Chlorodendrophyceae"

  inhouse_db$organism_3_class[inhouse_db$organism_4_order == "Cypriniformes"] <-
    "Actinopterygii"

  inhouse_db$organism_3_class[inhouse_db$organism_4_order == "Synechococcales"] <-
    "Cyanophyceae"

  inhouse_db$organism_3_class[inhouse_db$organism_4_order == "Patellariales"] <-
    "Dothideomycetes"

  inhouse_db$organism_3_class[inhouse_db$organism_4_order == "Pyramimonadales"] <-
    "Pyramimonadophyceae"

  inhouse_db$organism_3_class[inhouse_db$organism_4_order == "Suberitida"] <-
    "Demospongiae"

  inhouse_db$organism_3_class[inhouse_db$organism_4_order == "Tetractinellida"] <-
    "Demospongiae"

  inhouse_db$organism_3_class[inhouse_db$organism_4_order == "Gigartinales"] <-
    "Florideophyceae"

  inhouse_db$organism_3_class[inhouse_db$organism_4_order == "Corallinales"] <-
    "Florideophyceae"

  inhouse_db$organism_3_class[inhouse_db$organism_4_order == "Halymeniales"] <-
    "Florideophyceae"

  inhouse_db$organism_3_class[inhouse_db$organism_4_order == "Chlorellales"] <-
    "Trebouxiophyceae"

  inhouse_db$organism_3_class[inhouse_db$organism_4_order == "Sarcinochrysidales"] <-
    "Pelagophyceae"

  inhouse_db$organism_3_class[inhouse_db$organism_4_order == "Ceramiales"] <-
    "Florideophyceae"

  inhouse_db$organism_3_class[inhouse_db$organism_4_order == "Chroococcales"] <-
    "Cyanophyceae"

  inhouse_db$organism_3_class[inhouse_db$organism_4_order == "Nostocales"] <-
    "Cyanophyceae"

  inhouse_db$organism_2_phylum[inhouse_db$organism_3_class == "Cryptophyceae"] <-
    "Cryptophyta"

  inhouse_db$organism_2_phylum[inhouse_db$organism_3_class == "Florideophyceae"] <-
    "Rhodophyta"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_3_class == "Phaeophyceae"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_3_class == "Phaeophyceae"] <-
    "Chromista"
  inhouse_db$organism_2_phylum[inhouse_db$organism_3_class == "Phaeophyceae"] <-
    "Ochrophyta"

  # inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_3_class == "Flavobacteriia"] <- "y"
  inhouse_db$organism_3_class[inhouse_db$organism_3_class == "Flavobacteriia"] <-
    "Flavobacteria"

  # Because of NCBI: inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_3_class == "Magnoliopsida"] <-
  #   "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_3_class == "Magnoliopsida"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organism_3_class == "Magnoliopsida"] <-
    "Tracheophyta"

  # Because of NCBI: inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_3_class == "Magnoliopsida"] <-
  #   "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_3_class == "Magnoliopsida"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organism_3_class == "Magnoliopsida"] <-
    "Tracheophyta"

  # Because of NCBI: inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_3_class == "Cyanobacteria"] <-
  #   "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_3_class == "Cyanobacteria"] <-
    "Bacteria"
  inhouse_db$organism_2_phylum[inhouse_db$organism_3_class == "Cyanobacteria"] <-
    "Cyanophyceae"

  # double taxonomies -> no "y", just choose one
  # catalogue of Life as reference, then NCBI

  inhouse_db$organism_2_phylum[inhouse_db$organism_3_class == "Dinophyceae"] <-
    "Myzozoa"

  inhouse_db$organism_2_phylum[inhouse_db$organism_3_class == "Bryopsida"] <-
    "Bryophyta"

  inhouse_db$organism_2_phylum[inhouse_db$organism_3_class == "Sphagnopsida"] <-
    "Bryophyta"

  inhouse_db$organism_2_phylum[inhouse_db$organism_3_class == "Cycadopsida"] <-
    "Tracheophyta"

  inhouse_db$organism_2_phylum[inhouse_db$organism_3_class == "Cycadopsida"] <-
    "Tracheophyta"

  inhouse_db$organism_1_kingdom[inhouse_db$organism_3_class == "Raphidophyceae"] <-
    "Chromista"
  inhouse_db$organism_2_phylum[inhouse_db$organism_3_class == "Raphidophyceae"] <-
    "Ochrophyta"

  inhouse_db$organism_1_kingdom[inhouse_db$organism_3_class == "Pelagophyceae"] <-
    "Chromista"
  inhouse_db$organism_2_phylum[inhouse_db$organism_3_class == "Pelagophyceae"] <-
    "Ochrophyta"

  # inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_3_class == "Liliopsida"] <- "y"
  inhouse_db$organism_3_class[inhouse_db$organism_3_class == "Liliopsida"] <-
    "Magnoliopsida"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_2_phylum == "Evosea"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_2_phylum == "Evosea"] <-
    "Protozoa"

  # Because of NCBI: inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_2_phylum == "Rhodophyta"] <-
  #   "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_2_phylum == "Rhodophyta"] <-
    "Plantae"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_2_phylum == "Bacillariophyta"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_2_phylum == "Bacillariophyta"] <-
    "Chromista"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_2_phylum == "Tenericutes"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_2_phylum == "Tenericutes"] <-
    "Bacteria"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_2_phylum == "Dinophyta"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_2_phylum == "Dinophyta"] <-
    "Protozoa"

  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_2_phylum == "Sipuncula"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_2_phylum == "Sipuncula"] <-
    "Animalia"

  # inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_2_phylum == "Miozoa"] <- "y"
  inhouse_db$organism_2_phylum[inhouse_db$organism_2_phylum == "Miozoa"] <-
    "Myzozoa"

  # inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_2_phylum == "Myzozoa"] <- "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_2_phylum == "Myzozoa"] <-
    "Chromista"

  # inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_2_phylum == "Myzozoa"] <- "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_2_phylum == "Myzozoa"] <-
    "Chromista"

  # inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_2_phylum == "Myzozoa"] <- "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_2_phylum == "Myzozoa"] <-
    "Chromista"

  # inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_2_phylum == "Chlorophyta"] <- "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_2_phylum == "Chlorophyta"] <-
    "Plantae"

  # inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_2_phylum == "Ciliophora"] <- "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_2_phylum == "Ciliophora"] <-
    "Chromista"

  # inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_2_phylum == "Cyanobacteria"] <- "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_2_phylum == "Cyanobacteria"] <-
    "Bacteria"

  # inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_2_phylum == "Glaucophyta"] <- "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_2_phylum == "Glaucophyta"] <-
    "Plantae"

  # inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_2_phylum == "Ochrophyta"] <- "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_2_phylum == "Ochrophyta"] <-
    "Chromista"

  inhouse_db$organism_1_kingdom[inhouse_db$organism_2_phylum == "Cryptophyta"] <-
    "Chromista"

  # inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_2_phylum == "Tracheophyta"] <- "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_2_phylum == "Tracheophyta"] <-
    "Plantae"

  # inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_1_kingdom == "Protista"] <- "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_1_kingdom == "Protista"] <-
    "Protozoa"

  inhouse_db$organism_1_kingdom <-
    y_as_na(
      x = inhouse_db$organism_1_kingdom,
      y = ""
    )

  inhouse_db$organism_2_phylum <-
    y_as_na(
      x = inhouse_db$organism_2_phylum,
      y = ""
    )

  inhouse_db$organism_3_class <-
    y_as_na(
      x = inhouse_db$organism_3_class,
      y = ""
    )

  inhouse_db$organism_4_order <-
    y_as_na(
      x = inhouse_db$organism_4_order,
      y = ""
    )

  inhouse_db$organism_5_family <-
    y_as_na(
      x = inhouse_db$organism_5_family,
      y = ""
    )

  inhouse_db$organism_6_genus <-
    y_as_na(
      x = inhouse_db$organism_6_genus,
      y = ""
    )

  inhouse_db$organism_7_species <-
    y_as_na(
      x = inhouse_db$organism_7_species,
      y = ""
    )

  inhouse_db$organism_8_variety <-
    y_as_na(
      x = inhouse_db$organism_8_variety,
      y = ""
    )

  inhouse_db$organism_modified_taxonomy_manual <-
    y_as_na(
      x = inhouse_db$organism_modified_taxonomy_manual,
      y = ""
    )

  inhouse_db$organismCurated <-
    as.character(apply(inhouse_db[6:13], 1, function(x) {
      tail(na.omit(x), 1)
    }))

  organism_8_variety_cleaning <-
    inhouse_db %>%
    distinct(
      organism_1_kingdom,
      organism_2_phylum,
      organism_3_class,
      organism_4_order,
      organism_5_family,
      organism_6_genus,
      organism_7_species,
      organism_8_variety,
      .keep_all = TRUE
    ) %>%
    group_by(organism_8_variety) %>%
    filter(!is.na(organism_8_variety)) %>%
    add_count() %>%
    filter(n >= 2) %>%
    select(-n)

  cat(
    "you have",
    nrow(organism_8_variety_cleaning),
    "species with inconsistent upstream taxonomies",
    "\n"
  )

  inhouse_db_old <- inhouse_db %>%
    filter(!organismCurated %in% organism_8_variety_cleaning$organismCurated)

  inhouse_db_new <- inhouse_db %>%
    filter(organismCurated %in% organism_8_variety_cleaning$organismCurated) %>%
    select(
      organismOriginal,
      organismCleaned,
      organismCurated
    )

  if (nrow(inhouse_db_new) > 0) {
    inhouse_db_new <-
      left_join(
        inhouse_db_new,
        organism_8_variety_cleaning
      )
  }

  inhouse_db_organism_8_variety_clean <-
    rbind(inhouse_db_old, inhouse_db_new)

  organism_7_species_cleaning <-
    inhouse_db_organism_8_variety_clean %>%
    distinct(
      organism_1_kingdom,
      organism_2_phylum,
      organism_3_class,
      organism_4_order,
      organism_5_family,
      organism_6_genus,
      organism_7_species,
      .keep_all = TRUE
    ) %>%
    group_by(organism_7_species) %>%
    filter(!is.na(organism_7_species)) %>%
    add_count() %>%
    filter(n >= 2) %>%
    select(-n)

  cat(
    "you have",
    nrow(organism_7_species_cleaning),
    "species with inconsistent upstream taxonomies",
    "\n"
  )

  inhouse_db_old <- inhouse_db %>%
    filter(!organismCurated %in% organism_7_species_cleaning$organismCurated)

  inhouse_db_new <- inhouse_db %>%
    filter(organismCurated %in% organism_7_species_cleaning$organismCurated) %>%
    select(
      organismOriginal,
      organismCleaned,
      organismCurated
    )

  if (nrow(inhouse_db_new) > 0) {
    inhouse_db_new <-
      left_join(
        inhouse_db_new,
        organism_7_species_cleaning
      )
  }

  inhouse_db_organism_7_species_clean <-
    rbind(inhouse_db_old, inhouse_db_new)

  organism_6_genus_cleaning <-
    inhouse_db_organism_7_species_clean %>%
    distinct(
      organism_1_kingdom,
      organism_2_phylum,
      organism_3_class,
      organism_4_order,
      organism_5_family,
      organism_6_genus,
      .keep_all = TRUE
    ) %>%
    group_by(organism_6_genus) %>%
    filter(!is.na(organism_6_genus)) %>%
    add_count() %>%
    filter(n >= 2) %>%
    select(-n)

  cat(
    "you have",
    nrow(organism_6_genus_cleaning),
    "genera with inconsistent upstream taxonomies",
    "\n"
  )

  inhouse_db_old <- inhouse_db_organism_7_species_clean %>%
    filter(!organismCurated %in% organism_6_genus_cleaning$organismCurated)

  inhouse_db_new <- inhouse_db_organism_7_species_clean %>%
    filter(organismCurated %in% organism_6_genus_cleaning$organismCurated) %>%
    select(
      organismOriginal,
      organismCleaned,
      organismCurated
    )

  if (nrow(inhouse_db_new) > 0) {
    inhouse_db_new <-
      left_join(
        inhouse_db_new,
        organism_6_genus_cleaning
      )
  }

  inhouse_db_organism_6_genus_clean <-
    rbind(inhouse_db_old, inhouse_db_new)

  organism_5_family_cleaning <-
    inhouse_db_organism_6_genus_clean %>%
    distinct(
      organism_1_kingdom,
      organism_2_phylum,
      organism_3_class,
      organism_4_order,
      organism_5_family,
      .keep_all = TRUE
    ) %>%
    group_by(organism_5_family) %>%
    filter(!is.na(organism_5_family)) %>%
    add_count() %>%
    filter(n >= 2) %>%
    select(-n)

  cat(
    "you have",
    nrow(organism_5_family_cleaning),
    "families with inconsistent upstream taxonomies",
    "\n"
  )

  inhouse_db_old <- inhouse_db_organism_6_genus_clean %>%
    filter(!organismCurated %in% organism_5_family_cleaning$organismCurated)

  inhouse_db_new <- inhouse_db_organism_6_genus_clean %>%
    filter(organismCurated %in% organism_5_family_cleaning$organismCurated) %>%
    select(
      organismOriginal,
      organismCleaned,
      organismCurated
    )

  if (nrow(inhouse_db_new) > 0) {
    inhouse_db_new <-
      left_join(
        inhouse_db_new,
        organism_5_family_cleaning
      )
  }

  inhouse_db_organism_5_family_clean <-
    rbind(inhouse_db_old, inhouse_db_new)

  organism_4_order_cleaning <-
    inhouse_db_organism_5_family_clean %>%
    distinct(
      organism_1_kingdom,
      organism_2_phylum,
      organism_3_class,
      organism_4_order,
      .keep_all = TRUE
    ) %>%
    group_by(organism_4_order) %>%
    filter(!is.na(organism_4_order)) %>%
    add_count() %>%
    filter(n >= 2) %>%
    select(-n)

  cat(
    "you have",
    nrow(organism_4_order_cleaning),
    "orders with inconsistent upstream taxonomies",
    "\n"
  )

  inhouse_db_old <- inhouse_db_organism_5_family_clean %>%
    filter(!organismCurated %in% organism_4_order_cleaning$organismCurated)

  inhouse_db_new <- inhouse_db_organism_5_family_clean %>%
    filter(organismCurated %in% organism_4_order_cleaning$organismCurated) %>%
    select(
      organismOriginal,
      organismCleaned,
      organismCurated
    )

  if (nrow(inhouse_db_new) > 0) {
    inhouse_db_new <-
      left_join(
        inhouse_db_new,
        organism_4_order_cleaning
      )
  }

  inhouse_db_organism_4_order_clean <-
    rbind(inhouse_db_old, inhouse_db_new)

  organism_3_class_cleaning <-
    inhouse_db_organism_4_order_clean %>%
    distinct(organism_1_kingdom,
      organism_2_phylum,
      organism_3_class,
      .keep_all = TRUE
    ) %>%
    group_by(organism_3_class) %>%
    filter(!is.na(organism_3_class)) %>%
    add_count() %>%
    filter(n >= 2) %>%
    select(-n)

  cat(
    "you have",
    nrow(organism_3_class_cleaning),
    "classes with inconsistent upstream taxonomies",
    "\n"
  )

  inhouse_db_old <- inhouse_db_organism_4_order_clean %>%
    filter(!organismCurated %in% organism_3_class_cleaning$organismCurated)

  inhouse_db_new <- inhouse_db_organism_4_order_clean %>%
    filter(organismCurated %in% organism_3_class_cleaning$organismCurated) %>%
    select(
      organismOriginal,
      organismCleaned,
      organismCurated
    )

  if (nrow(inhouse_db_new) > 0) {
    inhouse_db_new <-
      left_join(inhouse_db_new, organism_3_class_cleaning)
  }

  inhouse_db_organism_3_class_clean <-
    rbind(inhouse_db_old, inhouse_db_new)

  organism_2_phylum_cleaning <-
    inhouse_db_organism_3_class_clean %>%
    distinct(organism_1_kingdom,
      organism_2_phylum,
      .keep_all = TRUE
    ) %>%
    group_by(organism_2_phylum) %>%
    filter(!is.na(organism_2_phylum)) %>%
    add_count() %>%
    filter(n >= 2) %>%
    select(-n)

  cat(
    "you have",
    nrow(organism_2_phylum_cleaning),
    "phyla with inconsistent upstream taxonomies",
    "\n"
  )

  inhouse_db_old <- inhouse_db_organism_3_class_clean %>%
    filter(!organismCurated %in% organism_2_phylum_cleaning$organismCurated)

  inhouse_db_new <- inhouse_db_organism_3_class_clean %>%
    filter(organismCurated %in% organism_2_phylum_cleaning$organismCurated) %>%
    select(
      organismOriginal,
      organismCleaned,
      organismCurated
    )

  if (nrow(inhouse_db_new) > 0) {
    inhouse_db_new <-
      left_join(
        inhouse_db_new,
        organism_2_phylum_cleaning
      )
  }

  inhouse_db_organism_2_phylum_clean <-
    rbind(inhouse_db_old, inhouse_db_new)

  organism_1_kingdom_cleaning <-
    inhouse_db_organism_2_phylum_clean %>%
    distinct(organism_1_kingdom,
      .keep_all = TRUE
    )

  organism_1_kingdom_cleaning_2 <-
    inhouse_db_organism_2_phylum_clean %>%
    filter(!is.na(organism_1_kingdom)) %>%
    distinct(organism_1_kingdom,
      .keep_all = TRUE
    )

  cat(
    "you have",
    nrow(organism_1_kingdom_cleaning_2),
    "different kingdoms represented",
    "\n"
  )

  inhouse_db_old <- inhouse_db_organism_2_phylum_clean %>%
    filter(!organismCurated %in% organism_1_kingdom_cleaning$organismCurated)

  inhouse_db_new <- inhouse_db_organism_2_phylum_clean %>%
    filter(organismCurated %in% organism_1_kingdom_cleaning$organismCurated)

  if (nrow(inhouse_db_new) > 0) {
    inhouse_db_new <-
      left_join(inhouse_db_new, organism_1_kingdom_cleaning)
  }

  inhouse_db_organism_1_kingdom_clean <-
    rbind(inhouse_db_old, inhouse_db_new)

  inhouse_db_organism_1_kingdom_clean$organism_1_kingdom <-
    gsub(
      pattern = "Not assigned",
      replacement = "",
      x = inhouse_db_organism_1_kingdom_clean$organism_1_kingdom
    )

  inhouse_db_organism_1_kingdom_clean$organism_2_phylum <-
    gsub(
      pattern = "Not assigned",
      replacement = "",
      x = inhouse_db_organism_1_kingdom_clean$organism_2_phylum
    )

  inhouse_db_organism_1_kingdom_clean$organism_3_class <-
    gsub(
      pattern = "Not assigned",
      replacement = "",
      x = inhouse_db_organism_1_kingdom_clean$organism_3_class
    )

  inhouse_db_organism_1_kingdom_clean$organism_4_order <-
    gsub(
      pattern = "Not assigned",
      replacement = "",
      x = inhouse_db_organism_1_kingdom_clean$organism_4_order
    )

  inhouse_db_organism_1_kingdom_clean$organism_5_family <-
    gsub(
      pattern = "Not assigned",
      replacement = "",
      x = inhouse_db_organism_1_kingdom_clean$organism_5_family
    )

  inhouse_db_organism_1_kingdom_clean$organism_6_genus <-
    gsub(
      pattern = "Not assigned",
      replacement = "",
      x = inhouse_db_organism_1_kingdom_clean$organism_6_genus
    )

  inhouse_db_organism_1_kingdom_clean$organism_6_genus <-
    gsub(
      pattern = "Fructus",
      replacement = "",
      x = inhouse_db_organism_1_kingdom_clean$organism_6_genus
    )

  inhouse_db_organism_1_kingdom_clean$organism_7_species <-
    gsub(
      pattern = "Not assigned",
      replacement = "",
      x = inhouse_db_organism_1_kingdom_clean$organism_7_species
    )

  inhouse_db_organism_1_kingdom_clean$organism_7_species <-
    gsub(
      pattern = "unidentified",
      replacement = "",
      x = inhouse_db_organism_1_kingdom_clean$organism_7_species
    )

  inhouse_db_organism_1_kingdom_clean$organismCurated <-
    gsub(
      pattern = "unidentified",
      replacement = "",
      x = inhouse_db_organism_1_kingdom_clean$organismCleaned
    )

  inhouse_db_organism_1_kingdom_clean$organismCurated <-
    gsub(
      pattern = "Fructus",
      replacement = "",
      x = inhouse_db_organism_1_kingdom_clean$organismCleaned
    )

  inhouse_db_organism_1_kingdom_clean$organismCurated <-
    gsub(
      pattern = "Radix",
      replacement = "",
      x = inhouse_db_organism_1_kingdom_clean$organismCleaned
    )

  # to avoid false genera
  inhouse_db_organism_1_kingdom_clean$organism_7_species[str_count(
    string = inhouse_db_organism_1_kingdom_clean$organism_7_species,
    pattern = "\\w+"
  ) == 1] <-
    ""

  inhouse_db_organism_1_kingdom_clean$organism_1_kingdom <-
    y_as_na(
      x = inhouse_db_organism_1_kingdom_clean$organism_1_kingdom,
      y = ""
    )

  inhouse_db_organism_1_kingdom_clean$organism_2_phylum <-
    y_as_na(
      x = inhouse_db_organism_1_kingdom_clean$organism_2_phylum,
      y = ""
    )

  inhouse_db_organism_1_kingdom_clean$organism_3_class <-
    y_as_na(
      x = inhouse_db_organism_1_kingdom_clean$organism_3_class,
      y = ""
    )

  inhouse_db_organism_1_kingdom_clean$organism_4_order <-
    y_as_na(
      x = inhouse_db_organism_1_kingdom_clean$organism_4_order,
      y = ""
    )

  inhouse_db_organism_1_kingdom_clean$organism_5_family <-
    y_as_na(
      x = inhouse_db_organism_1_kingdom_clean$organism_5_family,
      y = ""
    )

  inhouse_db_organism_1_kingdom_clean$organism_6_genus <-
    y_as_na(
      x = inhouse_db_organism_1_kingdom_clean$organism_6_genus,
      y = ""
    )

  inhouse_db_organism_1_kingdom_clean$organism_7_species <-
    y_as_na(
      x = inhouse_db_organism_1_kingdom_clean$organism_7_species,
      y = ""
    )

  inhouse_db_organism_1_kingdom_clean$organism_8_variety <-
    y_as_na(
      x = inhouse_db_organism_1_kingdom_clean$organism_8_variety,
      y = ""
    )

  inhouse_db_organism_1_kingdom_clean$organismCurated <-
    y_as_na(
      x = inhouse_db_organism_1_kingdom_clean$organismCurated,
      y = ""
    )

  inhouse_db_organism_1_kingdom_clean$organism_modified_taxonomy_manual <-
    y_as_na(
      x = inhouse_db_organism_1_kingdom_clean$organism_modified_taxonomy_manual,
      y = ""
    )

  inhouse_db_organism_1_kingdom_clean$organismCurated <-
    as.character(apply(inhouse_db_organism_1_kingdom_clean[6:13], 1, function(x) {
      tail(na.omit(x), 1)
    }))

  inhouse_db_organism_1_kingdom_clean <-
    inhouse_db_organism_1_kingdom_clean %>%
    mutate_all(as.character)

  inhouse_db_organism_1_kingdom_clean$organismCurated <-
    y_as_na(
      x = inhouse_db_organism_1_kingdom_clean$organismCurated,
      y = "character(0)"
    )

  inhouse_db_organism_1_kingdom_clean$organismCurated <-
    y_as_na(
      x = inhouse_db_organism_1_kingdom_clean$organismCurated,
      y = "NA"
    )

  return(inhouse_db_organism_1_kingdom_clean)
}