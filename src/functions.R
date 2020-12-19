#######################################################
####################   Functions   ####################
#######################################################

library(ChemmineR)
library(chorddiag)
library(collapsibleTree)
library(data.table)
library(ggalluvial)
library(ggfittext)
library(ggraph)
library(Hmisc)
library(jsonlite)
library(parallel)
library(pbmcapply)
library(plotly)
library(rcrossref)
library(readxl)
library(rentrez)
library(rvest)
library(splitstackshape)
library(stringdist)
library(stringi)
library(UpSetR)
library(webchem)
library(XML)
library(tidyverse)

#######################################################
#' Title
#'
#' @param path_to_db
#'
#' @return
#' @export
#'
#' @examples
db_loader <- function(path_to_db) {
  db <- read_delim(
    file = gzfile(path_to_db),
    col_types = cols(.default = "c"),
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  )
  return(db)
}

#######################################################

#' Title
#'
#' @param i
#'
#' @return
#' @export
#'
#' @examples
pubchem2inchi <- function(i) {
  tryCatch(
    {
      cpd <-
        data_translated_pubchem[i, "structure_original_numerical_pubchem"]
      url <-
        paste(
          "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
          cpd,
          "/property/InChI/txt",
          sep = ""
        )
      url <- gsub(
        pattern = "\\s",
        replacement = "%20",
        x = url
      )
      read_html(url) %>%
        html_text()
    },
    error = function(e) {
      NA
    }
  )
}

#######################################################

#' Title
#'
#' @param x
#' @param n
#'
#' @return
#' @export
#'
#' @examples
shift <- function(x, n) {
  c(x[-(seq(n))], rep(NA, n))
}

#######################################################

#' Title
#'
#' @param x
#' @param y
#'
#' @return
#' @export
#'
#' @examples
biocleaning <- function(x, y) {
  # selecting and adding row number
  df <- x %>%
    select(
      names.verbatim,
      names.verification
    ) %>%
    data.table() %>%
    select(
      names.verbatim,
      names.verification.bestResult.dataSourceTitle,
      names.verification.bestResult.matchedCanonicalSimple,
      names.verification.bestResult.classificationRank,
      names.verification.bestResult.classificationPath,
      names.verification.bestResult.isSynonym
    ) %>%
    mutate(nrow = row_number())

  # extracting preferred results data table
  ## as list of dataframes
  df2 <- x$names.verification$preferredResults

  # outputting row numbers
  rows <- df2 %>%
    data.table() %>%
    mutate(nrow = row_number()) %>%
    filter(. != "NULL") %>%
    select(nrow)

  ## as dataframe and adding row number
  df3 <- bind_rows(df2,
    .id = "id"
  )

  # selecting best result (with best score and best filled taxonomy)
  df4 <- df3 %>%
    rowwise() %>%
    mutate(
      kingdom = sum(as.numeric(
        stri_detect(
          str = classificationRank,
          fixed = c(
            "kingdom",
            "Kingdom",
            "regn."
          )
        )
      )),
      phylum = sum(as.numeric(
        stri_detect(
          str = classificationRank,
          fixed = c(
            "phylum",
            "Phylum",
            "phyl."
          )
        )
      )),
      class = sum(as.numeric(
        stri_detect(
          str = classificationRank,
          fixed = c(
            "class",
            "Class",
            "cl."
          )
        )
      )),
      order = sum(as.numeric(
        stri_detect(
          str = classificationRank,
          fixed = c(
            "order",
            "Order",
            "ord."
          )
        )
      )),
      family = sum(as.numeric(
        stri_detect(
          str = classificationRank,
          fixed = c(
            "family",
            "Family",
            "fam."
          )
        )
      )),
      genus = sum(as.numeric(
        stri_detect(
          str = classificationRank,
          fixed = c(
            "genus",
            "Genus"
          )
        )
      )),
      species = sum(as.numeric(
        stri_detect(
          str = classificationRank,
          fixed = c(
            "species",
            "Species",
            "spec.",
            "sp."
          )
        )
      ))
    ) %>%
    ungroup()

  df4$kingdom[df4$kingdom >= 1] <- 1
  df4$phylum[df4$phylum >= 1] <- 1
  df4$class[df4$class >= 1] <- 1
  df4$order[df4$order >= 1] <- 1
  df4$family[df4$family >= 1] <- 1
  df4$genus[df4$genus >= 1] <- 1
  df4$species[df4$species >= 1] <- 1

  # the synonym part is there to avoid the (actually)
  ## non-optimal output from Catalogue of Life in GNFinder
  ### (explained in https://github.com/gnames/gnfinder/issues/48)
  df5 <- df4 %>%
    mutate(n = rowSums(.[c(
      "kingdom",
      "phylum",
      "class",
      "order",
      "family",
      "genus",
      "species"
    )])) %>%
    group_by(id) %>%
    arrange(desc(n), !is.na(isSynonym)) %>%
    ungroup() %>%
    distinct(id,
      .keep_all = TRUE
    ) %>%
    arrange(as.numeric(id))

  df6 <- cbind(df5, rows)

  # adding row number
  df7 <- x$names.start %>%
    data.table() %>%
    mutate(nrow = row_number())

  colnames(df7)[1] <- "sum"

  # joining
  taxo <- right_join(df6, df7) %>%
    select(
      canonicalname = matchedCanonicalSimple,
      db_taxo = dataSourceTitle,
      taxonomy = classificationPath,
      rank = classificationRank,
      sum
    )

  # computing sum of characters to match with GNFinder results
  y$nchar <-
    nchar(x = y$organismTranslated)
  y[1, "sum"] <- nchar(colnames(y)[1]) + 1
  for (i in 2:nrow(y)) {
    y[i, "sum"] <- y[i - 1, "nchar"] + 1 + y[i - 1, "sum"]
  }

  # adding min and max to merge
  taxo <- taxo %>%
    mutate(
      value_min = sum,
      value_max = sum
    ) %>%
    data.table()

  # filtering non-empty values
  y_2 <- y %>%
    filter(!is.na(organismTranslated)) %>%
    mutate(value_min = sum)

  # filling sum values
  y_2$value_min <- as.numeric(y_2$value_min)
  y_2$value_max <- shift(y_2$sum, 1) - 1
  y_2[nrow(y_2), 5] <- y_2[nrow(y_2), 4] + 10000

  # transforming as data table (needed for next function)
  y_2 <- y_2 %>%
    data.table()

  # setting joining keys
  setkey(taxo, value_min, value_max)
  setkey(y_2, value_min, value_max)

  # joining
  pre_final_db <- foverlaps(
    taxo,
    y_2
  )

  # selecting
  final_db <- left_join(
    y,
    pre_final_db
  ) %>%
    select( # -namesverbatim,-nchar,-sum,-value_max,-value_min,-i.sum,-i.value_max,
      -i.value_min
    )

  return(final_db)
}

#######################################################

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
  a <- paste("\\b", dic$taxaLevel, "\\b", sep = "")
  b <- dic$taxaLevelStandard

  # selecting and splitting taxonomy and ranks
  df1 <- dfsel %>%
    select(
      identifier = 1,
      canonicalname,
      db_taxo,
      taxonomy,
      rank
    ) %>%
    cSplit(
      splitCols = "taxonomy",
      sep = "|"
    ) %>%
    cSplit(
      splitCols = "rank",
      sep = "|"
    ) %>%
    lapply(as.character) %>%
    as_tibble()

  # manipulating taxa
  df2 <- df1 %>%
    pivot_longer(
      cols = 4:ncol(.),
      names_to = c(".value", "level"),
      names_sep = "_",
      values_to = "taxonomy",
      values_drop_na = TRUE
    ) %>%
    distinct(identifier,
      canonicalname,
      level,
      .keep_all = TRUE
    )

  df2$rank <- stri_replace_all_regex(
    str = df2$rank,
    pattern = a,
    replacement = b,
    case_insensitive = FALSE,
    vectorize_all = FALSE
  )

  # removing false non-empty cells
  df2$identifier <- y_as_na(
    x = df2$identifier,
    y = ""
  )

  df2$rank <- y_as_na(
    x = df2$rank,
    y = ""
  )

  df2$taxonomy <- y_as_na(
    x = df2$taxonomy,
    y = ""
  )

  df2$identifier <- y_as_na(
    x = df2$identifier,
    y = "NA NA"
  )

  df2$rank <- y_as_na(
    x = df2$rank,
    y = ""
  )

  df2$rank <- ifelse(test = is.na(df2$rank),
    yes = "NA",
    no = df2$rank
  )

  colnames(df2)[3] <- "db_taxo"

  # manipulating taxa
  df3 <- df2 %>%
    pivot_wider(
      names_from = rank,
      values_from = taxonomy
    ) %>%
    select_if(
      names(.) %in%
        c(
          "identifier",
          "canonicalname",
          "db_taxo",
          "kingdom",
          "phylum",
          "class",
          "order",
          "family",
          "genus",
          "species"
        )
    )

  # pasting suffix to colnames to pivot then (the double pivot allows to tidy the data)
  colnames(df3)[4:ncol(df3)] <-
    paste("bio_", colnames(df3)[4:ncol(df3)], sep = "")

  # pivoting (long)
  df4 <- df3 %>%
    pivot_longer(
      cols = 4:ncol(.),
      names_to = c(".value", "level"),
      names_sep = "_",
      values_to = "taxonomy",
      values_drop_na = TRUE
    )

  # pivoting (wide)
  df5 <- df4 %>%
    group_by(canonicalname) %>%
    distinct(level, .keep_all = TRUE) %>%
    pivot_wider(
      names_from = level,
      values_from = bio
    ) %>%
    select_if(
      names(.) %in%
        c(
          "canonicalname",
          "db_taxo",
          "kingdom",
          "phylum",
          "class",
          "order",
          "family",
          "genus",
          "species"
        )
    )

  # adding taxa to initial df
  df6 <- left_join(dfsel, df5)

  return(df6)
}


#' Title
#'
#' @param x
#' @param y
#'
#' @return
#' @export
#'
#' @examples
y_as_na <- function(x, y) {
  if ("factor" %in% class(x)) {
    x <- as.character(x)
  } ## since ifelse wont work with factors
  ifelse(test = as.character(x) != y,
    yes = x,
    no = NA
  )
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
distinct_triplets <- function(x) {
  newdf <- x %>%
    filter(!is.na(referenceCleanedDoi) |
      !is.na(referenceOriginalExternal)) %>%
    filter(!is.na(inchikeySanitized) &
      !is.na(organismLowestTaxon)) %>%
    distinct(
      inchikeySanitized,
      referenceCleanedDoi,
      referenceOriginalExternal,
      organismLowestTaxon,
      .keep_all = TRUE
    ) %>%
    group_by(
      inchikeySanitized,
      referenceCleanedDoi,
      referenceOriginalExternal,
      organism_7_species
    ) %>%
    add_count() %>%
    ungroup() %>%
    filter(!is.na(organismLowestTaxon) |
      !n > 1) %>%
    select(-n) %>%
    group_by(
      inchikeySanitized,
      referenceCleanedDoi,
      referenceOriginalExternal,
      organism_6_genus
    ) %>%
    add_count() %>%
    ungroup() %>%
    filter(!is.na(organismLowestTaxon) |
      !n > 1) %>%
    select(-n) %>%
    group_by(
      inchikeySanitized,
      referenceCleanedDoi,
      referenceOriginalExternal,
      organism_5_family
    ) %>%
    add_count() %>%
    ungroup() %>%
    filter(!is.na(organismLowestTaxon) |
      !n > 1) %>%
    select(-n) %>%
    group_by(
      inchikeySanitized,
      referenceCleanedDoi,
      referenceOriginalExternal,
      organism_4_order
    ) %>%
    add_count() %>%
    ungroup() %>%
    filter(!is.na(organismLowestTaxon) |
      !n > 1) %>%
    select(-n) %>%
    group_by(
      inchikeySanitized,
      referenceCleanedDoi,
      referenceOriginalExternal,
      organism_3_class
    ) %>%
    add_count() %>%
    ungroup() %>%
    filter(!is.na(organismLowestTaxon) |
      !n > 1) %>%
    select(-n) %>%
    group_by(
      inchikeySanitized,
      referenceCleanedDoi,
      referenceOriginalExternal,
      organism_2_phylum
    ) %>%
    add_count() %>%
    ungroup() %>%
    filter(!is.na(organismLowestTaxon) |
      !n > 1) %>%
    select(-n) %>%
    group_by(
      inchikeySanitized,
      referenceCleanedDoi,
      referenceOriginalExternal,
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

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
tcm_standardizing <- function(x) {
  data_bulbus <- x %>%
    filter(grepl("^Bulbus", biologicalsource))
  data_bulbus$newbiologicalsource <-
    gsub("Bulbus", "\\1", data_bulbus$biologicalsource)
  data_bulbus$newbiologicalsource <-
    trimws(data_bulbus$newbiologicalsource)
  data_bulbus$newbiologicalsource <-
    capitalize(data_bulbus$newbiologicalsource)
  if (nrow(data_bulbus) != 0) {
    data_bulbus$newbiologicalsource <-
      paste(data_bulbus$newbiologicalsource, "bulbus")
  }

  data_caulis <- x %>%
    filter(grepl("^Caulis", biologicalsource))
  data_caulis$newbiologicalsource <-
    gsub("Caulis", "\\1", data_caulis$biologicalsource)
  data_caulis$newbiologicalsource <-
    trimws(data_caulis$newbiologicalsource)
  data_caulis$newbiologicalsource <-
    capitalize(data_caulis$newbiologicalsource)
  if (nrow(data_caulis) != 0) {
    data_caulis$newbiologicalsource <-
      paste(data_caulis$newbiologicalsource, "caulis")
  }

  data_caluis_et_folium <- x %>%
    filter(grepl("^Caluis_et_folium", biologicalsource))
  data_caluis_et_folium$newbiologicalsource <-
    gsub(
      "Caluis et folium",
      "\\1",
      data_caluis_et_folium$biologicalsource
    )
  data_caluis_et_folium$newbiologicalsource <-
    trimws(data_caluis_et_folium$newbiologicalsource)
  data_caluis_et_folium$newbiologicalsource <-
    capitalize(data_caluis_et_folium$newbiologicalsource)
  if (nrow(data_caluis_et_folium) != 0) {
    data_caluis_et_folium$newbiologicalsource <-
      paste(
        data_caluis_et_folium$newbiologicalsource,
        "caluis et folium"
      )
  }

  data_corolla <- x %>%
    filter(grepl("^Corolla", biologicalsource))
  data_corolla$newbiologicalsource <-
    gsub("Corolla", "\\1", data_corolla$biologicalsource)
  data_corolla$newbiologicalsource <-
    trimws(data_corolla$newbiologicalsource)
  data_corolla$newbiologicalsource <-
    capitalize(data_corolla$newbiologicalsource)
  if (nrow(data_corolla) != 0) {
    data_corolla$newbiologicalsource <-
      paste(data_corolla$newbiologicalsource, "corolla")
  }

  data_cortex <- x %>%
    filter(grepl("^Cortex", biologicalsource))
  data_cortex$newbiologicalsource <-
    gsub("Cortex", "\\1", data_cortex$biologicalsource)
  data_cortex$newbiologicalsource <-
    trimws(data_cortex$newbiologicalsource)
  data_cortex$newbiologicalsource <-
    capitalize(data_cortex$newbiologicalsource)
  if (nrow(data_cortex) != 0) {
    data_cortex$newbiologicalsource <-
      paste(data_cortex$newbiologicalsource, "cortex")
  }

  data_exocarpium <- x %>%
    filter(grepl("^Exocarpium", biologicalsource))
  data_exocarpium$newbiologicalsource <-
    gsub("Exocarpium", "\\1", data_exocarpium$biologicalsource)
  data_exocarpium$newbiologicalsource <-
    trimws(data_exocarpium$newbiologicalsource)
  data_exocarpium$newbiologicalsource <-
    capitalize(data_exocarpium$newbiologicalsource)
  if (nrow(data_exocarpium) != 0) {
    data_exocarpium$newbiologicalsource <-
      paste(data_exocarpium$newbiologicalsource, "exocarpium")
  }

  data_exocarpium_rubrum <- x %>%
    filter(grepl("^Exocarpium rubrum", biologicalsource))
  data_exocarpium_rubrum$newbiologicalsource <-
    gsub(
      "Exocarpium rubrum",
      "\\1",
      data_exocarpium_rubrum$biologicalsource
    )
  data_exocarpium_rubrum$newbiologicalsource <-
    trimws(data_exocarpium_rubrum$newbiologicalsource)
  data_exocarpium_rubrum$newbiologicalsource <-
    capitalize(data_exocarpium_rubrum$newbiologicalsource)
  if (nrow(data_exocarpium_rubrum) != 0) {
    data_exocarpium_rubrum$newbiologicalsource <-
      paste(
        data_exocarpium_rubrum$newbiologicalsource,
        "exocarpium rubrum"
      )
  }

  data_flos <- x %>%
    filter(grepl("^Flos", biologicalsource))
  data_flos$newbiologicalsource <-
    gsub("Flos", "\\1", data_flos$biologicalsource)
  data_flos$newbiologicalsource <-
    trimws(data_flos$newbiologicalsource)
  data_flos$newbiologicalsource <-
    capitalize(data_flos$newbiologicalsource)
  if (nrow(data_flos) != 0) {
    data_flos$newbiologicalsource <-
      paste(data_flos$newbiologicalsource, "flos")
  }

  data_folium <- x %>%
    filter(grepl("^Folium", biologicalsource))
  data_folium$newbiologicalsource <-
    gsub("Folium", "\\1", data_folium$biologicalsource)
  data_folium$newbiologicalsource <-
    trimws(data_folium$newbiologicalsource)
  data_folium$newbiologicalsource <-
    capitalize(data_folium$newbiologicalsource)
  if (nrow(data_folium) != 0) {
    data_folium$newbiologicalsource <-
      paste(data_folium$newbiologicalsource, "folium")
  }

  data_folium_et_cacumen <- x %>%
    filter(grepl("^Folium et cacumen", biologicalsource))
  data_folium_et_cacumen$newbiologicalsource <-
    gsub(
      "Folium et cacumen",
      "\\1",
      data_folium_et_cacumen$biologicalsource
    )
  data_folium_et_cacumen$newbiologicalsource <-
    trimws(data_folium_et_cacumen$newbiologicalsource)
  data_folium_et_cacumen$newbiologicalsource <-
    capitalize(data_folium_et_cacumen$newbiologicalsource)
  if (nrow(data_folium_et_cacumen) != 0) {
    data_folium_et_cacumen$newbiologicalsource <-
      paste(
        data_folium_et_cacumen$newbiologicalsource,
        "folium et cacumen"
      )
  }

  data_folium_et_caulis <- x %>%
    filter(grepl("^Folium et caulis", biologicalsource))
  data_folium_et_caulis$newbiologicalsource <-
    gsub(
      "Folium et caulis",
      "\\1",
      data_folium_et_caulis$biologicalsource
    )
  data_folium_et_caulis$newbiologicalsource <-
    trimws(data_folium_et_caulis$newbiologicalsource)
  data_folium_et_caulis$newbiologicalsource <-
    capitalize(data_folium_et_caulis$newbiologicalsource)
  if (nrow(data_folium_et_caulis) != 0) {
    data_folium_et_caulis$newbiologicalsource <-
      paste(
        data_folium_et_caulis$newbiologicalsource,
        "folium et caulis"
      )
  }

  data_fructus <- x %>%
    filter(grepl("^Fructus", biologicalsource))
  data_fructus$newbiologicalsource <-
    gsub("Fructus", "\\1", data_fructus$biologicalsource)
  data_fructus$newbiologicalsource <-
    trimws(data_fructus$newbiologicalsource)
  data_fructus$newbiologicalsource <-
    capitalize(data_fructus$newbiologicalsource)
  if (nrow(data_fructus) != 0) {
    data_fructus$newbiologicalsource <-
      paste(data_fructus$newbiologicalsource, "fructus")
  }

  data_fructus_germinatus <- x %>%
    filter(grepl("^Fructus germinatus", biologicalsource))
  data_fructus_germinatus$newbiologicalsource <-
    gsub(
      "Fructus germinatus",
      "\\1",
      data_fructus_germinatus$biologicalsource
    )
  data_fructus_germinatus$newbiologicalsource <-
    trimws(data_fructus_germinatus$newbiologicalsource)
  data_fructus_germinatus$newbiologicalsource <-
    capitalize(data_fructus_germinatus$newbiologicalsource)
  if (nrow(data_fructus_germinatus) != 0) {
    data_fructus_germinatus$newbiologicalsource <-
      paste(
        data_fructus_germinatus$newbiologicalsource,
        "fructus germinatus"
      )
  }

  data_fructus_immaturus <- x %>%
    filter(grepl("^Fructus immaturus", biologicalsource))
  data_fructus_immaturus$newbiologicalsource <-
    gsub(
      "Fructus immaturus",
      "\\1",
      data_fructus_immaturus$biologicalsource
    )
  data_fructus_immaturus$newbiologicalsource <-
    trimws(data_fructus_immaturus$newbiologicalsource)
  data_fructus_immaturus$newbiologicalsource <-
    capitalize(data_fructus_immaturus$newbiologicalsource)
  if (nrow(data_fructus_immaturus) != 0) {
    data_fructus_immaturus$newbiologicalsource <-
      paste(
        data_fructus_immaturus$newbiologicalsource,
        "fructus immaturus"
      )
  }

  data_fructus_retinervus <- x %>%
    filter(grepl("^Fructus retinervus", biologicalsource))
  data_fructus_retinervus$newbiologicalsource <-
    gsub(
      "Fructus retinervus",
      "\\1",
      data_fructus_retinervus$biologicalsource
    )
  data_fructus_retinervus$newbiologicalsource <-
    trimws(data_fructus_retinervus$newbiologicalsource)
  data_fructus_retinervus$newbiologicalsource <-
    capitalize(data_fructus_retinervus$newbiologicalsource)
  if (nrow(data_fructus_retinervus) != 0) {
    data_fructus_retinervus$newbiologicalsource <-
      paste(
        data_fructus_retinervus$newbiologicalsource,
        "fructus retinervus"
      )
  }

  data_fructus_rotundus <- x %>%
    filter(grepl("^Fructus rotundus", biologicalsource))
  data_fructus_rotundus$newbiologicalsource <-
    gsub(
      "Fructus rotundus",
      "\\1",
      data_fructus_rotundus$biologicalsource
    )
  data_fructus_rotundus$newbiologicalsource <-
    trimws(data_fructus_rotundus$newbiologicalsource)
  data_fructus_rotundus$newbiologicalsource <-
    capitalize(data_fructus_rotundus$newbiologicalsource)
  if (nrow(data_fructus_rotundus) != 0) {
    data_fructus_rotundus$newbiologicalsource <-
      paste(
        data_fructus_rotundus$newbiologicalsource,
        "fructus rotundus"
      )
  }

  data_herba <- x %>%
    filter(grepl("^Herba", biologicalsource))
  data_herba$newbiologicalsource <-
    gsub("Herba", "\\1", data_herba$biologicalsource)
  data_herba$newbiologicalsource <-
    trimws(data_herba$newbiologicalsource)
  data_herba$newbiologicalsource <-
    capitalize(data_herba$newbiologicalsource)
  if (nrow(data_herba) != 0) {
    data_herba$newbiologicalsource <-
      paste(data_herba$newbiologicalsource, "herba")
  }

  data_lignum <- x %>%
    filter(grepl("^Lignum", biologicalsource))
  data_lignum$newbiologicalsource <-
    gsub("Lignum", "\\1", data_lignum$biologicalsource)
  data_lignum$newbiologicalsource <-
    trimws(data_lignum$newbiologicalsource)
  data_lignum$newbiologicalsource <-
    capitalize(data_lignum$newbiologicalsource)
  if (nrow(data_lignum) != 0) {
    data_lignum$newbiologicalsource <-
      paste(data_lignum$newbiologicalsource, "lignum")
  }

  data_medulla <- x %>%
    filter(grepl("^Medulla", biologicalsource))
  data_medulla$newbiologicalsource <-
    gsub("Medulla", "\\1", data_medulla$biologicalsource)
  data_medulla$newbiologicalsource <-
    trimws(data_medulla$newbiologicalsource)
  data_medulla$newbiologicalsource <-
    capitalize(data_medulla$newbiologicalsource)
  if (nrow(data_medulla) != 0) {
    data_medulla$newbiologicalsource <-
      paste(data_medulla$newbiologicalsource, "medulla")
  }

  data_pericarpum <- x %>%
    filter(grepl("^Pericarpum", biologicalsource))
  data_pericarpum$newbiologicalsource <-
    gsub("Pericarpum", "\\1", data_pericarpum$biologicalsource)
  data_pericarpum$newbiologicalsource <-
    trimws(data_pericarpum$newbiologicalsource)
  data_pericarpum$newbiologicalsource <-
    capitalize(data_pericarpum$newbiologicalsource)
  if (nrow(data_pericarpum) != 0) {
    data_pericarpum$newbiologicalsource <-
      paste(data_pericarpum$newbiologicalsource, "pericarpum")
  }

  data_petiolus <- x %>%
    filter(grepl("^Petiolus", biologicalsource))
  data_petiolus$newbiologicalsource <-
    gsub("Petiolus", "\\1", data_petiolus$biologicalsource)
  data_petiolus$newbiologicalsource <-
    trimws(data_petiolus$newbiologicalsource)
  data_petiolus$newbiologicalsource <-
    capitalize(data_petiolus$newbiologicalsource)
  if (nrow(data_petiolus) != 0) {
    data_petiolus$newbiologicalsource <-
      paste(data_petiolus$newbiologicalsource, "petiolus")
  }

  data_pollen <- x %>%
    filter(grepl("^Pollen", biologicalsource))
  data_pollen$newbiologicalsource <-
    gsub("Pollen", "\\1", data_pollen$biologicalsource)
  data_pollen$newbiologicalsource <-
    trimws(data_pollen$newbiologicalsource)
  data_pollen$newbiologicalsource <-
    capitalize(data_pollen$newbiologicalsource)
  if (nrow(data_pollen) != 0) {
    data_pollen$newbiologicalsource <-
      paste(data_pollen$newbiologicalsource, "pollen")
  }

  data_radicis_cortex <- x %>%
    filter(grepl("^Radicis cortex", biologicalsource))
  data_radicis_cortex$newbiologicalsource <-
    gsub(
      "Radicis cortex",
      "\\1",
      data_radicis_cortex$biologicalsource
    )
  data_radicis_cortex$newbiologicalsource <-
    trimws(data_radicis_cortex$newbiologicalsource)
  data_radicis_cortex$newbiologicalsource <-
    capitalize(data_radicis_cortex$newbiologicalsource)
  if (nrow(data_radicis_cortex) != 0) {
    data_radicis_cortex$newbiologicalsource <-
      paste(data_radicis_cortex$newbiologicalsource, "radicis cortex")
  }

  data_radix <- x %>%
    filter(grepl("^Radix", biologicalsource))
  data_radix$newbiologicalsource <-
    gsub("Radix", "\\1", data_radix$biologicalsource)
  data_radix$newbiologicalsource <-
    trimws(data_radix$newbiologicalsource)
  data_radix$newbiologicalsource <-
    capitalize(data_radix$newbiologicalsource)
  if (nrow(data_radix) != 0) {
    data_radix$newbiologicalsource <-
      paste(data_radix$newbiologicalsource, "radix")
  }

  data_radix_et_rhizoma <- x %>%
    filter(grepl("^Radix et rhizoma", biologicalsource))
  data_radix_et_rhizoma$newbiologicalsource <-
    gsub(
      "Radix et rhizoma",
      "\\1",
      data_radix_et_rhizoma$biologicalsource
    )
  data_radix_et_rhizoma$newbiologicalsource <-
    trimws(data_radix_et_rhizoma$newbiologicalsource)
  data_radix_et_rhizoma$newbiologicalsource <-
    capitalize(data_radix_et_rhizoma$newbiologicalsource)
  if (nrow(data_radix_et_rhizoma) != 0) {
    data_radix_et_rhizoma$newbiologicalsource <-
      paste(
        data_radix_et_rhizoma$newbiologicalsource,
        "radix et rhizoma"
      )
  }

  data_radix_preparata <- x %>%
    filter(grepl("^Radix preparata", biologicalsource))
  data_radix_preparata$newbiologicalsource <-
    gsub(
      "Radix preparata",
      "\\1",
      data_radix_preparata$biologicalsource
    )
  data_radix_preparata$newbiologicalsource <-
    trimws(data_radix_preparata$newbiologicalsource)
  data_radix_preparata$newbiologicalsource <-
    capitalize(data_radix_preparata$newbiologicalsource)
  if (nrow(data_radix_preparata) != 0) {
    data_radix_preparata$newbiologicalsource <-
      paste(
        data_radix_preparata$newbiologicalsource,
        "radix preparata"
      )
  }

  data_ramulus <- x %>%
    filter(grepl("^Ramulus", biologicalsource))
  data_ramulus$newbiologicalsource <-
    gsub("Ramulus", "\\1", data_ramulus$biologicalsource)
  data_ramulus$newbiologicalsource <-
    trimws(data_ramulus$newbiologicalsource)
  data_ramulus$newbiologicalsource <-
    capitalize(data_ramulus$newbiologicalsource)
  if (nrow(data_ramulus) != 0) {
    data_ramulus$newbiologicalsource <-
      paste(data_ramulus$newbiologicalsource, "ramulus")
  }

  data_ramulus_cum_uncus <- x %>%
    filter(grepl("^Ramulus cum uncus", biologicalsource))
  data_ramulus_cum_uncus$newbiologicalsource <-
    gsub(
      "Ramulus cum uncus",
      "\\1",
      data_ramulus_cum_uncus$biologicalsource
    )
  data_ramulus_cum_uncus$newbiologicalsource <-
    trimws(data_ramulus_cum_uncus$newbiologicalsource)
  data_ramulus_cum_uncus$newbiologicalsource <-
    capitalize(data_ramulus_cum_uncus$newbiologicalsource)
  if (nrow(data_ramulus_cum_uncus) != 0) {
    data_ramulus_cum_uncus$newbiologicalsource <-
      paste(
        data_ramulus_cum_uncus$newbiologicalsource,
        "ramulus cum uncus"
      )
  }

  data_ramulus_et_folium <- x %>%
    filter(grepl("^Ramulus et folium", biologicalsource))
  data_ramulus_et_folium$newbiologicalsource <-
    gsub(
      "Ramulus et folium",
      "\\1",
      data_ramulus_et_folium$biologicalsource
    )
  data_ramulus_et_folium$newbiologicalsource <-
    trimws(data_ramulus_et_folium$newbiologicalsource)
  data_ramulus_et_folium$newbiologicalsource <-
    capitalize(data_ramulus_et_folium$newbiologicalsource)
  if (nrow(data_ramulus_et_folium) != 0) {
    data_ramulus_et_folium$newbiologicalsource <-
      paste(
        data_ramulus_et_folium$newbiologicalsource,
        "ramulus et folium"
      )
  }

  data_rhizoma <- x %>%
    filter(grepl("^Rhizoma", biologicalsource))
  data_rhizoma$newbiologicalsource <-
    gsub("Rhizoma", "\\1", data_rhizoma$biologicalsource)
  data_rhizoma$newbiologicalsource <-
    trimws(data_rhizoma$newbiologicalsource)
  data_rhizoma$newbiologicalsource <-
    capitalize(data_rhizoma$newbiologicalsource)
  if (nrow(data_rhizoma) != 0) {
    data_rhizoma$newbiologicalsource <-
      paste(data_rhizoma$newbiologicalsource, "rhizoma")
  }

  data_rhizoma_alba <- x %>%
    filter(grepl("^Rhizoma alba", biologicalsource))
  data_rhizoma_alba$newbiologicalsource <-
    gsub("Rhizoma alba", "\\1", data_rhizoma_alba$biologicalsource)
  data_rhizoma_alba$newbiologicalsource <-
    trimws(data_rhizoma_alba$newbiologicalsource)
  data_rhizoma_alba$newbiologicalsource <-
    capitalize(data_rhizoma_alba$newbiologicalsource)
  if (nrow(data_rhizoma_alba) != 0) {
    data_rhizoma_alba$newbiologicalsource <-
      paste(data_rhizoma_alba$newbiologicalsource, "rhizoma alba")
  }

  data_rhizoma_et_radix <- x %>%
    filter(grepl("^Rhizoma et radix", biologicalsource))
  data_rhizoma_et_radix$newbiologicalsource <-
    gsub(
      "Rhizoma et radix",
      "\\1",
      data_rhizoma_et_radix$biologicalsource
    )
  data_rhizoma_et_radix$newbiologicalsource <-
    trimws(data_rhizoma_et_radix$newbiologicalsource)
  data_rhizoma_et_radix$newbiologicalsource <-
    capitalize(data_rhizoma_et_radix$newbiologicalsource)
  if (nrow(data_rhizoma_et_radix) != 0) {
    data_rhizoma_et_radix$newbiologicalsource <-
      paste(
        data_rhizoma_et_radix$newbiologicalsource,
        "rhizoma et radix"
      )
  }

  data_semen <- x %>%
    filter(grepl("^Semen", biologicalsource))
  data_semen$newbiologicalsource <-
    gsub("Semen", "\\1", data_semen$biologicalsource)
  data_semen$newbiologicalsource <-
    trimws(data_semen$newbiologicalsource)
  data_semen$newbiologicalsource <-
    capitalize(data_semen$newbiologicalsource)
  if (nrow(data_semen) != 0) {
    data_semen$newbiologicalsource <-
      paste(data_semen$newbiologicalsource, "semen")
  }

  data_semen_germinatum <- x %>%
    filter(grepl("^Semen germinatum", biologicalsource))
  data_semen_germinatum$newbiologicalsource <-
    gsub(
      "Semen germinatum",
      "\\1",
      data_semen_germinatum$biologicalsource
    )
  data_semen_germinatum$newbiologicalsource <-
    trimws(data_semen_germinatum$newbiologicalsource)
  data_semen_germinatum$newbiologicalsource <-
    capitalize(data_semen_germinatum$newbiologicalsource)
  if (nrow(data_semen_germinatum) != 0) {
    data_semen_germinatum$newbiologicalsource <-
      paste(
        data_semen_germinatum$newbiologicalsource,
        "semen germinatum"
      )
  }

  data_spica <- x %>%
    filter(grepl("Spica ", biologicalsource, fixed = TRUE))
  data_spica$newbiologicalsource <-
    gsub("Spica", "\\1", data_spica$biologicalsource)
  data_spica$newbiologicalsource <-
    trimws(data_spica$newbiologicalsource)
  data_spica$newbiologicalsource <-
    capitalize(data_spica$newbiologicalsource)
  if (nrow(data_spica) != 0) {
    data_spica$newbiologicalsource <-
      paste(data_spica$newbiologicalsource, "spica")
  }

  data_stamen <- x %>%
    filter(grepl("^Stamen", biologicalsource))
  data_stamen$newbiologicalsource <-
    gsub("Stamen", "\\1", data_stamen$biologicalsource)
  data_stamen$newbiologicalsource <-
    trimws(data_stamen$newbiologicalsource)
  data_stamen$newbiologicalsource <-
    capitalize(data_stamen$newbiologicalsource)
  if (nrow(data_stamen) != 0) {
    data_stamen$newbiologicalsource <-
      paste(data_stamen$newbiologicalsource, "stamen")
  }

  data_stigma <- x %>%
    filter(grepl("Stigma ", biologicalsource, fixed = TRUE))
  data_stigma$newbiologicalsource <-
    gsub("Stigma", "\\1", data_stigma$biologicalsource)
  data_stigma$newbiologicalsource <-
    trimws(data_stigma$newbiologicalsource)
  data_stigma$newbiologicalsource <-
    capitalize(data_stigma$newbiologicalsource)
  if (nrow(data_stigma) != 0) {
    data_stigma$newbiologicalsource <-
      paste(data_stigma$newbiologicalsource, "stigma")
  }

  data_storax <- x %>%
    filter(grepl("^Storax", biologicalsource))
  data_storax$newbiologicalsource <-
    gsub("Storax", "\\1", data_storax$biologicalsource)
  data_storax$newbiologicalsource <-
    trimws(data_storax$newbiologicalsource)
  data_storax$newbiologicalsource <-
    capitalize(data_storax$newbiologicalsource)
  if (nrow(data_storax) != 0) {
    data_storax$newbiologicalsource <-
      paste(data_storax$newbiologicalsource, "storax")
  }

  data_thallus <- x %>%
    filter(grepl("^Thallus", biologicalsource, fixed = TRUE))
  data_thallus$newbiologicalsource <-
    gsub("Thallus", "\\1", data_thallus$biologicalsource)
  data_thallus$newbiologicalsource <-
    trimws(data_thallus$newbiologicalsource)
  data_thallus$newbiologicalsource <-
    capitalize(data_thallus$newbiologicalsource)
  if (nrow(data_thallus) != 0) {
    data_thallus$newbiologicalsource <-
      paste(data_thallus$newbiologicalsource, "thallus")
  }
  # not tuber

  x_2 <-
    rbind(
      data_bulbus,
      data_caulis,
      data_caluis_et_folium,
      data_corolla,
      data_cortex,
      data_exocarpium,
      data_exocarpium_rubrum,
      data_flos,
      data_folium,
      data_folium_et_cacumen,
      data_folium_et_caulis,
      data_fructus,
      data_fructus_germinatus,
      data_fructus_immaturus,
      data_fructus_retinervus,
      data_fructus_rotundus,
      data_herba,
      data_lignum,
      data_medulla,
      data_pericarpum,
      data_petiolus,
      data_pollen,
      data_radicis_cortex,
      data_radix,
      data_rhizoma_et_radix,
      data_radix_et_rhizoma,
      data_radix_preparata,
      data_ramulus,
      data_ramulus_cum_uncus,
      data_ramulus_et_folium,
      data_rhizoma,
      data_rhizoma_alba,
      data_rhizoma_et_radix,
      data_semen,
      data_semen_germinatum,
      data_spica,
      data_stamen,
      data_stigma,
      data_storax,
      data_thallus
    )

  x_3 <- left_join(x, x_2)

  x_3$newnewbiologicalsource <-
    apply(x_3, 1, function(x) {
      tail(na.omit(x), 1)
    })

  x_3 <- x_3 %>%
    select(
      biologicalsource,
      newnewbiologicalsource
    )

  x_4 <-
    left_join(x_3, food_names_list, by = c("newnewbiologicalsource" = "name"))

  x_4 <-
    left_join(x_4, tcm_names_list, by = c("newnewbiologicalsource" = "latin"))

  x_4 <-
    left_join(x_4,
      tcm_names_list,
      by = c("newnewbiologicalsource" = "common")
    )

  x_4$newnewnewbiologicalsource <-
    apply(x_4, 1, function(x) {
      tail(na.omit(x), 1)
    })

  x_4 <- x_4 %>%
    select(
      biologicalsource = biologicalsource.x,
      newnewnewbiologicalsource
    )

  data_standard_2 <- left_join(data_standard, x_4) %>%
    select(-biologicalsource) %>%
    select(
      biologicalsource = newnewnewbiologicalsource,
      everything()
    )

  data_standard_2
}

#######################################################

#' Title
#'
#' @param i
#'
#' @return
#' @export
#'
#' @examples
name2inchi <- function(i)
{
  tryCatch(
    {
      cpd <- dataTranslatedNominal[i, "nameSanitized"]
      url <-
        paste("https://cactus.nci.nih.gov/chemical/structure/",
          cpd,
          "/stdinchi",
          sep = ""
        )
      url <- gsub(
        pattern = "\\s",
        replacement = "%20",
        x = url
      )
      read_html(url) %>%
        html_text()
    },
    error = function(e) {
      NA
    }
  )
}
# {
#   tryCatch({
#     x <- cts_convert(
#       query = dataTranslatedNominal[i, "nameSanitized"],
#       from = "Chemical Name",
#       to = "InChI Code",
#       verbose = FALSE,
#       choices = 1
#     )
#     return(x)
#   }
#   , error = function(e) {
#     NA
#   })
# }

#######################################################

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
tcm_inverting <- function(x) {
  data_bulbus <- x %>%
    filter(grepl("^Bulbus", biologicalsource))
  data_bulbus$newbiologicalsource <-
    gsub("Bulbus", "\\1", data_bulbus$biologicalsource)
  data_bulbus$newbiologicalsource <-
    trimws(data_bulbus$newbiologicalsource)
  data_bulbus$newbiologicalsource <-
    capitalize(data_bulbus$newbiologicalsource)
  if (nrow(data_bulbus) != 0) {
    data_bulbus$newbiologicalsource <-
      paste(data_bulbus$newbiologicalsource, "bulbus")
  }

  data_caulis <- x %>%
    filter(grepl("^Caulis", biologicalsource))
  data_caulis$newbiologicalsource <-
    gsub("Caulis", "\\1", data_caulis$biologicalsource)
  data_caulis$newbiologicalsource <-
    trimws(data_caulis$newbiologicalsource)
  data_caulis$newbiologicalsource <-
    capitalize(data_caulis$newbiologicalsource)
  if (nrow(data_caulis) != 0) {
    data_caulis$newbiologicalsource <-
      paste(data_caulis$newbiologicalsource, "caulis")
  }

  data_caluis_et_folium <- x %>%
    filter(grepl("^Caluis_et_folium", biologicalsource))
  data_caluis_et_folium$newbiologicalsource <-
    gsub(
      "Caluis et folium",
      "\\1",
      data_caluis_et_folium$biologicalsource
    )
  data_caluis_et_folium$newbiologicalsource <-
    trimws(data_caluis_et_folium$newbiologicalsource)
  data_caluis_et_folium$newbiologicalsource <-
    capitalize(data_caluis_et_folium$newbiologicalsource)
  if (nrow(data_caluis_et_folium) != 0) {
    data_caluis_et_folium$newbiologicalsource <-
      paste(
        data_caluis_et_folium$newbiologicalsource,
        "caluis et folium"
      )
  }

  data_corolla <- x %>%
    filter(grepl("^Corolla", biologicalsource))
  data_corolla$newbiologicalsource <-
    gsub("Corolla", "\\1", data_corolla$biologicalsource)
  data_corolla$newbiologicalsource <-
    trimws(data_corolla$newbiologicalsource)
  data_corolla$newbiologicalsource <-
    capitalize(data_corolla$newbiologicalsource)
  if (nrow(data_corolla) != 0) {
    data_corolla$newbiologicalsource <-
      paste(data_corolla$newbiologicalsource, "corolla")
  }

  data_cortex <- x %>%
    filter(grepl("^Cortex", biologicalsource))
  data_cortex$newbiologicalsource <-
    gsub("Cortex", "\\1", data_cortex$biologicalsource)
  data_cortex$newbiologicalsource <-
    trimws(data_cortex$newbiologicalsource)
  data_cortex$newbiologicalsource <-
    capitalize(data_cortex$newbiologicalsource)
  if (nrow(data_cortex) != 0) {
    data_cortex$newbiologicalsource <-
      paste(data_cortex$newbiologicalsource, "cortex")
  }

  data_exocarpium <- x %>%
    filter(grepl("^Exocarpium", biologicalsource))
  data_exocarpium$newbiologicalsource <-
    gsub("Exocarpium", "\\1", data_exocarpium$biologicalsource)
  data_exocarpium$newbiologicalsource <-
    trimws(data_exocarpium$newbiologicalsource)
  data_exocarpium$newbiologicalsource <-
    capitalize(data_exocarpium$newbiologicalsource)
  if (nrow(data_exocarpium) != 0) {
    data_exocarpium$newbiologicalsource <-
      paste(data_exocarpium$newbiologicalsource, "exocarpium")
  }

  data_exocarpium_rubrum <- x %>%
    filter(grepl("^Exocarpium rubrum", biologicalsource))
  data_exocarpium_rubrum$newbiologicalsource <-
    gsub(
      "Exocarpium rubrum",
      "\\1",
      data_exocarpium_rubrum$biologicalsource
    )
  data_exocarpium_rubrum$newbiologicalsource <-
    trimws(data_exocarpium_rubrum$newbiologicalsource)
  data_exocarpium_rubrum$newbiologicalsource <-
    capitalize(data_exocarpium_rubrum$newbiologicalsource)
  if (nrow(data_exocarpium_rubrum) != 0) {
    data_exocarpium_rubrum$newbiologicalsource <-
      paste(
        data_exocarpium_rubrum$newbiologicalsource,
        "exocarpium rubrum"
      )
  }

  data_flos <- x %>%
    filter(grepl("^Flos", biologicalsource))
  data_flos$newbiologicalsource <-
    gsub("Flos", "\\1", data_flos$biologicalsource)
  data_flos$newbiologicalsource <-
    trimws(data_flos$newbiologicalsource)
  data_flos$newbiologicalsource <-
    capitalize(data_flos$newbiologicalsource)
  if (nrow(data_flos) != 0) {
    data_flos$newbiologicalsource <-
      paste(data_flos$newbiologicalsource, "flos")
  }

  data_folium <- x %>%
    filter(grepl("^Folium", biologicalsource))
  data_folium$newbiologicalsource <-
    gsub("Folium", "\\1", data_folium$biologicalsource)
  data_folium$newbiologicalsource <-
    trimws(data_folium$newbiologicalsource)
  data_folium$newbiologicalsource <-
    capitalize(data_folium$newbiologicalsource)
  if (nrow(data_folium) != 0) {
    data_folium$newbiologicalsource <-
      paste(data_folium$newbiologicalsource, "folium")
  }

  data_folium_et_cacumen <- x %>%
    filter(grepl("^Folium et cacumen", biologicalsource))
  data_folium_et_cacumen$newbiologicalsource <-
    gsub(
      "Folium et cacumen",
      "\\1",
      data_folium_et_cacumen$biologicalsource
    )
  data_folium_et_cacumen$newbiologicalsource <-
    trimws(data_folium_et_cacumen$newbiologicalsource)
  data_folium_et_cacumen$newbiologicalsource <-
    capitalize(data_folium_et_cacumen$newbiologicalsource)
  if (nrow(data_folium_et_cacumen) != 0) {
    data_folium_et_cacumen$newbiologicalsource <-
      paste(
        data_folium_et_cacumen$newbiologicalsource,
        "folium et cacumen"
      )
  }

  data_folium_et_caulis <- x %>%
    filter(grepl("^Folium et caulis", biologicalsource))
  data_folium_et_caulis$newbiologicalsource <-
    gsub(
      "Folium et caulis",
      "\\1",
      data_folium_et_caulis$biologicalsource
    )
  data_folium_et_caulis$newbiologicalsource <-
    trimws(data_folium_et_caulis$newbiologicalsource)
  data_folium_et_caulis$newbiologicalsource <-
    capitalize(data_folium_et_caulis$newbiologicalsource)
  if (nrow(data_folium_et_caulis) != 0) {
    data_folium_et_caulis$newbiologicalsource <-
      paste(
        data_folium_et_caulis$newbiologicalsource,
        "folium et caulis"
      )
  }

  data_fructus <- x %>%
    filter(grepl("^Fructus", biologicalsource))
  data_fructus$newbiologicalsource <-
    gsub("Fructus", "\\1", data_fructus$biologicalsource)
  data_fructus$newbiologicalsource <-
    trimws(data_fructus$newbiologicalsource)
  data_fructus$newbiologicalsource <-
    capitalize(data_fructus$newbiologicalsource)
  if (nrow(data_fructus) != 0) {
    data_fructus$newbiologicalsource <-
      paste(data_fructus$newbiologicalsource, "fructus")
  }

  data_fructus_germinatus <- x %>%
    filter(grepl("^Fructus germinatus", biologicalsource))
  data_fructus_germinatus$newbiologicalsource <-
    gsub(
      "Fructus germinatus",
      "\\1",
      data_fructus_germinatus$biologicalsource
    )
  data_fructus_germinatus$newbiologicalsource <-
    trimws(data_fructus_germinatus$newbiologicalsource)
  data_fructus_germinatus$newbiologicalsource <-
    capitalize(data_fructus_germinatus$newbiologicalsource)
  if (nrow(data_fructus_germinatus) != 0) {
    data_fructus_germinatus$newbiologicalsource <-
      paste(
        data_fructus_germinatus$newbiologicalsource,
        "fructus germinatus"
      )
  }

  data_fructus_immaturus <- x %>%
    filter(grepl("^Fructus immaturus", biologicalsource))
  data_fructus_immaturus$newbiologicalsource <-
    gsub(
      "Fructus immaturus",
      "\\1",
      data_fructus_immaturus$biologicalsource
    )
  data_fructus_immaturus$newbiologicalsource <-
    trimws(data_fructus_immaturus$newbiologicalsource)
  data_fructus_immaturus$newbiologicalsource <-
    capitalize(data_fructus_immaturus$newbiologicalsource)
  if (nrow(data_fructus_immaturus) != 0) {
    data_fructus_immaturus$newbiologicalsource <-
      paste(
        data_fructus_immaturus$newbiologicalsource,
        "fructus immaturus"
      )
  }

  data_fructus_retinervus <- x %>%
    filter(grepl("^Fructus retinervus", biologicalsource))
  data_fructus_retinervus$newbiologicalsource <-
    gsub(
      "Fructus retinervus",
      "\\1",
      data_fructus_retinervus$biologicalsource
    )
  data_fructus_retinervus$newbiologicalsource <-
    trimws(data_fructus_retinervus$newbiologicalsource)
  data_fructus_retinervus$newbiologicalsource <-
    capitalize(data_fructus_retinervus$newbiologicalsource)
  if (nrow(data_fructus_retinervus) != 0) {
    data_fructus_retinervus$newbiologicalsource <-
      paste(
        data_fructus_retinervus$newbiologicalsource,
        "fructus retinervus"
      )
  }

  data_fructus_rotundus <- x %>%
    filter(grepl("^Fructus rotundus", biologicalsource))
  data_fructus_rotundus$newbiologicalsource <-
    gsub(
      "Fructus rotundus",
      "\\1",
      data_fructus_rotundus$biologicalsource
    )
  data_fructus_rotundus$newbiologicalsource <-
    trimws(data_fructus_rotundus$newbiologicalsource)
  data_fructus_rotundus$newbiologicalsource <-
    capitalize(data_fructus_rotundus$newbiologicalsource)
  if (nrow(data_fructus_rotundus) != 0) {
    data_fructus_rotundus$newbiologicalsource <-
      paste(
        data_fructus_rotundus$newbiologicalsource,
        "fructus rotundus"
      )
  }

  data_herba <- x %>%
    filter(grepl("^Herba", biologicalsource))
  data_herba$newbiologicalsource <-
    gsub("Herba", "\\1", data_herba$biologicalsource)
  data_herba$newbiologicalsource <-
    trimws(data_herba$newbiologicalsource)
  data_herba$newbiologicalsource <-
    capitalize(data_herba$newbiologicalsource)
  if (nrow(data_herba) != 0) {
    data_herba$newbiologicalsource <-
      paste(data_herba$newbiologicalsource, "herba")
  }

  data_lignum <- x %>%
    filter(grepl("^Lignum", biologicalsource))
  data_lignum$newbiologicalsource <-
    gsub("Lignum", "\\1", data_lignum$biologicalsource)
  data_lignum$newbiologicalsource <-
    trimws(data_lignum$newbiologicalsource)
  data_lignum$newbiologicalsource <-
    capitalize(data_lignum$newbiologicalsource)
  if (nrow(data_lignum) != 0) {
    data_lignum$newbiologicalsource <-
      paste(data_lignum$newbiologicalsource, "lignum")
  }

  data_medulla <- x %>%
    filter(grepl("^Medulla", biologicalsource))
  data_medulla$newbiologicalsource <-
    gsub("Medulla", "\\1", data_medulla$biologicalsource)
  data_medulla$newbiologicalsource <-
    trimws(data_medulla$newbiologicalsource)
  data_medulla$newbiologicalsource <-
    capitalize(data_medulla$newbiologicalsource)
  if (nrow(data_medulla) != 0) {
    data_medulla$newbiologicalsource <-
      paste(data_medulla$newbiologicalsource, "medulla")
  }

  data_pericarpum <- x %>%
    filter(grepl("^Pericarpum", biologicalsource))
  data_pericarpum$newbiologicalsource <-
    gsub("Pericarpum", "\\1", data_pericarpum$biologicalsource)
  data_pericarpum$newbiologicalsource <-
    trimws(data_pericarpum$newbiologicalsource)
  data_pericarpum$newbiologicalsource <-
    capitalize(data_pericarpum$newbiologicalsource)
  if (nrow(data_pericarpum) != 0) {
    data_pericarpum$newbiologicalsource <-
      paste(data_pericarpum$newbiologicalsource, "pericarpum")
  }

  data_petiolus <- x %>%
    filter(grepl("^Petiolus", biologicalsource))
  data_petiolus$newbiologicalsource <-
    gsub("Petiolus", "\\1", data_petiolus$biologicalsource)
  data_petiolus$newbiologicalsource <-
    trimws(data_petiolus$newbiologicalsource)
  data_petiolus$newbiologicalsource <-
    capitalize(data_petiolus$newbiologicalsource)
  if (nrow(data_petiolus) != 0) {
    data_petiolus$newbiologicalsource <-
      paste(data_petiolus$newbiologicalsource, "petiolus")
  }

  data_pollen <- x %>%
    filter(grepl("^Pollen", biologicalsource))
  data_pollen$newbiologicalsource <-
    gsub("Pollen", "\\1", data_pollen$biologicalsource)
  data_pollen$newbiologicalsource <-
    trimws(data_pollen$newbiologicalsource)
  data_pollen$newbiologicalsource <-
    capitalize(data_pollen$newbiologicalsource)
  if (nrow(data_pollen) != 0) {
    data_pollen$newbiologicalsource <-
      paste(data_pollen$newbiologicalsource, "pollen")
  }

  data_radicis_cortex <- x %>%
    filter(grepl("^Radicis cortex", biologicalsource))
  data_radicis_cortex$newbiologicalsource <-
    gsub(
      "Radicis cortex",
      "\\1",
      data_radicis_cortex$biologicalsource
    )
  data_radicis_cortex$newbiologicalsource <-
    trimws(data_radicis_cortex$newbiologicalsource)
  data_radicis_cortex$newbiologicalsource <-
    capitalize(data_radicis_cortex$newbiologicalsource)
  if (nrow(data_radicis_cortex) != 0) {
    data_radicis_cortex$newbiologicalsource <-
      paste(data_radicis_cortex$newbiologicalsource, "radicis cortex")
  }

  data_radix <- x %>%
    filter(grepl("^Radix", biologicalsource))
  data_radix$newbiologicalsource <-
    gsub("Radix", "\\1", data_radix$biologicalsource)
  data_radix$newbiologicalsource <-
    trimws(data_radix$newbiologicalsource)
  data_radix$newbiologicalsource <-
    capitalize(data_radix$newbiologicalsource)
  if (nrow(data_radix) != 0) {
    data_radix$newbiologicalsource <-
      paste(data_radix$newbiologicalsource, "radix")
  }

  data_radix_et_rhizoma <- x %>%
    filter(grepl("^Radix et rhizoma", biologicalsource))
  data_radix_et_rhizoma$newbiologicalsource <-
    gsub(
      "Radix et rhizoma",
      "\\1",
      data_radix_et_rhizoma$biologicalsource
    )
  data_radix_et_rhizoma$newbiologicalsource <-
    trimws(data_radix_et_rhizoma$newbiologicalsource)
  data_radix_et_rhizoma$newbiologicalsource <-
    capitalize(data_radix_et_rhizoma$newbiologicalsource)
  if (nrow(data_radix_et_rhizoma) != 0) {
    data_radix_et_rhizoma$newbiologicalsource <-
      paste(
        data_radix_et_rhizoma$newbiologicalsource,
        "radix et rhizoma"
      )
  }

  data_radix_preparata <- x %>%
    filter(grepl("^Radix preparata", biologicalsource))
  data_radix_preparata$newbiologicalsource <-
    gsub(
      "Radix preparata",
      "\\1",
      data_radix_preparata$biologicalsource
    )
  data_radix_preparata$newbiologicalsource <-
    trimws(data_radix_preparata$newbiologicalsource)
  data_radix_preparata$newbiologicalsource <-
    capitalize(data_radix_preparata$newbiologicalsource)
  if (nrow(data_radix_preparata) != 0) {
    data_radix_preparata$newbiologicalsource <-
      paste(
        data_radix_preparata$newbiologicalsource,
        "radix preparata"
      )
  }

  data_ramulus <- x %>%
    filter(grepl("^Ramulus", biologicalsource))
  data_ramulus$newbiologicalsource <-
    gsub("Ramulus", "\\1", data_ramulus$biologicalsource)
  data_ramulus$newbiologicalsource <-
    trimws(data_ramulus$newbiologicalsource)
  data_ramulus$newbiologicalsource <-
    capitalize(data_ramulus$newbiologicalsource)
  if (nrow(data_ramulus) != 0) {
    data_ramulus$newbiologicalsource <-
      paste(data_ramulus$newbiologicalsource, "ramulus")
  }

  data_ramulus_cum_uncus <- x %>%
    filter(grepl("^Ramulus cum uncus", biologicalsource))
  data_ramulus_cum_uncus$newbiologicalsource <-
    gsub(
      "Ramulus cum uncus",
      "\\1",
      data_ramulus_cum_uncus$biologicalsource
    )
  data_ramulus_cum_uncus$newbiologicalsource <-
    trimws(data_ramulus_cum_uncus$newbiologicalsource)
  data_ramulus_cum_uncus$newbiologicalsource <-
    capitalize(data_ramulus_cum_uncus$newbiologicalsource)
  if (nrow(data_ramulus_cum_uncus) != 0) {
    data_ramulus_cum_uncus$newbiologicalsource <-
      paste(
        data_ramulus_cum_uncus$newbiologicalsource,
        "ramulus cum uncus"
      )
  }

  data_ramulus_et_folium <- x %>%
    filter(grepl("^Ramulus et folium", biologicalsource))
  data_ramulus_et_folium$newbiologicalsource <-
    gsub(
      "Ramulus et folium",
      "\\1",
      data_ramulus_et_folium$biologicalsource
    )
  data_ramulus_et_folium$newbiologicalsource <-
    trimws(data_ramulus_et_folium$newbiologicalsource)
  data_ramulus_et_folium$newbiologicalsource <-
    capitalize(data_ramulus_et_folium$newbiologicalsource)
  if (nrow(data_ramulus_et_folium) != 0) {
    data_ramulus_et_folium$newbiologicalsource <-
      paste(
        data_ramulus_et_folium$newbiologicalsource,
        "ramulus et folium"
      )
  }

  data_rhizoma <- x %>%
    filter(grepl("^Rhizoma", biologicalsource))
  data_rhizoma$newbiologicalsource <-
    gsub("Rhizoma", "\\1", data_rhizoma$biologicalsource)
  data_rhizoma$newbiologicalsource <-
    trimws(data_rhizoma$newbiologicalsource)
  data_rhizoma$newbiologicalsource <-
    capitalize(data_rhizoma$newbiologicalsource)
  if (nrow(data_rhizoma) != 0) {
    data_rhizoma$newbiologicalsource <-
      paste(data_rhizoma$newbiologicalsource, "rhizoma")
  }

  data_rhizoma_alba <- x %>%
    filter(grepl("^Rhizoma alba", biologicalsource))
  data_rhizoma_alba$newbiologicalsource <-
    gsub("Rhizoma alba", "\\1", data_rhizoma_alba$biologicalsource)
  data_rhizoma_alba$newbiologicalsource <-
    trimws(data_rhizoma_alba$newbiologicalsource)
  data_rhizoma_alba$newbiologicalsource <-
    capitalize(data_rhizoma_alba$newbiologicalsource)
  if (nrow(data_rhizoma_alba) != 0) {
    data_rhizoma_alba$newbiologicalsource <-
      paste(data_rhizoma_alba$newbiologicalsource, "rhizoma alba")
  }

  data_rhizoma_et_radix <- x %>%
    filter(grepl("^Rhizoma et radix", biologicalsource))
  data_rhizoma_et_radix$newbiologicalsource <-
    gsub(
      "Rhizoma et radix",
      "\\1",
      data_rhizoma_et_radix$biologicalsource
    )
  data_rhizoma_et_radix$newbiologicalsource <-
    trimws(data_rhizoma_et_radix$newbiologicalsource)
  data_rhizoma_et_radix$newbiologicalsource <-
    capitalize(data_rhizoma_et_radix$newbiologicalsource)
  if (nrow(data_rhizoma_et_radix) != 0) {
    data_rhizoma_et_radix$newbiologicalsource <-
      paste(
        data_rhizoma_et_radix$newbiologicalsource,
        "rhizoma et radix"
      )
  }

  data_semen <- x %>%
    filter(grepl("^Semen", biologicalsource))
  data_semen$newbiologicalsource <-
    gsub("Semen", "\\1", data_semen$biologicalsource)
  data_semen$newbiologicalsource <-
    trimws(data_semen$newbiologicalsource)
  data_semen$newbiologicalsource <-
    capitalize(data_semen$newbiologicalsource)
  if (nrow(data_semen) != 0) {
    data_semen$newbiologicalsource <-
      paste(data_semen$newbiologicalsource, "semen")
  }

  data_semen_germinatum <- x %>%
    filter(grepl("^Semen germinatum", biologicalsource))
  data_semen_germinatum$newbiologicalsource <-
    gsub(
      "Semen germinatum",
      "\\1",
      data_semen_germinatum$biologicalsource
    )
  data_semen_germinatum$newbiologicalsource <-
    trimws(data_semen_germinatum$newbiologicalsource)
  data_semen_germinatum$newbiologicalsource <-
    capitalize(data_semen_germinatum$newbiologicalsource)
  if (nrow(data_semen_germinatum) != 0) {
    data_semen_germinatum$newbiologicalsource <-
      paste(
        data_semen_germinatum$newbiologicalsource,
        "semen germinatum"
      )
  }

  data_spica <- x %>%
    filter(grepl("Spica ", biologicalsource, fixed = TRUE))
  data_spica$newbiologicalsource <-
    gsub("Spica", "\\1", data_spica$biologicalsource)
  data_spica$newbiologicalsource <-
    trimws(data_spica$newbiologicalsource)
  data_spica$newbiologicalsource <-
    capitalize(data_spica$newbiologicalsource)
  if (nrow(data_spica) != 0) {
    data_spica$newbiologicalsource <-
      paste(data_spica$newbiologicalsource, "spica")
  }

  data_stamen <- x %>%
    filter(grepl("^Stamen", biologicalsource))
  data_stamen$newbiologicalsource <-
    gsub("Stamen", "\\1", data_stamen$biologicalsource)
  data_stamen$newbiologicalsource <-
    trimws(data_stamen$newbiologicalsource)
  data_stamen$newbiologicalsource <-
    capitalize(data_stamen$newbiologicalsource)
  if (nrow(data_stamen) != 0) {
    data_stamen$newbiologicalsource <-
      paste(data_stamen$newbiologicalsource, "stamen")
  }

  data_stigma <- x %>%
    filter(grepl("Stigma ", biologicalsource, fixed = TRUE))
  data_stigma$newbiologicalsource <-
    gsub("Stigma", "\\1", data_stigma$biologicalsource)
  data_stigma$newbiologicalsource <-
    trimws(data_stigma$newbiologicalsource)
  data_stigma$newbiologicalsource <-
    capitalize(data_stigma$newbiologicalsource)
  if (nrow(data_stigma) != 0) {
    data_stigma$newbiologicalsource <-
      paste(data_stigma$newbiologicalsource, "stigma")
  }

  data_storax <- x %>%
    filter(grepl("^Storax", biologicalsource))
  data_storax$newbiologicalsource <-
    gsub("Storax", "\\1", data_storax$biologicalsource)
  data_storax$newbiologicalsource <-
    trimws(data_storax$newbiologicalsource)
  data_storax$newbiologicalsource <-
    capitalize(data_storax$newbiologicalsource)
  if (nrow(data_storax) != 0) {
    data_storax$newbiologicalsource <-
      paste(data_storax$newbiologicalsource, "storax")
  }

  data_thallus <- x %>%
    filter(grepl("^Thallus", biologicalsource, fixed = TRUE))
  data_thallus$newbiologicalsource <-
    gsub("Thallus", "\\1", data_thallus$biologicalsource)
  data_thallus$newbiologicalsource <-
    trimws(data_thallus$newbiologicalsource)
  data_thallus$newbiologicalsource <-
    capitalize(data_thallus$newbiologicalsource)
  if (nrow(data_thallus) != 0) {
    data_thallus$newbiologicalsource <-
      paste(data_thallus$newbiologicalsource, "thallus")
  }
  # not tuber

  x_2 <-
    rbind(
      data_bulbus,
      data_caulis,
      data_caluis_et_folium,
      data_corolla,
      data_cortex,
      data_exocarpium,
      data_exocarpium_rubrum,
      data_flos,
      data_folium,
      data_folium_et_cacumen,
      data_folium_et_caulis,
      data_fructus,
      data_fructus_germinatus,
      data_fructus_immaturus,
      data_fructus_retinervus,
      data_fructus_rotundus,
      data_herba,
      data_lignum,
      data_medulla,
      data_pericarpum,
      data_petiolus,
      data_pollen,
      data_radicis_cortex,
      data_radix,
      data_rhizoma_et_radix,
      data_radix_et_rhizoma,
      data_radix_preparata,
      data_ramulus,
      data_ramulus_cum_uncus,
      data_ramulus_et_folium,
      data_rhizoma,
      data_rhizoma_alba,
      data_rhizoma_et_radix,
      data_semen,
      data_semen_germinatum,
      data_spica,
      data_stamen,
      data_stigma,
      data_storax,
      data_thallus
    )

  x_3 <- left_join(x, x_2)

  x_3$newnewbiologicalsource <-
    apply(x_3, 1, function(x) {
      tail(na.omit(x), 1)
    })

  x_4 <- x_3 %>%
    select(-biologicalsource, -newbiologicalsource) %>%
    mutate(biologicalsource = newnewbiologicalsource) %>%
    select(-newnewbiologicalsource)

  return(x_4)
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
tcm_cleaning <- function(x) {
  data <- x

  data$newbiologicalsource <-
    gsub(" bulbus", "", data$latin, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" caulis", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" caulis et folium", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" corolla", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" cortex", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" exocarpium", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" exocarpium rubrum", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" flos", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" folium", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" folium et cacumen", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" folium et caulis", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" fructus", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" fructus germinatus", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" fructus immaturus", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" fructus retinervus", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" fructus rotundus", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" herba", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" lignum", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" medulla", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" pericarpum", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" petiolus", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" pollen", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" radicis cortex", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" radix", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" radix et rhizoma", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" radix preparata", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" ramulus", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" ramulus cum uncus", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" ramus et folium", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" rhizoma", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" rhizoma alba", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" rhizoma et radix", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" semen", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" semen germinatum", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" spica", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" stamen", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" stigma", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" thallus", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" et", "", data$newbiologicalsource, fixed = TRUE)

  # not tuber
  data_final <- data %>%
    mutate(latin = newbiologicalsource) %>%
    select(latin, common, biologicalsource)

  return(data_final)
}


#######################################################

#' Title
#'
#' @param data_selected
#' @param db
#' @param structure_field
#'
#' @return
#' @export
#'
#' @examples
standardizing_original <- function(data_selected,
                                   db,
                                   structure_field) {
  data_selected[setdiff(
    c(
      "name",
      "biologicalsource",
      "reference"
    ),
    names(data_selected)
  )] <- NA

  data_standard <- data.frame(data_selected) %>%
    mutate(database = db) %>%
    select(
      database,
      name,
      all_of(structure_field),
      biologicalsource,
      reference
    ) %>%
    distinct_at(vars(
      all_of(structure_field),
      biologicalsource
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

  return(data_standard)
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
preparing_name <- function(x) {
  x$nameSanitized <- x$structureOriginalNominal
  x$nameSanitized <- gsub("α", "alpha", x$nameSanitized)
  x$nameSanitized <- gsub("Α", "alpha", x$nameSanitized)
  x$nameSanitized <- gsub("β", "beta", x$nameSanitized)
  x$nameSanitized <- gsub("Β", "beta", x$nameSanitized)
  x$nameSanitized <- gsub("γ", "gamma", x$nameSanitized)
  x$nameSanitized <- gsub("Γ", "gamma", x$nameSanitized)
  x$nameSanitized <- gsub("δ", "delta", x$nameSanitized)
  x$nameSanitized <- gsub("Δ", "delta", x$nameSanitized)
  x$nameSanitized <- gsub("ε", "epsilon", x$nameSanitized)
  x$nameSanitized <- gsub("Ε", "epsilon", x$nameSanitized)
  x$nameSanitized <- gsub("- ", "-", x$nameSanitized)
  x$nameSanitized <- gsub("–", "-", x$nameSanitized)
  x$nameSanitized <- gsub("\\) ", "\\)", x$nameSanitized)
  x$nameSanitized <- trimws(x$nameSanitized)

  return(x)
}

#######################################################
#' Title
#'
#' @param x
#' @param no_rows_per_frame
#' @param path_to_store
#'
#' @return
#' @export
#'
#' @examples
split_data_table <- function(x, no_rows_per_frame, path_to_store) {
  split_vec <- seq(1, nrow(x), no_rows_per_frame)

  for (split_cut in split_vec) {
    sample <- x[split_cut:(split_cut + (no_rows_per_frame - 1))]
    write.table(
      sample,
      paste(
        path_to_store,
        "translatedOrganismGnfinderUntil_",
        as.integer(split_cut + (no_rows_per_frame - 1)),
        ".tsv",
        sep = ""
      ),
      row.names = FALSE,
      quote = FALSE,
      sep = "\t",
      fileEncoding = "UTF-8"
    )
  }
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
tcm_pharmacopoeia_cleaning <- function(x) {
  data <- x

  data$newbiologicalsource <-
    gsub("Bulbus ", "", data$organismTranslated, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Cacumen ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Caulis ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Corolla ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Cortex ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Exocarpium ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Exocarpium rubrum ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Flos ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Folium ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Folium ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Fructus ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Fructus germinatus ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Fructus immaturus ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Fructus retinervus ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Fructus rotundus ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Herba ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Lignum ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Medulla ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Pericarpum ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Petiolus ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Pollen ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Radicis cortex ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Radix ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Radix et rhizoma ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Ramulus ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Ramus ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Rhizoma ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Semen ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Semen germinatum ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Spica ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Stamen ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Stigma ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Thallus ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Uncis ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("et ", "", data$newbiologicalsource, fixed = TRUE)

  # not tuber
  data_final <- data %>%
    select(-organismTranslated) %>%
    mutate(organismTranslated = newbiologicalsource) %>%
    select(-newbiologicalsource)

  return(data_final)
}

#######################################################

#' Title
#'
#' @param num
#'
#' @return
#' @export
#'
#' @examples
gnfinder_cleaning <- function(num) {
  inpath_organism_f <- paste(
    pathTranslatedOrganismDistinct,
    "translatedOrganismGnfinderUntil_",
    num,
    ".tsv",
    sep = ""
  )

  inpath_gnfinder_f <-
    paste(
      pathSanitizedOrganismDirJson,
      "sanitizedOrganismGnfinderUntil_",
      num,
      ".json",
      sep = ""
    )

  outpath_f <-
    paste(pathSanitizedOrganismDirTsv,
      "sanitizedOrganismUntil_",
      num,
      ".tsv.zip",
      sep = ""
    )

  gnfound <- data.frame(fromJSON(
    txt = inpath_gnfinder_f,
    simplifyDataFrame = TRUE
  ))

  data_bio <- read_delim(
    file = inpath_organism_f,
    delim = "\t",
    escape_double = FALSE,
    trim_ws = FALSE
  ) %>%
    mutate_all(as.character)

  data_bio_clean <- biocleaning(
    x = gnfound,
    y = data_bio
  )

  write.table(
    x = data_bio_clean,
    file = gzfile(
      description = outpath_f,
      compression = 9,
      encoding = "UTF-8"
    ),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t",
    fileEncoding = "UTF-8"
  )
}

#######################################################

#' Title
#'
#' @param dfsel
#'
#' @return
#' @export
#'
#' @examples
taxo_cleaning_auto <- function(dfsel) {
  test <- dfsel %>%
    filter(!is.na(organismSanitized)) %>%
    distinct(organismSanitized,
      organism_database,
      .keep_all = TRUE
    )

  test$organism_1_kingdom <-
    y_as_na(
      x = test$organism_1_kingdom,
      y = "Not assigned"
    )

  test$organism_2_phylum <-
    y_as_na(
      x = test$organism_2_phylum,
      y = "Not assigned"
    )

  test$organism_3_class <-
    y_as_na(
      x = test$organism_3_class,
      y = "Not assigned"
    )

  test$organism_4_order <-
    y_as_na(
      x = test$organism_4_order,
      y = "Not assigned"
    )

  test$organism_5_family <-
    y_as_na(
      x = test$organism_5_family,
      y = "Not assigned"
    )

  test$organism_6_genus <-
    y_as_na(
      x = test$organism_6_genus,
      y = "Not assigned"
    )

  test$organism_7_species <-
    y_as_na(
      x = test$organism_7_species,
      y = "Not assigned"
    )

  species <- test %>%
    filter(!is.na(organism_7_species)) %>%
    arrange(match(
      x = organism_database,
      table = c(
        "Catalogue of Life",
        "NCBI"
      )
    )) %>%
    group_by(
      organism_1_kingdom,
      organism_7_species
    ) %>%
    mutate(
      organism_2_phylum = na.locf(
        object = organism_2_phylum,
        na.rm = FALSE,
        fromLast = TRUE
      ),
      organism_3_class = na.locf(
        object = organism_3_class,
        na.rm = FALSE,
        fromLast = TRUE
      ),
      organism_4_order = na.locf(
        object = organism_4_order,
        na.rm = FALSE,
        fromLast = TRUE
      ),
      organism_5_family = na.locf(
        object = organism_5_family,
        na.rm = FALSE,
        fromLast = TRUE
      ),
      organism_6_genus = na.locf(
        object = organism_6_genus,
        na.rm = FALSE,
        fromLast = TRUE
      ),
    ) %>%
    ungroup() %>%
    distinct(organism_7_species,
      .keep_all = TRUE
    ) %>%
    select(-organismOriginal, -organismTranslated, -organismSanitized)

  species_fill <- test %>%
    filter(!is.na(organism_7_species)) %>%
    select(
      organismOriginal,
      organismTranslated,
      organismSanitized,
      organism_7_species
    )

  species_full <- left_join(species_fill, species)

  unspecies <- test %>%
    filter(is.na(organism_7_species))

  genus_1 <- rbind(species_full, unspecies)

  genus <- genus_1 %>%
    filter(!is.na(organism_6_genus)) %>%
    arrange(match(
      x = organism_database,
      table = c(
        "Catalogue of Life",
        "NCBI"
      )
    )) %>%
    group_by(
      organism_1_kingdom,
      organism_6_genus
    ) %>%
    mutate(
      organism_2_phylum = na.locf(
        object = organism_2_phylum,
        na.rm = FALSE,
        fromLast = TRUE
      ),
      organism_3_class = na.locf(
        object = organism_3_class,
        na.rm = FALSE,
        fromLast = TRUE
      ),
      organism_4_order = na.locf(
        object = organism_4_order,
        na.rm = FALSE,
        fromLast = TRUE
      ),
      organism_5_family = na.locf(
        object = organism_5_family,
        na.rm = FALSE,
        fromLast = TRUE
      )
    ) %>%
    ungroup() %>%
    distinct(organism_6_genus,
      .keep_all = TRUE
    ) %>%
    select(
      -organismOriginal,
      -organismTranslated,
      -organismSanitized,
      -organism_7_species
    )

  genus_fill <- genus_1 %>%
    filter(!is.na(organism_6_genus)) %>%
    select(
      organismOriginal,
      organismTranslated,
      organismSanitized,
      organism_7_species,
      organism_6_genus
    )

  genus_full <- left_join(genus_fill, genus)

  ungenus <- genus_1 %>%
    filter(is.na(organism_6_genus))

  family_1 <- rbind(genus_full, ungenus)

  family <- family_1 %>%
    filter(!is.na(organism_5_family)) %>%
    arrange(match(
      x = organism_database,
      table = c(
        "Catalogue of Life",
        "NCBI"
      )
    )) %>%
    group_by(
      organism_1_kingdom,
      organism_5_family
    ) %>%
    mutate(
      organism_2_phylum = na.locf(
        object = organism_2_phylum,
        na.rm = FALSE,
        fromLast = TRUE
      ),
      organism_3_class = na.locf(
        object = organism_3_class,
        na.rm = FALSE,
        fromLast = TRUE
      ),
      organism_4_order = na.locf(
        object = organism_4_order,
        na.rm = FALSE,
        fromLast = TRUE
      )
    ) %>%
    ungroup() %>%
    distinct(organism_5_family,
      .keep_all = TRUE
    ) %>%
    select(
      -organismOriginal,
      -organismTranslated,
      -organismSanitized,
      -organism_7_species,
      -organism_6_genus
    )

  family_fill <- family_1 %>%
    filter(!is.na(organism_5_family)) %>%
    select(
      organismOriginal,
      organismTranslated,
      organismSanitized,
      organism_7_species,
      organism_6_genus,
      organism_5_family
    )

  family_full <- left_join(family_fill, family)

  unfamily <- family_1 %>%
    filter(is.na(organism_5_family))

  order_1 <- rbind(family_full, unfamily)

  order <- order_1 %>%
    filter(!is.na(organism_4_order)) %>%
    arrange(match(
      x = organism_database,
      table = c(
        "Catalogue of Life",
        "NCBI"
      )
    )) %>%
    group_by(
      organism_1_kingdom,
      organism_4_order
    ) %>%
    mutate(
      organism_2_phylum = na.locf(
        object = organism_2_phylum,
        na.rm = FALSE,
        fromLast = TRUE
      ),
      organism_3_class = na.locf(
        object = organism_3_class,
        na.rm = FALSE,
        fromLast = TRUE
      )
    ) %>%
    ungroup() %>%
    ungroup() %>%
    distinct(organism_4_order,
      .keep_all = TRUE
    ) %>%
    select(
      -organismOriginal,
      -organismTranslated,
      -organismSanitized,
      -organism_7_species,
      -organism_6_genus,
      -organism_5_family
    )

  order_fill <- order_1 %>%
    filter(!is.na(organism_4_order)) %>%
    select(
      organismOriginal,
      organismTranslated,
      organismSanitized,
      organism_7_species,
      organism_6_genus,
      organism_5_family,
      organism_4_order
    )

  order_full <- left_join(order_fill, order)

  unorder <- order_1 %>%
    filter(is.na(organism_4_order))

  class_1 <- rbind(order_full, unorder)

  class <- class_1 %>%
    filter(!is.na(organism_3_class)) %>%
    arrange(match(
      x = organism_database,
      table = c(
        "Catalogue of Life",
        "NCBI"
      )
    )) %>%
    group_by(
      organism_1_kingdom,
      organism_3_class
    ) %>%
    mutate(organism_2_phylum = na.locf(
      object = organism_2_phylum,
      na.rm = FALSE,
      fromLast = TRUE
    )) %>%
    ungroup() %>%
    distinct(organism_3_class,
      .keep_all = TRUE
    ) %>%
    select(
      -organismOriginal,
      -organismTranslated,
      -organismSanitized,
      -organism_7_species,
      -organism_6_genus,
      -organism_5_family,
      -organism_4_order,
    )

  class_fill <- class_1 %>%
    filter(!is.na(organism_3_class)) %>%
    select(
      organismOriginal,
      organismTranslated,
      organismSanitized,
      organism_7_species,
      organism_6_genus,
      organism_5_family,
      organism_4_order,
      organism_3_class
    )

  class_full <- left_join(class_fill, class)

  unclass <- class_1 %>%
    filter(is.na(organism_3_class))

  phylum_1 <- rbind(class_full, unclass)

  phylum <- phylum_1 %>%
    filter(!is.na(organism_2_phylum)) %>%
    arrange(match(
      x = organism_database,
      table = c(
        "Catalogue of Life",
        "NCBI"
      )
    )) %>%
    group_by(organism_2_phylum) %>%
    mutate(organism_1_kingdom = na.locf(
      object = organism_1_kingdom,
      na.rm = FALSE,
      fromLast = TRUE
    )) %>%
    ungroup() %>%
    distinct(organism_2_phylum,
      .keep_all = TRUE
    ) %>%
    select(
      -organismOriginal,
      -organismTranslated,
      -organismSanitized,
      -organism_7_species,
      -organism_6_genus,
      -organism_5_family,
      -organism_4_order,
      -organism_3_class
    )

  phylum_fill <- phylum_1 %>%
    filter(!is.na(organism_2_phylum)) %>%
    select(
      organismOriginal,
      organismTranslated,
      organismSanitized,
      organism_7_species,
      organism_6_genus,
      organism_5_family,
      organism_4_order,
      organism_3_class,
      organism_2_phylum
    )

  phylum_full <- left_join(phylum_fill, phylum)

  unphylum <- phylum_1 %>%
    filter(is.na(organism_2_phylum))

  kingdom_1 <- rbind(phylum_full, unphylum)

  kingdom_tojoin <- kingdom_1 %>%
    select(-organismOriginal, -organismTranslated)

  tojoin <- dfsel

  newdf <- left_join(tojoin,
    kingdom_tojoin,
    by = c("organismSanitized" = "organismSanitized")
  ) %>%
    mutate(organism_modified_taxonomy_auto = if_else(
      condition = paste(
        organism_1_kingdom.x,
        organism_2_phylum.x,
        organism_3_class.x,
        organism_4_order.x,
        organism_5_family.x,
        organism_6_genus.x,
        organism_7_species.x
      ) ==
        paste(
          organism_1_kingdom.y,
          organism_2_phylum.y,
          organism_3_class.y,
          organism_4_order.y,
          organism_5_family.y,
          organism_6_genus.y,
          organism_7_species.y
        ),
      true = "y",
      false = ""
    )) %>%
    select(
      organismOriginal,
      organismTranslated,
      organismSanitized,
      organism_database = organism_database.y,
      organism_1_kingdom = organism_1_kingdom.y,
      organism_2_phylum = organism_2_phylum.y,
      organism_3_class = organism_3_class.y,
      organism_4_order = organism_4_order.y,
      organism_5_family = organism_5_family.y,
      organism_6_genus = organism_6_genus.y,
      organism_7_species = organism_7_species.y,
      organism_modified_taxonomy_auto
    ) %>%
    distinct(organismOriginal,
      organismTranslated,
      organismSanitized,
      .keep_all = TRUE
    )

  newdf$organism_modified_taxonomy_auto <-
    y_as_na(
      x = newdf$organism_modified_taxonomy_auto,
      y = ""
    )

  return(newdf)
}

#######################################################

#' Title
#'
#' @param X
#'
#' @return
#' @export
#'
#' @examples
getref <- function(X) {
  tryCatch(
    {
      cr_works(
        query = X,
        sort = "score",
        order = "desc",
        limit = 1
      )
    },
    error = function(e) {
      NA
    }
  )
}

#######################################################

#' Title
#'
#' @param X
#'
#' @return
#' @export
#'
#' @examples
getrefPubmed <- function(X) {
  tryCatch(
    {
      df <- entrez_summary(db = "pubmed", id = X)

      translatedDoi <-
        ifelse(test = "doi" %in% df[["articleids"]][, 1],
          yes = trimws(df[["articleids"]][["value"]][[which(df[["articleids"]] == "doi")]]),
          no = NA
        )
      translatedJournal <- df[["fulljournalname"]]
      translatedTitle <- df[["title"]]
      translatedAuthor <- df[["sortfirstauthor"]]
      translatedDate <- df[["pubdate"]]

      ids <-
        data.frame(
          translatedDoi,
          translatedJournal,
          translatedTitle,
          translatedAuthor,
          translatedDate
        )
      return(ids)
    },
    error = function(e) {
      data.frame(
        "translatedDoi" = NA,
        "translatedJournal" = NA,
        "translatedTitle" = NA,
        "translatedAuthor" = NA,
        "translatedDate" = NA
      )
    }
  )
}

#######################################################

#' Title
#'
#' @param X
#'
#' @return
#' @export
#'
#' @examples
getrefDoi <- function(X) {
  tryCatch(
    {
      cr_works(dois = X)
    },
    error = function(e) {
      NA
    }
  )
}

#######################################################
