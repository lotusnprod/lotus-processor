source("paths.R")
library(data.table)
library(stringi)
source("r/homemadeShift.R")

#' Title
#'
#' @param gnfound
#' @param names
#' @param names_quotes
#' @param organismCol
#'
#' @return
#' @export
#'
#' @examples
biocleaning <- function(gnfound, names, names_quotes, organismCol) {
  log_debug("Biocleaning")
  log_debug("Biocleaning: finished creating dataframe")
  # extracting preferred results data table
  ## as list of dataframes
  df2 <- gnfound$names.verification$preferredResults

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
      )),
      variety = sum(as.numeric(
        stri_detect(
          str = classificationRank,
          fixed = c(
            "variety",
            "varietas",
            "var"
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
  df4$variety[df4$variety >= 1] <- 1

  df4[setdiff(
    x = c("isSynonym"),
    y = names(df4)
  )] <- NA

  # the synonym part is there to avoid the (actually)
  ## non-optimal output from Catalogue of Life in GNFinder
  ### (explained in https://github.com/gnames/gnfinder/issues/48)
  df5a <- df4 %>%
    mutate(
      n = rowSums(.[c(
        "kingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
        "variety"
      )]),
      id = as.integer(id)
    ) %>%
    group_by(id) %>%
    arrange(desc(n), !is.na(isSynonym)) %>%
    ungroup() %>%
    distinct(id,
      .keep_all = TRUE
    ) %>%
    arrange(id) %>%
    select(id)

  df5b <- df4 %>%
    mutate(
      n = rowSums(.[c(
        "kingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
        "variety"
      )]),
      id = as.integer(id)
    ) %>%
    group_by(id) %>%
    arrange(desc(n), !is.na(isSynonym)) %>%
    ungroup() %>%
    arrange(id)

  df5b[setdiff(
    x = c("classificationIds"),
    y = names(df5b)
  )] <- NA

  df6a <- bind_cols(df5a, rows)
  df6b <- left_join(df6a, df5b) %>%
    filter(!is.na(classificationIds))

  if (nrow(df6b) == 0) {
    df6b[1, colnames(df6b)] <- NA
  }

  # adding row number
  df7 <- gnfound$names.start %>%
    data.table() %>%
    mutate(nrow = row_number())

  colnames(df7)[1] <- "sum"

  # joining
  taxo <- right_join(df6b, df7)

  taxo[setdiff(
    x = c(
      "matchedCanonicalFull",
      "currentCanonicalFull",
      "taxonId",
      "dataSourceTitle",
      "classificationPath",
      "classificationRank",
      "classificationIds"
    ),
    y = names(taxo)
  )] <- NA

  taxo <- taxo %>%
    select(
      canonicalname = matchedCanonicalFull,
      canonicalnameCurrent = currentCanonicalFull,
      taxonId,
      dbTaxo = dataSourceTitle,
      taxonomy = classificationPath,
      rank = classificationRank,
      ids = classificationIds,
      sum
    )

  log_debug("Biocleaning: finished joining")

  dbQuality <- gnfound$names.verification$dataSourceQuality
  dbTaxo <- gnfound$names.verification$bestResult$dataSourceTitle

  dfQuality <- data.frame(dbTaxo, dbQuality) %>%
    distinct(dbTaxo, .keep_all = TRUE)

  taxoEnhanced <- left_join(taxo, dfQuality)

  # computing sum of characters to match with GNFinder results
  if (organismCol == "organismOriginal") {
    names_quotes$nchar <-
      nchar(x = names_quotes$`"organismOriginal"`)
  }

  if (organismCol == "organismInterim") {
    names_quotes$nchar <-
      nchar(x = names_quotes$`"organismInterim"`)
  }

  names_quotes[1, "sum"] <- nchar(colnames(names_quotes)[1]) + 1
  for (i in 2:nrow(names_quotes)) {
    names_quotes[i, "sum"] <- names_quotes[i - 1, "nchar"] + 1 + names_quotes[i - 1, "sum"]
  }

  # adding min and max to merge
  taxoEnhanced <- taxoEnhanced %>%
    mutate(
      value_min = sum,
      value_max = sum
    ) %>%
    data.table()

  # filtering non-empty values
  y_2 <- names_quotes %>%
    mutate(value_min = sum)

  # filling sum values
  y_2$value_min <- as.numeric(y_2$value_min)
  y_2$value_max <- homemadeShift(y_2$sum, 1) - 1
  y_2[nrow(y_2), 5] <- y_2[nrow(y_2), 4] + 10000

  # transforming as data table (needed for next function)
  y_2 <- y_2 %>%
    bind_cols(., data_bio) %>%
    select(organismOriginal, nchar, sum, value_min, value_max) %>%
    data.table()

  # setting joining keys
  setkey(taxoEnhanced, value_min, value_max)
  setkey(y_2, value_min, value_max)
  # joining
  pre_final_db <- foverlaps(
    taxoEnhanced,
    na.omit(y_2)
  )

  # selecting
  final_db <- left_join(
    names,
    pre_final_db
  ) %>%
    select(
      -i.sum,
      -i.value_max,
      -i.value_min
    )

  return(final_db)
}
