#######################################################
####################   Functions   ####################
#######################################################

library(Hmisc)
library(plotly)
library(ChemmineR)
library(chorddiag)
library(collapsibleTree)
library(data.table)
library(digest)
library(dplyr)
library(eulerr)
library(ggraph)
library(igraph)
library(jsonlite)
library(parallel)
library(pbmcapply)
library(purrr)
library(RColorBrewer)
library(rcrossref)
library(readr)
library(readxl)
library(rentrez)
library(reticulate)
library(rvest)
library(splitstackshape)
library(stringi)
library(stringr)
library(taxize)
library(tidyverse)
library(UpSetR)
library(webchem)
library(XML)
library(zoo)
library(tidyr)

#######################################################
#######################################################

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
#######################################################

pubchem2inchi <- function(i)
{
  tryCatch({
    cpd <-
      data_translated_pubchem[i, "structure_original_numerical_pubchem"]
    url <-
      paste(
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
        cpd,
        "/property/InChI/txt",
        sep = ""
      )
    url <- gsub(pattern = "\\s",
                replacement = "%20",
                x = url)
    read_html(url) %>%
      html_text()
  }
  , error = function(e) {
    NA
  })
}

#######################################################
#######################################################

shift <- function(x, n) {
  c(x[-(seq(n))], rep(NA, n))
}

#######################################################
#######################################################

biocleaning <- function(x, y, organismCol)
{
  #selecting and adding row number
  df <- x %>%
    select(names.verbatim,
           names.verification) %>%
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
  
  #extracting preferred results data table
  ##as list of dataframes
  df2 <- x$names.verification$preferredResults
  
  #outputting row numbers
  rows <- df2 %>%
    data.table() %>%
    mutate(nrow = row_number()) %>%
    filter(. != "NULL") %>% select(nrow)
  
  ##as dataframe and adding row number
  df3 <- bind_rows(df2,
                   .id = "id")
  
  #selecting best result (with best score and best filled taxonomy)
  df4 <- df3 %>%
    rowwise() %>%
    mutate(
      kingdom = sum(as.numeric(
        stri_detect(
          str = classificationRank,
          fixed = c("kingdom",
                    "Kingdom",
                    "regn.")
        )
      )),
      phylum =  sum(as.numeric(
        stri_detect(str = classificationRank,
                    fixed = c("phylum",
                              "Phylum",
                              "phyl."))
      )),
      class =  sum(as.numeric(
        stri_detect(str = classificationRank,
                    fixed = c("class",
                              "Class",
                              "cl."))
      )),
      order =  sum(as.numeric(
        stri_detect(str = classificationRank,
                    fixed = c("order",
                              "Order",
                              "ord."))
      )),
      family = sum(as.numeric(
        stri_detect(str = classificationRank,
                    fixed = c("family",
                              "Family",
                              "fam."))
      )),
      genus =  sum(as.numeric(
        stri_detect(str = classificationRank,
                    fixed = c("genus",
                              "Genus"))
      )),
      species = sum(as.numeric(
        stri_detect(
          str = classificationRank,
          fixed = c("species",
                    "Species",
                    "spec.",
                    "sp.")
        )
      )),
      variety = sum(as.numeric(
        stri_detect(
          str = classificationRank,
          fixed = c("variety",
                    "varietas",
                    "var")
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
  
  #the synonym part is there to avoid the (actually)
  ##non-optimal output from Catalogue of Life in GNFinder
  ###(explained in https://github.com/gnames/gnfinder/issues/48)
  df5 <- df4 %>%
    mutate(n = rowSums(.[c("kingdom",
                           "phylum",
                           "class",
                           "order",
                           "family",
                           "genus",
                           "species",
                           "variety")])) %>%
    group_by(id) %>%
    arrange(desc(n), !is.na(isSynonym)) %>%
    ungroup() %>%
    distinct(id,
             .keep_all = TRUE) %>%
    arrange(as.numeric(id))
  
  df6 <- cbind(df5, rows)
  
  #adding row number
  df7 <- x$names.start %>%
    data.table() %>%
    mutate(nrow = row_number())
  
  colnames(df7)[1] <- "sum"
  
  #joining
  taxo <- right_join(df6, df7) %>%
    select(
      canonicalname = matchedCanonicalFull,
      taxonId,
      dbTaxo = dataSourceTitle,
      taxonomy = classificationPath,
      rank = classificationRank,
      sum
    )
  
  dbQuality <- x$names.verification$dataSourceQuality
  dbTaxo <- x$names.verification$bestResult$dataSourceTitle
  
  dfQuality <- data.frame(dbTaxo, dbQuality) %>%
    distinct(dbTaxo, .keep_all = TRUE)
  
  taxoEnhanced <- left_join(taxo, dfQuality)
  
  #computing sum of characters to match with GNFinder results
  if (organismCol == "organismOriginal")
    y$nchar <-
    nchar(x = y$organismOriginal)
  
  if (organismCol == "organismInterim")
  y$nchar <-
    nchar(x = y$organismInterim)
  
  y[1, "sum"] <- nchar(colnames(y)[1]) + 1
  for (i in 2:nrow(y)) {
    y[i, "sum"] <- y[i - 1, "nchar"] + 1 + y[i - 1, "sum"]
  }
  
  #adding min and max to merge
  taxoEnhanced <- taxoEnhanced %>%
    mutate(value_min = sum,
           value_max = sum) %>%
    data.table()
  
  #filtering non-empty values
  y_2 <- y %>%
    mutate(value_min = sum)
  
  #filling sum values
  y_2$value_min <- as.numeric(y_2$value_min)
  y_2$value_max <- shift(y_2$sum, 1) - 1
  y_2[nrow(y_2), 5] <- y_2[nrow(y_2), 4] + 10000
  
  #transforming as data table (needed for next function)
  y_2 <- y_2 %>%
    data.table()
  
  #setting joining keys
  setkey(taxoEnhanced, value_min, value_max)
  setkey(y_2, value_min, value_max)
  
  #joining
  pre_final_db <- foverlaps(taxoEnhanced,
                            y_2)
  
  #selecting
  final_db <- left_join(y,
                        pre_final_db) %>%
    select(# -namesverbatim,-nchar,-sum,
      -value_max,
      -value_min,
      -i.sum,
      -i.value_max,
      -i.value_min)
  
  return(final_db)
}

#######################################################
#######################################################

manipulating_taxo <- function(dfsel, dic) {
  #creating variables for replacement by dictionary
  a <- paste("\\b", dic$taxaLevel, "\\b", sep = "")
  b <- dic$taxaLevelStandard
  
  #selecting and splitting taxonomy and ranks
  df1 <- dfsel %>%
    select(identifier = 1,
           canonicalname,
           dbTaxo,
           taxonomy,
           rank) %>%
    cSplit(splitCols = "taxonomy",
           sep = "|") %>%
    cSplit(splitCols = "rank",
           sep = "|") %>%
    lapply(as.character) %>%
    as_tibble()
  
  #manipulating taxa
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
             .keep_all = TRUE)
  
  df2$rank <- stri_replace_all_regex(
    str = df2$rank,
    pattern = a,
    replacement = b,
    case_insensitive = FALSE,
    vectorize_all = FALSE
  )
  
  #removing false non-empty cells
  df2$identifier <- y_as_na(x = df2$identifier,
                            y = "")
  
  df2$rank <- y_as_na(x = df2$rank,
                      y = "")
  
  df2$taxonomy <- y_as_na(x = df2$taxonomy,
                          y = "")
  
  df2$identifier <- y_as_na(x = df2$identifier,
                            y = "NA NA")
  
  df2$rank <- y_as_na(x = df2$rank,
                      y = "")
  
  df2$rank <- ifelse(test = is.na(df2$rank),
                     yes = "NA",
                     no = df2$rank)
  
  colnames(df2)[3] <- "dbTaxo"
  
  #manipulating taxa
  df3 <- df2 %>%
    pivot_wider(names_from = rank,
                values_from = taxonomy) %>%
    select_if(
      names(.) %in%
        c(
          "identifier",
          "canonicalname",
          "dbTaxo",
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
  
  #pasting suffix to colnames to pivot then (the double pivot allows to tidy the data)
  colnames(df3)[4:ncol(df3)] <-
    paste("bio_", colnames(df3)[4:ncol(df3)], sep = "")
  
  #pivoting (long)
  df4 <- df3 %>%
    pivot_longer(
      cols = 4:ncol(.),
      names_to = c(".value", "level"),
      names_sep = "_",
      values_to = "taxonomy",
      values_drop_na = TRUE
    )
  
  #pivoting (wide)
  df5 <- df4 %>%
    group_by(canonicalname) %>%
    distinct(level, .keep_all = TRUE) %>%
    pivot_wider(names_from = level,
                values_from = bio) %>%
    select_if(
      names(.) %in%
        c(
          "canonicalname",
          "dbTaxo",
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
  
  #adding taxa to initial df
  df6 <- left_join(dfsel, df5)
  
  return(df6)
}

#######################################################
#######################################################

# biofilling <- function(x)
# {
#   #mutating to char to avoid mismatches
#   x <- x %>%
#     mutate_all(as.character)
#   
#   #avoiding false non-empty taxonomies
#   x$taxonomy <- y_as_na(x = x$taxonomy,
#                         y = "")
#   
#   #filtering taxonomies to fill
#   df <- x
#   
#   #creating lines if none to fill
#   if (nrow(df) == 0)
#     x[nrow(x) + 1, 1] <- "NA na"
#   
#   if (nrow(df) == 0)
#     df[nrow(df) + 1, 1] <- "NA na"
#   
#   if (nrow(df) == 0)
#     x[nrow(x) + 1, 2] <- "NA na"
#   
#   if (nrow(df) == 0)
#     df[nrow(df) + 1, 2] <- "NA na"
#   
#   #manipulating taxonomies (see function above)
#   df2 <- manipulating_taxo(dfsel = df,
#                            dic = taxaRanksDictionary)
#   
#   #filtering non-empty taxonomies
#   df2_full <- df2 %>%
#     filter(
#       !is.na(kingdom) |
#         !is.na(phylum) |
#         !is.na(class) |
#         !is.na(order) |
#         !is.na(family) |
#         !is.na(genus)
#     ) %>%
#     select(-rank, -taxonomy)
#   
#   #filtering empty taxonomies
#   df3 <- df2 %>%
#     filter(
#       is.na(kingdom) &
#         is.na(phylum) &
#         is.na(class) &
#         is.na(order) &
#         is.na(family) &
#         is.na(genus)
#     ) %>%
#     distinct(organismTranslated,
#              canonicalname) %>%
#     mutate_at(
#       .vars = vars(canonicalname),
#       .funs = function(x) {
#         gsub(pattern = "[^ -~]",
#              replacement = "",
#              x = x)
#       }
#     ) #because of problems otherwise (like Abelia  ×  grandiflora)
#   
#   #outputting distinct empty taxonomies to avoid iterations
#   df3_dis <- df3 %>%
#     filter(!is.na(canonicalname)) %>%
#     distinct(canonicalname)
#   
#   if (nrow(df3_dis) == 0)
#     df3_dis[nrow(df3_dis) + 1, 1] <- "NA na"
#   
#   if (nrow(df3_dis) == 0)
#     df3_dis[nrow(df3_dis) + 1, 1] <- "NA na"
#   
#   if (nrow(df3_dis) == 0)
#     df3_dis[nrow(df3_dis) + 1, 2] <- "NA na"
#   
#   if (nrow(df3_dis) == 0)
#     df3_dis[nrow(df3_dis) + 1, 2] <- "NA na"
#   
#   #running gnresolve on identifiers with empty taxonomies
#   df4 <- gnr_resolve(
#     names = df3_dis$canonicalname,
#     data_source_ids = NULL,
#     resolve_once = FALSE,
#     with_context = TRUE,
#     canonical = TRUE,
#     cap_first = FALSE,
#     best_match_only = FALSE,
#     http = "post",
#     fields = "all"
#   )
#   
#   #selecting best result (with best score and best filled taxonomy)
#   df5 <- df4 %>%
#     rowwise() %>%
#     mutate(
#       kingdom = sum(as.numeric(
#         stri_detect(
#           str = classification_path_ranks,
#           fixed = c("kingdom", "Kingdom", "regn.")
#         )
#       )),
#       phylum =  sum(as.numeric(
#         stri_detect(str = classification_path_ranks,
#                     fixed = c("phylum", "Phylum", "phyl."))
#       )),
#       class =  sum(as.numeric(
#         stri_detect(str = classification_path_ranks,
#                     fixed = c("class", "Class", "cl."))
#       )),
#       order =  sum(as.numeric(
#         stri_detect(str = classification_path_ranks,
#                     fixed = c("order", "Order", "ord."))
#       )),
#       family = sum(as.numeric(
#         stri_detect(str = classification_path_ranks,
#                     fixed = c("family", "Family", "fam."))
#       )),
#       genus =  sum(as.numeric(
#         stri_detect(str = classification_path_ranks,
#                     fixed = c("genus", "Genus"))
#       )),
#       species = sum(as.numeric(
#         stri_detect(
#           str = classification_path_ranks,
#           fixed = c("species", "Species", "spec.", "sp.")
#         )
#       )),
#       variety = sum(as.numeric(
#         stri_detect(
#           str = classification_path_ranks,
#           fixed = c("variety",
#                     "varietas",
#                     "var")
#         )
#       ))
#     ) %>%
#     ungroup()
#   
#   df5$kingdom[df5$kingdom >= 1] <- 1
#   df5$phylum[df5$phylum >= 1] <- 1
#   df5$class[df5$class >= 1] <- 1
#   df5$order[df5$order >= 1] <- 1
#   df5$family[df5$family >= 1] <- 1
#   df5$genus[df5$genus >= 1] <- 1
#   df5$species[df5$species >= 1] <- 1
#   df5$variety[df5$variety >= 1] <- 1
#   
#   df5 <- df5 %>%
#     mutate(n = rowSums(.[c("kingdom",
#                            "phylum",
#                            "class",
#                            "order",
#                            "family",
#                            "genus",
#                            "species",
#                            "variety")])) %>%
#     group_by(user_supplied_name) %>%
#     arrange(desc(score), desc(n)) %>%
#     ungroup() %>%
#     distinct(user_supplied_name,
#              .keep_all = TRUE) %>%
#     select(
#       user_supplied_name,
#       canonicalname = matched_name2,
#       taxonId = taxon_id,
#       dbTaxo = data_source_title,
#       taxonomy = classification_path,
#       rank = classification_path_ranks
#     )
#   
#   #manipulating taxa
#   df6 <- manipulating_taxo(dfsel = df5,
#                            dic = taxaRanksDictionary) %>%
#     select(-user_supplied_name)
#   
#   #joining
#   df7 <-
#     left_join(df3,
#               df6)
#   
#   #harmonizing columns (for rbind)
#   df7[setdiff(
#     x = c(
#       "organismTranslated",
#       "canonicalname",
#       "dbTaxo",
#       "kingdom",
#       "phylum",
#       "class",
#       "order",
#       "family",
#       "genus",
#       "species",
#       "variety",
#       "nchar", 
#       "sum"
#     ),
#     y = names(df7)
#   )] <- NA
#   
#   #filtering non-empty taxonomies
#   df7_full <- df7 %>%
#     filter(
#       !is.na(kingdom) |
#         !is.na(phylum) |
#         !is.na(class) |
#         !is.na(order) |
#         !is.na(family) |
#         !is.na(genus)
#     ) %>%
#     select(-rank, -taxonomy)
#   
#   #filtering empty taxonomies
#   df8 <- df7 %>%
#     filter(
#       is.na(kingdom) &
#         is.na(phylum) &
#         is.na(class) &
#         is.na(order) &
#         is.na(family) &
#         is.na(genus)
#     ) %>%
#     select(organismTranslated) %>%
#     mutate_at(
#       .vars = vars(organismTranslated),
#       .funs = function(x) {
#         gsub(pattern = "[^ -~]",
#              replacement =  "",
#              x =  x)
#       }
#     ) %>%
#     mutate_at(
#       .vars = vars(organismTranslated),
#       .funs = function(x) {
#         gsub(pattern = ";",
#              replacement = "",
#              x = x)
#       }
#     ) %>%
#     mutate(query = organismTranslated)
#   
#   #removing disturbing words
#   ##creating variables for replacement by dictionary
#   c <- paste("\\b", blacklistDictionary$blackName, "\\b", sep = "")
#   d <- blacklistDictionary$replacement
#   
#   df8$query <-
#     stri_replace_all_regex(
#       str = df8$query,
#       pattern = c,
#       replacement = d,
#       case_insensitive = TRUE,
#       vectorize_all = FALSE
#     )
#   
#   df8$query <- gsub(pattern = "FALSE",
#                     replacement = "",
#                     x = df8$query)
#   
#   df8$query <- trimws(df8$query)
#   
#   #outputting distinct empty taxonomies to avoid iterations
#   df8_dis <- df8 %>%
#     distinct(query)
#   
#   #running gnresolve on whole identifier cleaned
#   ##if less than 500 empty entries
#   if (nrow(df8_dis) <= 500)
#     df9 <- gnr_resolve(
#       names = df8_dis$query,
#       data_source_ids = NULL,
#       resolve_once = FALSE,
#       with_context = TRUE,
#       canonical = TRUE,
#       cap_first = FALSE,
#       best_match_only = FALSE,
#       http = "post",
#       fields = "all"
#     )
#   
#   ##if more than 500 empty entries slicing in multiple pieces because of gnr_resolve structure
#   ###if more than 1000 empty entries
#   if (nrow(df8_dis) > 1000)
#     df8_dis_a <- slice(df8_dis, 1:500)
#   
#   if (nrow(df8_dis) > 1000)
#     df8_dis_b <- slice(df8_dis, 501:1000)
#   
#   if (nrow(df8_dis) > 1000)
#     df8_dis_c <- slice(df8_dis, 1001:nrow(df8_dis))
#   
#   if (nrow(df8_dis) > 1000)
#     df9_a <- gnr_resolve(
#       names = df8_dis_a$query,
#       data_source_ids = NULL,
#       resolve_once = FALSE,
#       with_context = TRUE,
#       canonical = TRUE,
#       cap_first = FALSE,
#       best_match_only = FALSE,
#       http = "post",
#       fields = "all"
#     )
#   
#   if (nrow(df8_dis) > 1000)
#     df9_b <- gnr_resolve(
#       names = df8_dis_b$query,
#       data_source_ids = NULL,
#       resolve_once = FALSE,
#       with_context = TRUE,
#       canonical = TRUE,
#       cap_first = FALSE,
#       best_match_only = FALSE,
#       http = "post",
#       fields = "all"
#     )
#   
#   if (nrow(df8_dis) > 1000)
#     df9_c <- gnr_resolve(
#       names = df8_dis_c$query,
#       data_source_ids = NULL,
#       resolve_once = FALSE,
#       with_context = TRUE,
#       canonical = TRUE,
#       cap_first = FALSE,
#       best_match_only = FALSE,
#       http = "post",
#       fields = "all"
#     )
#   
#   if (nrow(df8_dis) > 1000)
#     df9_a[setdiff(x = names(df4), y = names(df9_a))] <- NA
#   
#   if (nrow(df8_dis) > 1000)
#     df9_b[setdiff(x = names(df4), y = names(df9_b))] <- NA
#   
#   if (nrow(df8_dis) > 1000)
#     df9_c[setdiff(x = names(df4), y = names(df9_c))] <- NA
#   
#   if (nrow(df8_dis) > 1000)
#     df9 <- rbind(df9_a, df9_b, df9_c)
#   
#   ###if between 500 and 1000 empty entries
#   if (nrow(df8_dis) <= 1000 &
#       nrow(df8_dis) > 500)
#     df8_dis_a <- slice(df8_dis, 1:500)
#   
#   if (nrow(df8_dis) <= 1000 &
#       nrow(df8_dis) > 500)
#     df8_dis_b <- slice(df8_dis, 501:nrow(df8_dis))
#   
#   if (nrow(df8_dis) <= 1000 &
#       nrow(df8_dis) > 500)
#     df9_a <- gnr_resolve(
#       names = df8_dis_a$query,
#       data_source_ids = NULL,
#       resolve_once = FALSE,
#       with_context = TRUE,
#       canonical = TRUE,
#       cap_first = FALSE,
#       best_match_only = FALSE,
#       http = "post",
#       fields = "all"
#     )
#   
#   if (nrow(df8_dis) <= 1000 &
#       nrow(df8_dis) > 500)
#     df9_b <- gnr_resolve(
#       names = df8_dis_b$query,
#       data_source_ids = NULL,
#       resolve_once = FALSE,
#       with_context = TRUE,
#       canonical = TRUE,
#       cap_first = FALSE,
#       best_match_only = FALSE,
#       http = "post",
#       fields = "all"
#     )
#   if (nrow(df8_dis) <= 1000 &
#       nrow(df8_dis) > 500)
#     df9_b[setdiff(x = names(df4), y = names(df9_b))] <- NA
#   
#   if (nrow(df8_dis) <= 1000 &
#       nrow(df8_dis) > 500)
#     df9 <- rbind(df9_a, df9_b)
#   
#   #selecting best result (with best score and best filled taxonomy)
#   df10 <- df9 %>%
#     rowwise() %>%
#     mutate(
#       kingdom = sum(as.numeric(
#         stri_detect(
#           str = classification_path_ranks,
#           fixed = c("kingdom", "Kingdom", "regn.")
#         )
#       )),
#       phylum =  sum(as.numeric(
#         stri_detect(str = classification_path_ranks,
#                     fixed = c("phylum", "Phylum", "phyl."))
#       )),
#       class =  sum(as.numeric(
#         stri_detect(str = classification_path_ranks,
#                     fixed = c("class", "Class", "cl."))
#       )),
#       order =  sum(as.numeric(
#         stri_detect(str = classification_path_ranks,
#                     fixed = c("order", "Order", "ord."))
#       )),
#       family = sum(as.numeric(
#         stri_detect(str = classification_path_ranks,
#                     fixed = c("family", "Family", "fam."))
#       )),
#       genus =  sum(as.numeric(
#         stri_detect(str = classification_path_ranks,
#                     fixed = c("genus", "Genus"))
#       )),
#       species = sum(as.numeric(
#         stri_detect(
#           str = classification_path_ranks,
#           fixed = c("species", "Species", "spec.", "sp.")
#         )
#       )),
#       variety = sum(as.numeric(
#         stri_detect(
#           str = classification_path_ranks,
#           fixed = c("variety",
#                     "varietas",
#                     "var")
#         )
#       ))
#     ) %>%
#     ungroup()
#   
#   df10$kingdom[df10$kingdom >= 1] <- 1
#   df10$phylum[df10$phylum >= 1] <- 1
#   df10$class[df10$class >= 1] <- 1
#   df10$order[df10$order >= 1] <- 1
#   df10$family[df10$family >= 1] <- 1
#   df10$genus[df10$genus >= 1] <- 1
#   df10$species[df10$species >= 1] <- 1
#   df10$variety[df10$variety >= 1] <- 1
#   
#   df10 <- df10 %>%
#     mutate(n = rowSums(.[c("kingdom",
#                            "phylum",
#                            "class",
#                            "order",
#                            "family",
#                            "genus",
#                            "species",
#                            "variety")])) %>%
#     group_by(user_supplied_name) %>%
#     arrange(desc(score), desc(n)) %>%
#     ungroup() %>%
#     distinct(user_supplied_name,
#              .keep_all = TRUE) %>%
#     select(
#       user_supplied_name,
#       canonicalname = matched_name2,
#       taxonId = taxon_id,
#       dbTaxo = data_source_title,
#       taxonomy = classification_path,
#       rank = classification_path_ranks
#     )
#   
#   #manipulating taxa
#   df11 <- manipulating_taxo(dfsel = df10,
#                             dic = taxaRanksDictionary)
#   
#   #joining
#   df12 <-
#     left_join(df8,
#               df11,
#               by = c("query" = "user_supplied_name"))
#   
#   #harmonizing columns (for rbind later on)
#   df12[setdiff(x = names(df7), y = names(df12))] <- NA
#   
#   #filtering non-empty taxonomies
#   df12_full <- df12 %>%
#     filter(
#       !is.na(kingdom) |
#         !is.na(phylum) |
#         !is.na(class) |
#         !is.na(order) |
#         !is.na(family) |
#         !is.na(genus)
#     ) %>%
#     select(-rank, -taxonomy, -query)
#   
#   #filtering empty taxonomies and removing non-UTF8 characters (else cause bugs when running gnr_resolve)
#   df13 <- df12 %>%
#     filter(
#       is.na(kingdom) &
#         is.na(phylum) &
#         is.na(class) &
#         is.na(order) &
#         is.na(family) &
#         is.na(genus)
#     ) %>%
#     select(organismTranslated,
#            query) %>%
#     mutate_at(
#       .vars = vars(organismTranslated, query),
#       .funs = function(x) {
#         gsub(pattern = "[^ -~]",
#              replacement =  "",
#              x =  x)
#       }
#     ) %>%
#     mutate_at(
#       .vars = vars(organismTranslated, query),
#       .funs = function(x) {
#         gsub(pattern = ";",
#              replacement = "",
#              x = x)
#       }
#     ) %>%
#     mutate(query = word(string = query, start = 1))
#   
#   
#   #outputting distinct empty taxonomies to avoid iterations
#   df13_dis <- df13 %>%
#     distinct(query)
#   
#   if (nrow(df13_dis) == 0)
#     df13_dis[nrow(df13_dis) + 1, 1] <- "NA na"
#   
#   if (nrow(df13_dis) == 0)
#     df13_dis[nrow(df13_dis) + 1, 1] <- "NA na"
#   
#   if (nrow(df13_dis) == 0)
#     df13_dis[nrow(df13_dis) + 1, 2] <- "NA na"
#   
#   if (nrow(df13_dis) == 0)
#     df13_dis[nrow(df13_dis) + 1, 2] <- "NA na"
#   
#   #running gnresolve on first string of identifier cleaned
#   ##if less than 500 empty entries
#   if (nrow(df13_dis) <= 500)
#     df14 <- gnr_resolve(
#       names = df13_dis$query,
#       data_source_ids = NULL,
#       resolve_once = FALSE,
#       with_context = TRUE,
#       canonical = TRUE,
#       cap_first = FALSE,
#       best_match_only = FALSE,
#       http = "post",
#       fields = "all"
#     )
#   
#   ##if more than 1000 empty entries
#   if (nrow(df13_dis) > 1000)
#     df13_dis_a <- slice(df13_dis, 1:500)
#   
#   if (nrow(df13_dis) > 1000)
#     df13_dis_b <- slice(df13_dis, 501:1000)
#   
#   if (nrow(df13_dis) > 1000)
#     df13_dis_c <- slice(df13_dis, 1001:nrow(df13_dis))
#   
#   if (nrow(df13_dis) > 1000)
#     df14_a <- gnr_resolve(
#       names = df13_dis_a$query,
#       data_source_ids = NULL,
#       resolve_once = FALSE,
#       with_context = TRUE,
#       canonical = TRUE,
#       cap_first = FALSE,
#       best_match_only = FALSE,
#       http = "post",
#       fields = "all"
#     )
#   
#   if (nrow(df13_dis) > 1000)
#     df14_b <- gnr_resolve(
#       names = df13_dis_b$query,
#       data_source_ids = NULL,
#       resolve_once = FALSE,
#       with_context = TRUE,
#       canonical = TRUE,
#       cap_first = FALSE,
#       best_match_only = FALSE,
#       http = "post",
#       fields = "all"
#     )
#   
#   if (nrow(df13_dis) > 1000)
#     df14_c <- gnr_resolve(
#       names = df13_dis_c$query,
#       data_source_ids = NULL,
#       resolve_once = FALSE,
#       with_context = TRUE,
#       canonical = TRUE,
#       cap_first = FALSE,
#       best_match_only = FALSE,
#       http = "post",
#       fields = "all"
#     )
#   
#   if (nrow(df13_dis) > 1000)
#     df14_a[setdiff(x = names(df4), y = names(df14_a))] <- NA
#   
#   if (nrow(df13_dis) > 1000)
#     df14_b[setdiff(x = names(df14_a), y = names(df14_b))] <- NA
#   
#   if (nrow(df13_dis) > 1000)
#     df14_c[setdiff(x = names(df14_b), y = names(df14_c))] <- NA
#   
#   if (nrow(df13_dis) > 1000)
#     df14 <- rbind(df14_a, df14_b, df14_c)
#   
#   ##if between 500 and 1000 empty entries
#   if (nrow(df13_dis) <= 1000 &
#       nrow(df13_dis) > 500)
#     df13_dis_a <- slice(df13_dis, 1:500)
#   
#   if (nrow(df13_dis) <= 1000 &
#       nrow(df13_dis) > 500)
#     df13_dis_b <- slice(df13_dis, 501:nrow(df13_dis))
#   
#   if (nrow(df13_dis) <= 1000 &
#       nrow(df13_dis) > 500)
#     df14_a <- gnr_resolve(
#       names = df13_dis_a$query,
#       data_source_ids = NULL,
#       resolve_once = FALSE,
#       with_context = TRUE,
#       canonical = TRUE,
#       cap_first = FALSE,
#       best_match_only = FALSE,
#       http = "post",
#       fields = "all"
#     )
#   if (nrow(df13_dis) <= 1000 &
#       nrow(df13_dis) > 500)
#     df14_b <- gnr_resolve(
#       names = df13_dis_b$query,
#       data_source_ids = NULL,
#       resolve_once = FALSE,
#       with_context = TRUE,
#       canonical = TRUE,
#       cap_first = FALSE,
#       best_match_only = FALSE,
#       http = "post",
#       fields = "all"
#     )
#   
#   if (nrow(df13_dis) <= 1000 &
#       nrow(df13_dis) > 500)
#     df14_b[setdiff(x = names(df14_a), y = names(df14_b))] <- NA
#   
#   if (nrow(df13_dis) <= 1000 &
#       nrow(df13_dis) > 500)
#     df14 <- rbind(df14_a, df14_b)
#   
#   #selecting best result (with best score and best filled taxonomy)
#   df15 <- df14 %>%
#     rowwise() %>%
#     mutate(
#       kingdom = sum(as.numeric(
#         stri_detect(
#           str = classification_path_ranks,
#           fixed = c("kingdom", "Kingdom", "regn.")
#         )
#       )),
#       phylum =  sum(as.numeric(
#         stri_detect(str = classification_path_ranks,
#                     fixed = c("phylum", "Phylum", "phyl."))
#       )),
#       class =  sum(as.numeric(
#         stri_detect(str = classification_path_ranks,
#                     fixed = c("class", "Class", "cl."))
#       )),
#       order =  sum(as.numeric(
#         stri_detect(str = classification_path_ranks,
#                     fixed = c("order", "Order", "ord."))
#       )),
#       family = sum(as.numeric(
#         stri_detect(str = classification_path_ranks,
#                     fixed = c("family", "Family", "fam."))
#       )),
#       genus =  sum(as.numeric(
#         stri_detect(str = classification_path_ranks,
#                     fixed = c("genus", "Genus"))
#       )),
#       species = sum(as.numeric(
#         stri_detect(
#           str = classification_path_ranks,
#           fixed = c("species", "Species", "spec.", "sp.")
#         )
#       )),
#       variety = sum(as.numeric(
#         stri_detect(
#           str = classification_path_ranks,
#           fixed = c("variety",
#                     "varietas",
#                     "var")
#         )
#       ))
#     ) %>%
#     ungroup()
#   
#   df15$kingdom[df15$kingdom >= 1] <- 1
#   df15$phylum[df15$phylum >= 1] <- 1
#   df15$class[df15$class >= 1] <- 1
#   df15$order[df15$order >= 1] <- 1
#   df15$family[df15$family >= 1] <- 1
#   df15$genus[df15$genus >= 1] <- 1
#   df15$species[df15$species >= 1] <- 1
#   df15$variety[df15$variety >= 1] <- 1
#   
#   df15 <- df15 %>%
#     mutate(n = rowSums(.[c("kingdom",
#                            "phylum",
#                            "class",
#                            "order",
#                            "family",
#                            "genus",
#                            "species",
#                            "variety")])) %>%
#     group_by(user_supplied_name) %>%
#     arrange(desc(score), desc(n)) %>%
#     ungroup() %>%
#     distinct(user_supplied_name,
#              .keep_all = TRUE) %>%
#     select(
#       user_supplied_name,
#       canonicalname = matched_name2,
#       taxonId = taxon_id,
#       dbTaxo = data_source_title,
#       taxonomy = classification_path,
#       rank = classification_path_ranks
#     )
#   
#   #manipulating taxa
#   df16 <- manipulating_taxo(dfsel = df15,
#                             dic = taxaRanksDictionary)
#   
#   #joining
#   df17 <-
#     left_join(df13,
#               df16,
#               by = c("query" = "user_supplied_name"))
#   
#   #harmonizing columns (for rbind later on)
#   df17[setdiff(x = names(df12), y = names(df17))] <- NA
#   
#   #filtering non-empty taxonomies
#   df17_full <- df17 %>%
#     filter(
#       !is.na(kingdom) |
#         !is.na(phylum) |
#         !is.na(class) |
#         !is.na(order) |
#         !is.na(family) |
#         !is.na(genus)
#     ) %>%
#     select(-rank, -taxonomy, -query)
#   
#   #filtering empty taxonomies
#   df17_empty <- df17 %>%
#     filter(
#       is.na(kingdom) &
#         is.na(phylum) &
#         is.na(class) &
#         is.na(order) &
#         is.na(family) &
#         is.na(genus)
#     ) %>%
#     select(-rank, -taxonomy, -query)
#   
#   #selecting joining variable
#   y <- x %>%
#     select(organismTranslated)
#   
#   #joining (part 1)
#   prefinal_df <- left_join(y, df2_full)
#   
#   prefinal_df_2 <- rbind(df7_full,
#                          df12_full,
#                          df17_full,
#                          df17_empty)
#   
#   dbQuality <- df2_full %>%
#     select(dbTaxo, dbQuality) %>% 
#     distinct(dbTaxo, .keep_all = TRUE)
#   
#   prefinal_df_3 <- left_join(prefinal_df_2, dbQuality)
#   
#   #joining (rest) and removing NA when multiple organisms found
#   final_df <- rbind(prefinal_df,
#                     prefinal_df_3) %>%
#     distinct(organismTranslated,
#              canonicalname,
#              .keep_all = TRUE) %>%
#     group_by(organismTranslated) %>%
#     add_count() %>%
#     ungroup() %>%
#     filter(!is.na(canonicalname) |
#              !n > 1) %>%
#     select(-n) %>%
#     distinct(organismTranslated,
#              canonicalname,
#              .keep_all = TRUE)
#   
#   return(final_df)
# }

#######################################################
#######################################################

y_as_na <- function(x, y)
{
  if ("factor" %in% class(x))
    x <- as.character(x) ## since ifelse wont work with factors
  ifelse(test = as.character(x) != y,
         yes = x,
         no = NA)
}

#######################################################
#######################################################

taxo_cleaning_manual <- function(dfsel)
{
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
        organism_7_species
      ),
      .funs = function(x) {
        gsub(pattern = "^c\\(|\\)$",
             replacement =  "",
             x =  x)
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
        organism_7_species
      ),
      .funs = function(x) {
        gsub(pattern = "\"",
             replacement = "",
             x = x)
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
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Allomyrina dichotoma"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Allomyrina dichotoma"] <-
    "Animalia"
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Allomyrina dichotoma"] <-
    "Arthropoda"
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Allomyrina dichotoma"] <-
    "Insecta"
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Allomyrina dichotoma"] <-
    "Coleoptera"
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Allomyrina dichotoma"] <-
    "Scarabaeidae"
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Allomyrina dichotoma"] <-
    "Trypoxylus"
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Allomyrina dichotoma"] <-
    "Trypoxylus dichotomus"
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Agaricus pattersonae"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Agaricus pattersonae"] <-
    "Fungi"
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Agaricus pattersonae"] <-
    "Basidiomycota"
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Agaricus pattersonae"] <-
    "Agaricomycetes"
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Agaricus pattersonae"] <-
    "Agaricales"
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Agaricus pattersonae"] <-
    "Agaricaceae"
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Agaricus pattersonae"] <-
    "Agaricus"
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Agaricus pattersonae"] <-
    "Agaricus pattersoniae"
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Melaphis chinensis"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Melaphis chinensis"] <-
    "Animalia"
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Melaphis chinensis"] <-
    "Arthropoda"
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Melaphis chinensis"] <-
    "Insecta"
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Melaphis chinensis"] <-
    "Hemiptera"
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Melaphis chinensis"] <-
    "Aphididae"
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Melaphis chinensis"] <-
    "Schlechtendalia"
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Melaphis chinensis"] <-
    "Schlechtendalia chinensis"
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Chloroclysta truncata"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Chloroclysta truncata"] <-
    "Animalia"
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Chloroclysta truncata"] <-
    "Arthropoda"
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Chloroclysta truncata"] <-
    "Insecta"
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Chloroclysta truncata"] <-
    "Lepidoptera"
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Chloroclysta truncata"] <-
    "Geometridae"
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Chloroclysta truncata"] <-
    "Dysstroma"
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Chloroclysta truncata"] <-
    "Dysstroma truncata"
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Lindenbergia urticaefolia"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Lindenbergia urticaefolia"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Lindenbergia urticaefolia"] <-
    "Tracheophyta"
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Lindenbergia urticaefolia"] <-
    "Magnoliopsida"
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Lindenbergia urticaefolia"] <-
    "Lamiales"
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Lindenbergia urticaefolia"] <-
    "Orobanchaceae"
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Lindenbergia urticaefolia"] <-
    "Lindenbergia"
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Lindenbergia urticaefolia"] <-
    "Lindenbergia urticifolia"
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Tetraselmis chui"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Tetraselmis chui"] <-
    "Bacteria"
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Tetraselmis chui"] <-
    "Chlorophyta"
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Tetraselmis chui"] <-
    "Chlorodendrophyceae"
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Tetraselmis chui"] <-
    "Chlorodendrales"
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Tetraselmis chui"] <-
    "Chlorodendraceae"
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Tetraselmis chui"] <-
    "Tetraselmis"
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Tetraselmis chui"] <-
    "Tetraselmis chuii"
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Cyanospira rippkae"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Cyanospira rippkae"] <-
    "Bacteria"
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Cyanospira rippkae"] <-
    "Cyanobacteria"
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Cyanospira rippkae"] <-
    "Cyanophyceae"
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Cyanospira rippkae"] <-
    "Nostocales"
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Cyanospira rippkae"] <-
    "Aphanizomenonaceae"
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Cyanospira rippkae"] <-
    "Cyanospira"
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Cyanospira rippkae"] <-
    "Cyanospira rippkae"
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Nicandra physaloides"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Nicandra physaloides"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Nicandra physaloides"] <-
    "Tracheophyta"
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Nicandra physaloides"] <-
    "Magnoliopsida"
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Nicandra physaloides"] <-
    "Solanales"
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Nicandra physaloides"] <-
    "Solanaceae"
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Nicandra physaloides"] <-
    "Nicandra"
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Nicandra physaloides"] <-
    "Nicandra physalodes"
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Salvia shannoni"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Salvia shannoni"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Salvia shannoni"] <-
    "Tracheophyta"
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Salvia shannoni"] <-
    "Magnoliopsida"
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Salvia shannoni"] <-
    "Lamiales"
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Salvia shannoni"] <-
    "Lamiaceae"
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Salvia shannoni"] <-
    "Salvia"
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Salvia shannoni"] <-
    ""
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Streptomyces tsukubaensis"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Streptomyces tsukubaensis"] <-
    "Bacteria"
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Streptomyces tsukubaensis"] <-
    "Actinobacteria"
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Streptomyces tsukubaensis"] <-
    "Actinobacteria"
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Streptomyces tsukubaensis"] <-
    "Actinomycetales"
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Streptomyces tsukubaensis"] <-
    "Streptomycetaceae"
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Streptomyces tsukubaensis"] <-
    "Streptomyces"
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Streptomyces tsukubaensis"] <-
    "Streptomyces tsukubensis"
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Evea brasiliensis"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Evea brasiliensis"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Evea brasiliensis"] <-
    "Tracheophyta"
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Evea brasiliensis"] <-
    "Magnoliopsida"
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Evea brasiliensis"] <-
    "Malpighiales"
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Evea brasiliensis"] <-
    "Euphorbiaceae"
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Evea brasiliensis"] <-
    "Hevea"
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Evea brasiliensis"] <-
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
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Chinensis"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Chinensis"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Chinensis"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Chinensis"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Chinensis"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Chinensis"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Chinensis"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Chinensis"] <-
    ""
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Sinensis"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Sinensis"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Sinensis"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Sinensis"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Sinensis"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Sinensis"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Sinensis"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Sinensis"] <-
    ""
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Ootheca"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Ootheca"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Ootheca"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Ootheca"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Ootheca"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Ootheca"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Ootheca"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Ootheca"] <-
    ""
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Uncultured"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Uncultured"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Uncultured"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Uncultured"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Uncultured"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Uncultured"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Uncultured"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Uncultured"] <-
    ""
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Stigma"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Stigma"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Stigma"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Stigma"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Stigma"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Stigma"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Stigma"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Stigma"] <-
    ""
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Spica"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Spica"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Spica"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Spica"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Spica"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Spica"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Spica"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Spica"] <-
    ""
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Semen"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Semen"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Semen"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Semen"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Semen"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Semen"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Semen"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Semen"] <-
    ""
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Rotundus"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Rotundus"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Rotundus"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Rotundus"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Rotundus"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Rotundus"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Rotundus"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Rotundus"] <-
    ""
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Rhizoma"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Rhizoma"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Rhizoma"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Rhizoma"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Rhizoma"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Rhizoma"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Rhizoma"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Rhizoma"] <-
    ""
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Ramulus"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Ramulus"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Ramulus"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Ramulus"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Ramulus"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Ramulus"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Ramulus"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Ramulus"] <-
    ""
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Radix"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Radix"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Radix"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Radix"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Radix"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Radix"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Radix"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Radix"] <-
    ""
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Pollen"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Pollen"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Pollen"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Pollen"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Pollen"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Pollen"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Pollen"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Pollen"] <-
    ""
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Lignum"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Lignum"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Lignum"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Lignum"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Lignum"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Lignum"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Lignum"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Lignum"] <-
    ""
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Fructus"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Fructus"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Fructus"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Fructus"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Fructus"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Fructus"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Fructus"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Fructus"] <-
    ""
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Flos"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Flos"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Flos"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Flos"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Flos"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Flos"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Flos"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Flos"] <-
    ""
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Corolla"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Corolla"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Corolla"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Corolla"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Corolla"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Corolla"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Corolla"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Corolla"] <-
    ""
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Cacumen"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Cacumen"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Cacumen"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Cacumen"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Cacumen"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Cacumen"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Cacumen"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Cacumen"] <-
    ""
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Bulbus"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Bulbus"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Bulbus"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Bulbus"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Bulbus"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Bulbus"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Bulbus"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Bulbus"] <-
    ""
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Megaleia"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Megaleia"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Megaleia"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Megaleia"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Megaleia"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Megaleia"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Megaleia"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Megaleia"] <-
    ""
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Candidatus"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Candidatus"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Candidatus"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Candidatus"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Candidatus"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Candidatus"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Candidatus"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Candidatus"] <-
    ""
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Tasmanian"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Tasmanian"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Tasmanian"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Tasmanian"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Tasmanian"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Tasmanian"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Tasmanian"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Tasmanian"] <-
    ""
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Asian"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Asian"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Asian"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Asian"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Asian"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Asian"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Asian"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Asian"] <-
    ""
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Mammalian"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Mammalian"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Mammalian"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Mammalian"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Mammalian"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Mammalian"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Mammalian"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Mammalian"] <-
    ""
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Red"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Red"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Red"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Red"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Red"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Red"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Red"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Red"] <-
    ""
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Turkey"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Turkey"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Turkey"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Turkey"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Turkey"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Turkey"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Turkey"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Turkey"] <-
    ""
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Synthetis"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Synthetis"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Synthetis"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Synthetis"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Synthetis"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Synthetis"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Synthetis"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Synthetis"] <-
    ""
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Pagellus erythrinus"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Pagellus erythrinus"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Pagellus erythrinus"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Pagellus erythrinus"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Pagellus erythrinus"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Pagellus erythrinus"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Pagellus erythrinus"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Pagellus erythrinus"] <-
    ""
  #comes from Becker translation
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Lagenorhynchus obliquidens"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Lagenorhynchus obliquidens"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Lagenorhynchus obliquidens"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Lagenorhynchus obliquidens"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Lagenorhynchus obliquidens"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Lagenorhynchus obliquidens"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Lagenorhynchus obliquidens"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Lagenorhynchus obliquidens"] <-
    ""
  #comes from Lag translation
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Actinomycetales bacterium	"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Actinomycetales bacterium	"] <-
    "Bacteria"
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Actinomycetales bacterium	"] <-
    "Actinobacteria"
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Actinomycetales bacterium	"] <-
    "Actinobacteria"
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Actinomycetales bacterium	"] <-
    "Actinomycetales"
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Actinomycetales bacterium	"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Actinomycetales bacterium	"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Actinomycetales bacterium	"] <-
    ""
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Galla"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Galla"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Galla"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Galla"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Galla"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Galla"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Galla"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Galla"] <-
    ""
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Red Sea bacterium KT-2K1"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Red Sea bacterium KT-2K1"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Red Sea bacterium KT-2K1"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Red Sea bacterium KT-2K1"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Red Sea bacterium KT-2K1"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Red Sea bacterium KT-2K1"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Red Sea bacterium KT-2K1"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Red Sea bacterium KT-2K1"] <-
    ""
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Green Pelican GFP transformation vector"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Green Pelican GFP transformation vector"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Green Pelican GFP transformation vector"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Green Pelican GFP transformation vector"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Green Pelican GFP transformation vector"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Green Pelican GFP transformation vector"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Green Pelican GFP transformation vector"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Green Pelican GFP transformation vector"] <-
    ""
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Peripatoides"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Peripatoides"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Peripatoides"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Peripatoides"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Peripatoides"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Peripatoides"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Peripatoides"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Peripatoides"] <-
    ""
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Bacterium MPBA1"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Bacterium MPBA1"] <-
    "Bacteria"
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Bacterium MPBA1"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Bacterium MPBA1"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Bacterium MPBA1"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Bacterium MPBA1"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Bacterium MPBA1"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Bacterium MPBA1"] <-
    ""
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Rhizobiaceae bacterium"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Rhizobiaceae bacterium"] <-
    "Bacteria"
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Rhizobiaceae bacterium"] <-
    "Proteobacteria"
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Rhizobiaceae bacterium"] <-
    "Alphaproteobacteria"
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Rhizobiaceae bacterium"] <-
    "Rhizobiales"
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Rhizobiaceae bacterium"] <-
    "Rhizobiaceae"
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Rhizobiaceae bacterium"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Rhizobiaceae bacterium"] <-
    ""
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Cyanophyta"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Cyanophyta"] <-
    "Bacteria"
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Cyanophyta"] <-
    "Cyanobacteria"
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Cyanophyta"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Cyanophyta"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Cyanophyta"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Cyanophyta"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Cyanophyta"] <-
    ""
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Candida"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Candida"] <-
    "Fungi"
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Candida"] <-
    "Ascomycota"
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Candida"] <-
    "Saccharomycetes"
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Candida"] <-
    "Saccharomycetales"
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Candida"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Candida"] <-
    "Candida"
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Candida"] <-
    ""
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Microsorium"] <-
    "y"
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Microsorium"] <-
    "Microsorum"
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Pseudoeurotium"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Pseudoeurotium"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Pseudoeurotium"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Pseudoeurotium"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Pseudoeurotium"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Pseudoeurotium"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Pseudoeurotium"] <-
    "Pseudeurotium"
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Pseudoeurotium"] <-
    ""
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Lanea"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Lanea"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Lanea"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Lanea"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Lanea"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Lanea"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Lanea"] <-
    "Lannea"
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Lanea"] <-
    ""
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Aspidium"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Aspidium"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Aspidium"] <-
    "Tracheophyta"
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Aspidium"] <-
    "Polypodiopsida"
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Aspidium"] <-
    "Polypodiales"
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Aspidium"] <-
    "Dryopteridaceae"
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Aspidium"] <-
    "Polystichum"
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Aspidium"] <-
    ""
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Iris"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Iris"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Iris"] <-
    "Tracheophyta"
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Iris"] <-
    "Liliopsida"
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Iris"] <-
    "Asparagales"
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Iris"] <-
    "Iridaceae"
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Iris"] <-
    "Iris"
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Iris"] <-
    ""
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Plantae"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Plantae"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Plantae"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Plantae"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Plantae"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Plantae"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Plantae"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Plantae"] <-
    ""
  
  #sadly can be multiple kingdoms, therefore droped
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Algae"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Algae"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Algae"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Algae"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Algae"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Algae"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Algae"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Algae"] <-
    ""
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Fungi"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Fungi"] <-
    "Fungi"
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Fungi"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Fungi"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Fungi"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Fungi"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Fungi"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Fungi"] <-
    ""
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Anaerobic"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Anaerobic"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Anaerobic"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Anaerobic"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Anaerobic"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Anaerobic"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Anaerobic"] <-
    ""
  inhouse_db$organism_7_species[inhouse_db$organismSanitized == "Anaerobic"] <-
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
  
  #examples mismatched genera
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
  
  #example
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_7_species == "Solanum etuberosum"] <-
    "y"
  inhouse_db$organism_7_species[inhouse_db$organism_7_species == "Solanum etuberosum"] <-
    "Solanum tuberosum"
  
  #example
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
  
  #mismatched genus
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
  
  #example
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
  
  #example_2
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
  
  #double taxonomies -> no "y", just choose one
  #catalogue of Life as reference, then NCBI
  
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
  #no organism_4_order
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
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organismSanitized == "Mayodendron"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organismSanitized == "Mayodendron"] <-
    ""
  inhouse_db$organism_2_phylum[inhouse_db$organismSanitized == "Mayodendron"] <-
    ""
  inhouse_db$organism_3_class[inhouse_db$organismSanitized == "Mayodendron"] <-
    ""
  inhouse_db$organism_4_order[inhouse_db$organismSanitized == "Mayodendron"] <-
    ""
  inhouse_db$organism_5_family[inhouse_db$organismSanitized == "Mayodendron"] <-
    ""
  inhouse_db$organism_6_genus[inhouse_db$organismSanitized == "Mayodendron"] <-
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
  #no order
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
  
  #double taxonomies -> no "y", just choose one
  #catalogue of Life as reference, then NCBI
  
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
  
  #double taxonomies -> no "y", just choose one
  #catalogue of Life as reference, then NCBI
  
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
  
  #debate
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
  
  #not even in the initial two
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
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_4_order == "Ranunculales"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_4_order == "Ranunculales"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organism_4_order == "Ranunculales"] <-
    "Tracheophyta"
  inhouse_db$organism_3_class[inhouse_db$organism_4_order == "Ranunculales"] <-
    "Magnoliopsida"
  
  inhouse_db$organism_modified_taxonomy_manual[inhouse_db$organism_4_order == "Brassicales"] <-
    "y"
  inhouse_db$organism_1_kingdom[inhouse_db$organism_4_order == "Brassicales"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organism_4_order == "Brassicales"] <-
    "Tracheophyta"
  inhouse_db$organism_3_class[inhouse_db$organism_4_order == "Brassicales"] <-
    "Magnoliopsida"
  
  inhouse_db$organism_1_kingdom[inhouse_db$organism_4_order == "Asparagales"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organism_4_order == "Asparagales"] <-
    "Tracheophyta"
  inhouse_db$organism_3_class[inhouse_db$organism_4_order == "Asparagales"] <-
    "Magnoliopsida"
  
  inhouse_db$organism_1_kingdom[inhouse_db$organism_4_order == "Alismatales"] <-
    "Plantae"
  inhouse_db$organism_2_phylum[inhouse_db$organism_4_order == "Alismatales"] <-
    "Tracheophyta"
  inhouse_db$organism_3_class[inhouse_db$organism_4_order == "Alismatales"] <-
    "Magnoliopsida"
  
  #double taxonomies -> no "y", just choose one
  #catalogue of Life as reference, then NCBI
  
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
  
  #double taxonomies -> no "y", just choose one
  #catalogue of Life as reference, then NCBI
  
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
    y_as_na(x = inhouse_db$organism_1_kingdom,
            y = "")
  
  inhouse_db$organism_2_phylum <-
    y_as_na(x = inhouse_db$organism_2_phylum,
            y = "")
  
  inhouse_db$organism_3_class <-
    y_as_na(x = inhouse_db$organism_3_class,
            y = "")
  
  inhouse_db$organism_4_order <-
    y_as_na(x = inhouse_db$organism_4_order,
            y = "")
  
  inhouse_db$organism_5_family <-
    y_as_na(x = inhouse_db$organism_5_family,
            y = "")
  
  inhouse_db$organism_6_genus <-
    y_as_na(x = inhouse_db$organism_6_genus,
            y = "")
  
  inhouse_db$organism_7_species <-
    y_as_na(x = inhouse_db$organism_7_species,
            y = "")
  
  inhouse_db$organism_modified_taxonomy_manual <-
    y_as_na(x = inhouse_db$organism_modified_taxonomy_manual,
            y = "")
  
  inhouse_db$organismCurated <-
    as.character(apply(inhouse_db[5:11], 1, function(x)
      tail(na.omit(x), 1)))
  
  organism_7_species_cleaning <-
    inhouse_db %>%
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
    select(organismOriginal,
           organismTranslated,
           organismSanitized,
           organismCurated)
  
  if (nrow(inhouse_db_new) > 0)
    inhouse_db_new <-
    left_join(inhouse_db_new,
              organism_7_species_cleaning)
  
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
    select(organismOriginal,
           organismTranslated,
           organismSanitized,
           organismCurated)
  
  if (nrow(inhouse_db_new) > 0)
    inhouse_db_new <-
    left_join(inhouse_db_new,
              organism_6_genus_cleaning)
  
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
    select(organismOriginal,
           organismTranslated,
           organismSanitized,
           organismCurated)
  
  if (nrow(inhouse_db_new) > 0)
    inhouse_db_new <-
    left_join(inhouse_db_new,
              organism_5_family_cleaning)
  
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
    select(organismOriginal,
           organismTranslated,
           organismSanitized,
           organismCurated)
  
  if (nrow(inhouse_db_new) > 0)
    inhouse_db_new <-
    left_join(inhouse_db_new,
              organism_4_order_cleaning)
  
  inhouse_db_organism_4_order_clean <-
    rbind(inhouse_db_old, inhouse_db_new)
  
  organism_3_class_cleaning <-
    inhouse_db_organism_4_order_clean %>%
    distinct(organism_1_kingdom,
             organism_2_phylum,
             organism_3_class,
             .keep_all = TRUE) %>%
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
    select(organismOriginal,
           organismTranslated,
           organismSanitized,
           organismCurated)
  
  if (nrow(inhouse_db_new) > 0)
    inhouse_db_new <-
    left_join(inhouse_db_new, organism_3_class_cleaning)
  
  inhouse_db_organism_3_class_clean <-
    rbind(inhouse_db_old, inhouse_db_new)
  
  organism_2_phylum_cleaning <-
    inhouse_db_organism_3_class_clean %>%
    distinct(organism_1_kingdom,
             organism_2_phylum,
             .keep_all = TRUE) %>%
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
    filter(organismCurated %in% organism_2_phylum_cleaning$organismCurated)  %>%
    select(organismOriginal,
           organismTranslated,
           organismSanitized,
           organismCurated)
  
  if (nrow(inhouse_db_new) > 0)
    inhouse_db_new <-
    left_join(inhouse_db_new,
              organism_2_phylum_cleaning)
  
  inhouse_db_organism_2_phylum_clean <-
    rbind(inhouse_db_old, inhouse_db_new)
  
  organism_1_kingdom_cleaning <-
    inhouse_db_organism_2_phylum_clean %>%
    distinct(organism_1_kingdom,
             .keep_all = TRUE)
  
  organism_1_kingdom_cleaning_2 <-
    inhouse_db_organism_2_phylum_clean %>%
    filter(!is.na(organism_1_kingdom)) %>%
    distinct(organism_1_kingdom,
             .keep_all = TRUE)
  
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
  
  if (nrow(inhouse_db_new) > 0)
    inhouse_db_new <-
    left_join(inhouse_db_new, organism_1_kingdom_cleaning)
  
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
      x = inhouse_db_organism_1_kingdom_clean$organismSanitized
    )
  
  inhouse_db_organism_1_kingdom_clean$organismCurated <-
    gsub(
      pattern = "Fructus",
      replacement = "",
      x = inhouse_db_organism_1_kingdom_clean$organismSanitized
    )
  
  inhouse_db_organism_1_kingdom_clean$organismCurated <-
    gsub(
      pattern = "Radix",
      replacement = "",
      x = inhouse_db_organism_1_kingdom_clean$organismSanitized
    )
  
  #to avoid false genera
  inhouse_db_organism_1_kingdom_clean$organism_7_species[str_count(string = inhouse_db_organism_1_kingdom_clean$organism_7_species,
                                                                   pattern = "\\w+") == 1] <-
    ""
  
  inhouse_db_organism_1_kingdom_clean$organism_1_kingdom <-
    y_as_na(x = inhouse_db_organism_1_kingdom_clean$organism_1_kingdom,
            y = "")
  
  inhouse_db_organism_1_kingdom_clean$organism_2_phylum <-
    y_as_na(x = inhouse_db_organism_1_kingdom_clean$organism_2_phylum,
            y = "")
  
  inhouse_db_organism_1_kingdom_clean$organism_3_class <-
    y_as_na(x = inhouse_db_organism_1_kingdom_clean$organism_3_class,
            y = "")
  
  inhouse_db_organism_1_kingdom_clean$organism_4_order <-
    y_as_na(x = inhouse_db_organism_1_kingdom_clean$organism_4_order,
            y = "")
  
  inhouse_db_organism_1_kingdom_clean$organism_5_family <-
    y_as_na(x = inhouse_db_organism_1_kingdom_clean$organism_5_family,
            y = "")
  
  inhouse_db_organism_1_kingdom_clean$organism_6_genus <-
    y_as_na(x = inhouse_db_organism_1_kingdom_clean$organism_6_genus,
            y = "")
  
  inhouse_db_organism_1_kingdom_clean$organism_7_species <-
    y_as_na(x = inhouse_db_organism_1_kingdom_clean$organism_7_species,
            y = "")
  
  inhouse_db_organism_1_kingdom_clean$organismCurated <-
    y_as_na(x = inhouse_db_organism_1_kingdom_clean$organismCurated,
            y = "")
  
  inhouse_db_organism_1_kingdom_clean$organism_modified_taxonomy_manual <-
    y_as_na(x = inhouse_db_organism_1_kingdom_clean$organism_modified_taxonomy_manual,
            y = "")
  
  inhouse_db_organism_1_kingdom_clean$organismCurated <-
    as.character(apply(inhouse_db_organism_1_kingdom_clean[5:11], 1, function(x)
      tail(na.omit(x), 1)))
  
  inhouse_db_organism_1_kingdom_clean <-
    inhouse_db_organism_1_kingdom_clean %>%
    mutate_all(as.character)
  
  inhouse_db_organism_1_kingdom_clean$organismCurated <-
    y_as_na(x = inhouse_db_organism_1_kingdom_clean$organismCurated,
            y = "character(0)")
  
  inhouse_db_organism_1_kingdom_clean$organismCurated <-
    y_as_na(x = inhouse_db_organism_1_kingdom_clean$organismCurated,
            y = "NA")
  
  return(inhouse_db_organism_1_kingdom_clean)
}

#######################################################
#######################################################

plantcycompiling <- function(x)
{
  path <- paste(x,
                "/",
                file,
                sep = "")
  
  data <- read.delim(file = path)
  
  data$X. <- as.character(data$X.)
  
  data$X. <- sub(" - ", " ", data$X.)
  
  data_2 <- data %>%
    cSplit(splitCols = "X.",
           sep = " ") %>%
    select(1:2) %>% 
    tibble()
  
  colnames(data_2)[1] <- "A"
  colnames(data_2)[2] <- "B"
  
  data_2$A <- as.character(data_2$A)
  
  data_3 <- data_2 %>%
    filter(A == "INCHI" |
             A == "UNIQUE-ID")
  
  data_4 <- data_3 %>%
    filter(A == "UNIQUE-ID" &
             lead(A, n = 1) == "INCHI" |
             A == "INCHI" &
             lag(A, n = 1) == "UNIQUE-ID")
  
  data_5 <- data_4 %>%
    pivot_wider(names_from = A,
                values_from = B) %>%
    unnest()
  
  data$X. <- as.character(data$X.)
  
  test <- data %>%
    filter(., grepl(pattern = "# Organism: ",
                    x = X.))
  
  test_2 <-
    data.frame(gsub(pattern = "# Organism: ",
                    replacement = "",
                    test)) %>% cSplit(splitCols = 1,
                                      sep = " ")
  
  colnames(test_2)[1] <- "A"
  
  test_2$A <- as.character(test_2$A)
  
  colnames(test_2)[2] <- "B"
  
  test_2$B <- as.character(test_2$B)
  
  biologicalsource <- paste(test_2$A,
                            test_2$B)
  
  colnames(data_5)[1] <- "uniqueid"
  
  colnames(data_5)[2] <- "inchi"
  
  data_5$biologicalsource <- biologicalsource
  
  data_standard <- data_5
  
  biologicalsource_2 <- gsub(pattern = " ",
                             replacement = "_",
                             x = biologicalsource)
  
  outpath <- file.path(pathDataExternalDbSourcePlantcyc,
                       paste(biologicalsource_2,
                             ".tsv.zip",
                             sep = ""))
  
  write.table(
    x = data_standard,
    file = gzfile(description = outpath,
                  compression = 9,
                  encoding = "UTF-8"),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t",
    fileEncoding = "UTF-8"
  )
}

#######################################################
#######################################################

distinct_biosources <- function(x)
{
  newdf  <- x %>%
    filter(!is.na(organism_lowertaxon)) %>%
    distinct(organism_lowertaxon,
             .keep_all = TRUE) %>%
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
#######################################################

distinct_pairs <- function(x)
{
  newdf  <- x %>%
    filter(
      !is.na(structureCurated) &
        !is.na(organism_lowertaxon) &
        !is.na(referenceOriginal) |
        # this will have to be adapted later on
        database == "dnp_1"
    ) %>% # this will have to be adapted later on
    distinct(structureCurated,
             organism_lowertaxon,
             .keep_all = TRUE) %>%
    group_by(structureCurated,
             organism_7_species) %>%
    add_count() %>%
    ungroup() %>%
    filter(!is.na(organism_lowertaxon) |
             !n > 1) %>%
    select(-n) %>%
    group_by(structureCurated,
             organism_6_genus) %>%
    add_count() %>%
    ungroup() %>%
    filter(!is.na(organism_lowertaxon) |
             !n > 1) %>%
    select(-n) %>%
    group_by(structureCurated,
             organism_5_family) %>%
    add_count() %>%
    ungroup() %>%
    filter(!is.na(organism_lowertaxon) |
             !n > 1) %>%
    select(-n) %>%
    group_by(structureCurated,
             organism_4_order) %>%
    add_count() %>%
    ungroup() %>%
    filter(!is.na(organism_lowertaxon) |
             !n > 1) %>%
    select(-n) %>%
    group_by(structureCurated,
             organism_3_class) %>%
    add_count() %>%
    ungroup() %>%
    filter(!is.na(organism_lowertaxon) |
             !n > 1) %>%
    select(-n) %>%
    group_by(structureCurated,
             organism_2_phylum) %>%
    add_count() %>%
    ungroup() %>%
    filter(!is.na(organism_lowertaxon) |
             !n > 1) %>%
    select(-n) %>%
    group_by(structureCurated,
             organism_1_kingdom) %>%
    add_count() %>%
    ungroup() %>%
    filter(!is.na(organism_lowertaxon) |
             !n > 1) %>%
    select(-n)
  
  return(newdf)
  
}

#######################################################
#######################################################

tcm_standardizing <- function(x)
{
  data_bulbus <- x %>%
    filter(grepl("^Bulbus", biologicalsource))
  data_bulbus$newbiologicalsource <-
    gsub("Bulbus", "\\1", data_bulbus$biologicalsource)
  data_bulbus$newbiologicalsource <-
    trimws(data_bulbus$newbiologicalsource)
  data_bulbus$newbiologicalsource <-
    capitalize(data_bulbus$newbiologicalsource)
  if (nrow(data_bulbus) != 0)
    data_bulbus$newbiologicalsource <-
    paste(data_bulbus$newbiologicalsource, "bulbus")
  
  data_caulis <- x %>%
    filter(grepl("^Caulis", biologicalsource))
  data_caulis$newbiologicalsource <-
    gsub("Caulis", "\\1", data_caulis$biologicalsource)
  data_caulis$newbiologicalsource <-
    trimws(data_caulis$newbiologicalsource)
  data_caulis$newbiologicalsource <-
    capitalize(data_caulis$newbiologicalsource)
  if (nrow(data_caulis) != 0)
    data_caulis$newbiologicalsource <-
    paste(data_caulis$newbiologicalsource, "caulis")
  
  data_caluis_et_folium <- x %>%
    filter(grepl("^Caluis_et_folium", biologicalsource))
  data_caluis_et_folium$newbiologicalsource <-
    gsub("Caluis et folium",
         "\\1",
         data_caluis_et_folium$biologicalsource)
  data_caluis_et_folium$newbiologicalsource <-
    trimws(data_caluis_et_folium$newbiologicalsource)
  data_caluis_et_folium$newbiologicalsource <-
    capitalize(data_caluis_et_folium$newbiologicalsource)
  if (nrow(data_caluis_et_folium) != 0)
    data_caluis_et_folium$newbiologicalsource <-
    paste(data_caluis_et_folium$newbiologicalsource,
          "caluis et folium")
  
  data_corolla <- x %>%
    filter(grepl("^Corolla", biologicalsource))
  data_corolla$newbiologicalsource <-
    gsub("Corolla", "\\1", data_corolla$biologicalsource)
  data_corolla$newbiologicalsource <-
    trimws(data_corolla$newbiologicalsource)
  data_corolla$newbiologicalsource <-
    capitalize(data_corolla$newbiologicalsource)
  if (nrow(data_corolla) != 0)
    data_corolla$newbiologicalsource <-
    paste(data_corolla$newbiologicalsource, "corolla")
  
  data_cortex <- x %>%
    filter(grepl("^Cortex", biologicalsource))
  data_cortex$newbiologicalsource <-
    gsub("Cortex", "\\1", data_cortex$biologicalsource)
  data_cortex$newbiologicalsource <-
    trimws(data_cortex$newbiologicalsource)
  data_cortex$newbiologicalsource <-
    capitalize(data_cortex$newbiologicalsource)
  if (nrow(data_cortex) != 0)
    data_cortex$newbiologicalsource <-
    paste(data_cortex$newbiologicalsource, "cortex")
  
  data_exocarpium <- x %>%
    filter(grepl("^Exocarpium", biologicalsource))
  data_exocarpium$newbiologicalsource <-
    gsub("Exocarpium", "\\1", data_exocarpium$biologicalsource)
  data_exocarpium$newbiologicalsource <-
    trimws(data_exocarpium$newbiologicalsource)
  data_exocarpium$newbiologicalsource <-
    capitalize(data_exocarpium$newbiologicalsource)
  if (nrow(data_exocarpium) != 0)
    data_exocarpium$newbiologicalsource <-
    paste(data_exocarpium$newbiologicalsource, "exocarpium")
  
  data_exocarpium_rubrum <- x %>%
    filter(grepl("^Exocarpium rubrum", biologicalsource))
  data_exocarpium_rubrum$newbiologicalsource <-
    gsub("Exocarpium rubrum",
         "\\1",
         data_exocarpium_rubrum$biologicalsource)
  data_exocarpium_rubrum$newbiologicalsource <-
    trimws(data_exocarpium_rubrum$newbiologicalsource)
  data_exocarpium_rubrum$newbiologicalsource <-
    capitalize(data_exocarpium_rubrum$newbiologicalsource)
  if (nrow(data_exocarpium_rubrum) != 0)
    data_exocarpium_rubrum$newbiologicalsource <-
    paste(data_exocarpium_rubrum$newbiologicalsource,
          "exocarpium rubrum")
  
  data_flos <- x %>%
    filter(grepl("^Flos", biologicalsource))
  data_flos$newbiologicalsource <-
    gsub("Flos", "\\1", data_flos$biologicalsource)
  data_flos$newbiologicalsource <-
    trimws(data_flos$newbiologicalsource)
  data_flos$newbiologicalsource <-
    capitalize(data_flos$newbiologicalsource)
  if (nrow(data_flos) != 0)
    data_flos$newbiologicalsource <-
    paste(data_flos$newbiologicalsource, "flos")
  
  data_folium <- x %>%
    filter(grepl("^Folium", biologicalsource))
  data_folium$newbiologicalsource <-
    gsub("Folium", "\\1", data_folium$biologicalsource)
  data_folium$newbiologicalsource <-
    trimws(data_folium$newbiologicalsource)
  data_folium$newbiologicalsource <-
    capitalize(data_folium$newbiologicalsource)
  if (nrow(data_folium) != 0)
    data_folium$newbiologicalsource <-
    paste(data_folium$newbiologicalsource, "folium")
  
  data_folium_et_cacumen <- x %>%
    filter(grepl("^Folium et cacumen", biologicalsource))
  data_folium_et_cacumen$newbiologicalsource <-
    gsub("Folium et cacumen",
         "\\1",
         data_folium_et_cacumen$biologicalsource)
  data_folium_et_cacumen$newbiologicalsource <-
    trimws(data_folium_et_cacumen$newbiologicalsource)
  data_folium_et_cacumen$newbiologicalsource <-
    capitalize(data_folium_et_cacumen$newbiologicalsource)
  if (nrow(data_folium_et_cacumen) != 0)
    data_folium_et_cacumen$newbiologicalsource <-
    paste(data_folium_et_cacumen$newbiologicalsource,
          "folium et cacumen")
  
  data_folium_et_caulis <- x %>%
    filter(grepl("^Folium et caulis", biologicalsource))
  data_folium_et_caulis$newbiologicalsource <-
    gsub("Folium et caulis",
         "\\1",
         data_folium_et_caulis$biologicalsource)
  data_folium_et_caulis$newbiologicalsource <-
    trimws(data_folium_et_caulis$newbiologicalsource)
  data_folium_et_caulis$newbiologicalsource <-
    capitalize(data_folium_et_caulis$newbiologicalsource)
  if (nrow(data_folium_et_caulis) != 0)
    data_folium_et_caulis$newbiologicalsource <-
    paste(data_folium_et_caulis$newbiologicalsource,
          "folium et caulis")
  
  data_fructus <- x %>%
    filter(grepl("^Fructus", biologicalsource))
  data_fructus$newbiologicalsource <-
    gsub("Fructus", "\\1", data_fructus$biologicalsource)
  data_fructus$newbiologicalsource <-
    trimws(data_fructus$newbiologicalsource)
  data_fructus$newbiologicalsource <-
    capitalize(data_fructus$newbiologicalsource)
  if (nrow(data_fructus) != 0)
    data_fructus$newbiologicalsource <-
    paste(data_fructus$newbiologicalsource, "fructus")
  
  data_fructus_germinatus <- x %>%
    filter(grepl("^Fructus germinatus", biologicalsource))
  data_fructus_germinatus$newbiologicalsource <-
    gsub("Fructus germinatus",
         "\\1",
         data_fructus_germinatus$biologicalsource)
  data_fructus_germinatus$newbiologicalsource <-
    trimws(data_fructus_germinatus$newbiologicalsource)
  data_fructus_germinatus$newbiologicalsource <-
    capitalize(data_fructus_germinatus$newbiologicalsource)
  if (nrow(data_fructus_germinatus) != 0)
    data_fructus_germinatus$newbiologicalsource <-
    paste(data_fructus_germinatus$newbiologicalsource,
          "fructus germinatus")
  
  data_fructus_immaturus <- x %>%
    filter(grepl("^Fructus immaturus", biologicalsource))
  data_fructus_immaturus$newbiologicalsource <-
    gsub("Fructus immaturus",
         "\\1",
         data_fructus_immaturus$biologicalsource)
  data_fructus_immaturus$newbiologicalsource <-
    trimws(data_fructus_immaturus$newbiologicalsource)
  data_fructus_immaturus$newbiologicalsource <-
    capitalize(data_fructus_immaturus$newbiologicalsource)
  if (nrow(data_fructus_immaturus) != 0)
    data_fructus_immaturus$newbiologicalsource <-
    paste(data_fructus_immaturus$newbiologicalsource,
          "fructus immaturus")
  
  data_fructus_retinervus <- x %>%
    filter(grepl("^Fructus retinervus", biologicalsource))
  data_fructus_retinervus$newbiologicalsource <-
    gsub("Fructus retinervus",
         "\\1",
         data_fructus_retinervus$biologicalsource)
  data_fructus_retinervus$newbiologicalsource <-
    trimws(data_fructus_retinervus$newbiologicalsource)
  data_fructus_retinervus$newbiologicalsource <-
    capitalize(data_fructus_retinervus$newbiologicalsource)
  if (nrow(data_fructus_retinervus) != 0)
    data_fructus_retinervus$newbiologicalsource <-
    paste(data_fructus_retinervus$newbiologicalsource,
          "fructus retinervus")
  
  data_fructus_rotundus <- x %>%
    filter(grepl("^Fructus rotundus", biologicalsource))
  data_fructus_rotundus$newbiologicalsource <-
    gsub("Fructus rotundus",
         "\\1",
         data_fructus_rotundus$biologicalsource)
  data_fructus_rotundus$newbiologicalsource <-
    trimws(data_fructus_rotundus$newbiologicalsource)
  data_fructus_rotundus$newbiologicalsource <-
    capitalize(data_fructus_rotundus$newbiologicalsource)
  if (nrow(data_fructus_rotundus) != 0)
    data_fructus_rotundus$newbiologicalsource <-
    paste(data_fructus_rotundus$newbiologicalsource,
          "fructus rotundus")
  
  data_herba <- x %>%
    filter(grepl("^Herba", biologicalsource))
  data_herba$newbiologicalsource <-
    gsub("Herba", "\\1", data_herba$biologicalsource)
  data_herba$newbiologicalsource <-
    trimws(data_herba$newbiologicalsource)
  data_herba$newbiologicalsource <-
    capitalize(data_herba$newbiologicalsource)
  if (nrow(data_herba) != 0)
    data_herba$newbiologicalsource <-
    paste(data_herba$newbiologicalsource, "herba")
  
  data_lignum <- x %>%
    filter(grepl("^Lignum", biologicalsource))
  data_lignum$newbiologicalsource <-
    gsub("Lignum", "\\1", data_lignum$biologicalsource)
  data_lignum$newbiologicalsource <-
    trimws(data_lignum$newbiologicalsource)
  data_lignum$newbiologicalsource <-
    capitalize(data_lignum$newbiologicalsource)
  if (nrow(data_lignum) != 0)
    data_lignum$newbiologicalsource <-
    paste(data_lignum$newbiologicalsource, "lignum")
  
  data_medulla <- x %>%
    filter(grepl("^Medulla", biologicalsource))
  data_medulla$newbiologicalsource <-
    gsub("Medulla", "\\1", data_medulla$biologicalsource)
  data_medulla$newbiologicalsource <-
    trimws(data_medulla$newbiologicalsource)
  data_medulla$newbiologicalsource <-
    capitalize(data_medulla$newbiologicalsource)
  if (nrow(data_medulla) != 0)
    data_medulla$newbiologicalsource <-
    paste(data_medulla$newbiologicalsource, "medulla")
  
  data_pericarpum <- x %>%
    filter(grepl("^Pericarpum", biologicalsource))
  data_pericarpum$newbiologicalsource <-
    gsub("Pericarpum", "\\1", data_pericarpum$biologicalsource)
  data_pericarpum$newbiologicalsource <-
    trimws(data_pericarpum$newbiologicalsource)
  data_pericarpum$newbiologicalsource <-
    capitalize(data_pericarpum$newbiologicalsource)
  if (nrow(data_pericarpum) != 0)
    data_pericarpum$newbiologicalsource <-
    paste(data_pericarpum$newbiologicalsource, "pericarpum")
  
  data_petiolus <- x %>%
    filter(grepl("^Petiolus", biologicalsource))
  data_petiolus$newbiologicalsource <-
    gsub("Petiolus", "\\1", data_petiolus$biologicalsource)
  data_petiolus$newbiologicalsource <-
    trimws(data_petiolus$newbiologicalsource)
  data_petiolus$newbiologicalsource <-
    capitalize(data_petiolus$newbiologicalsource)
  if (nrow(data_petiolus) != 0)
    data_petiolus$newbiologicalsource <-
    paste(data_petiolus$newbiologicalsource, "petiolus")
  
  data_pollen <- x %>%
    filter(grepl("^Pollen", biologicalsource))
  data_pollen$newbiologicalsource <-
    gsub("Pollen", "\\1", data_pollen$biologicalsource)
  data_pollen$newbiologicalsource <-
    trimws(data_pollen$newbiologicalsource)
  data_pollen$newbiologicalsource <-
    capitalize(data_pollen$newbiologicalsource)
  if (nrow(data_pollen) != 0)
    data_pollen$newbiologicalsource <-
    paste(data_pollen$newbiologicalsource, "pollen")
  
  data_radicis_cortex <- x %>%
    filter(grepl("^Radicis cortex", biologicalsource))
  data_radicis_cortex$newbiologicalsource <-
    gsub("Radicis cortex",
         "\\1",
         data_radicis_cortex$biologicalsource)
  data_radicis_cortex$newbiologicalsource <-
    trimws(data_radicis_cortex$newbiologicalsource)
  data_radicis_cortex$newbiologicalsource <-
    capitalize(data_radicis_cortex$newbiologicalsource)
  if (nrow(data_radicis_cortex) != 0)
    data_radicis_cortex$newbiologicalsource <-
    paste(data_radicis_cortex$newbiologicalsource, "radicis cortex")
  
  data_radix <- x %>%
    filter(grepl("^Radix", biologicalsource))
  data_radix$newbiologicalsource <-
    gsub("Radix", "\\1", data_radix$biologicalsource)
  data_radix$newbiologicalsource <-
    trimws(data_radix$newbiologicalsource)
  data_radix$newbiologicalsource <-
    capitalize(data_radix$newbiologicalsource)
  if (nrow(data_radix) != 0)
    data_radix$newbiologicalsource <-
    paste(data_radix$newbiologicalsource, "radix")
  
  data_radix_et_rhizoma <- x %>%
    filter(grepl("^Radix et rhizoma", biologicalsource))
  data_radix_et_rhizoma$newbiologicalsource <-
    gsub("Radix et rhizoma",
         "\\1",
         data_radix_et_rhizoma$biologicalsource)
  data_radix_et_rhizoma$newbiologicalsource <-
    trimws(data_radix_et_rhizoma$newbiologicalsource)
  data_radix_et_rhizoma$newbiologicalsource <-
    capitalize(data_radix_et_rhizoma$newbiologicalsource)
  if (nrow(data_radix_et_rhizoma) != 0)
    data_radix_et_rhizoma$newbiologicalsource <-
    paste(data_radix_et_rhizoma$newbiologicalsource,
          "radix et rhizoma")
  
  data_radix_preparata <- x %>%
    filter(grepl("^Radix preparata", biologicalsource))
  data_radix_preparata$newbiologicalsource <-
    gsub("Radix preparata",
         "\\1",
         data_radix_preparata$biologicalsource)
  data_radix_preparata$newbiologicalsource <-
    trimws(data_radix_preparata$newbiologicalsource)
  data_radix_preparata$newbiologicalsource <-
    capitalize(data_radix_preparata$newbiologicalsource)
  if (nrow(data_radix_preparata) != 0)
    data_radix_preparata$newbiologicalsource <-
    paste(data_radix_preparata$newbiologicalsource,
          "radix preparata")
  
  data_ramulus <- x %>%
    filter(grepl("^Ramulus", biologicalsource))
  data_ramulus$newbiologicalsource <-
    gsub("Ramulus", "\\1", data_ramulus$biologicalsource)
  data_ramulus$newbiologicalsource <-
    trimws(data_ramulus$newbiologicalsource)
  data_ramulus$newbiologicalsource <-
    capitalize(data_ramulus$newbiologicalsource)
  if (nrow(data_ramulus) != 0)
    data_ramulus$newbiologicalsource <-
    paste(data_ramulus$newbiologicalsource, "ramulus")
  
  data_ramulus_cum_uncus <- x %>%
    filter(grepl("^Ramulus cum uncus", biologicalsource))
  data_ramulus_cum_uncus$newbiologicalsource <-
    gsub("Ramulus cum uncus",
         "\\1",
         data_ramulus_cum_uncus$biologicalsource)
  data_ramulus_cum_uncus$newbiologicalsource <-
    trimws(data_ramulus_cum_uncus$newbiologicalsource)
  data_ramulus_cum_uncus$newbiologicalsource <-
    capitalize(data_ramulus_cum_uncus$newbiologicalsource)
  if (nrow(data_ramulus_cum_uncus) != 0)
    data_ramulus_cum_uncus$newbiologicalsource <-
    paste(data_ramulus_cum_uncus$newbiologicalsource,
          "ramulus cum uncus")
  
  data_ramulus_et_folium <- x %>%
    filter(grepl("^Ramulus et folium", biologicalsource))
  data_ramulus_et_folium$newbiologicalsource <-
    gsub("Ramulus et folium",
         "\\1",
         data_ramulus_et_folium$biologicalsource)
  data_ramulus_et_folium$newbiologicalsource <-
    trimws(data_ramulus_et_folium$newbiologicalsource)
  data_ramulus_et_folium$newbiologicalsource <-
    capitalize(data_ramulus_et_folium$newbiologicalsource)
  if (nrow(data_ramulus_et_folium) != 0)
    data_ramulus_et_folium$newbiologicalsource <-
    paste(data_ramulus_et_folium$newbiologicalsource,
          "ramulus et folium")
  
  data_rhizoma <- x %>%
    filter(grepl("^Rhizoma", biologicalsource))
  data_rhizoma$newbiologicalsource <-
    gsub("Rhizoma", "\\1", data_rhizoma$biologicalsource)
  data_rhizoma$newbiologicalsource <-
    trimws(data_rhizoma$newbiologicalsource)
  data_rhizoma$newbiologicalsource <-
    capitalize(data_rhizoma$newbiologicalsource)
  if (nrow(data_rhizoma) != 0)
    data_rhizoma$newbiologicalsource <-
    paste(data_rhizoma$newbiologicalsource, "rhizoma")
  
  data_rhizoma_alba <- x %>%
    filter(grepl("^Rhizoma alba", biologicalsource))
  data_rhizoma_alba$newbiologicalsource <-
    gsub("Rhizoma alba", "\\1", data_rhizoma_alba$biologicalsource)
  data_rhizoma_alba$newbiologicalsource <-
    trimws(data_rhizoma_alba$newbiologicalsource)
  data_rhizoma_alba$newbiologicalsource <-
    capitalize(data_rhizoma_alba$newbiologicalsource)
  if (nrow(data_rhizoma_alba) != 0)
    data_rhizoma_alba$newbiologicalsource <-
    paste(data_rhizoma_alba$newbiologicalsource, "rhizoma alba")
  
  data_rhizoma_et_radix <- x %>%
    filter(grepl("^Rhizoma et radix", biologicalsource))
  data_rhizoma_et_radix$newbiologicalsource <-
    gsub("Rhizoma et radix",
         "\\1",
         data_rhizoma_et_radix$biologicalsource)
  data_rhizoma_et_radix$newbiologicalsource <-
    trimws(data_rhizoma_et_radix$newbiologicalsource)
  data_rhizoma_et_radix$newbiologicalsource <-
    capitalize(data_rhizoma_et_radix$newbiologicalsource)
  if (nrow(data_rhizoma_et_radix) != 0)
    data_rhizoma_et_radix$newbiologicalsource <-
    paste(data_rhizoma_et_radix$newbiologicalsource,
          "rhizoma et radix")
  
  data_semen <- x %>%
    filter(grepl("^Semen", biologicalsource))
  data_semen$newbiologicalsource <-
    gsub("Semen", "\\1", data_semen$biologicalsource)
  data_semen$newbiologicalsource <-
    trimws(data_semen$newbiologicalsource)
  data_semen$newbiologicalsource <-
    capitalize(data_semen$newbiologicalsource)
  if (nrow(data_semen) != 0)
    data_semen$newbiologicalsource <-
    paste(data_semen$newbiologicalsource, "semen")
  
  data_semen_germinatum <- x %>%
    filter(grepl("^Semen germinatum", biologicalsource))
  data_semen_germinatum$newbiologicalsource <-
    gsub("Semen germinatum",
         "\\1",
         data_semen_germinatum$biologicalsource)
  data_semen_germinatum$newbiologicalsource <-
    trimws(data_semen_germinatum$newbiologicalsource)
  data_semen_germinatum$newbiologicalsource <-
    capitalize(data_semen_germinatum$newbiologicalsource)
  if (nrow(data_semen_germinatum) != 0)
    data_semen_germinatum$newbiologicalsource <-
    paste(data_semen_germinatum$newbiologicalsource,
          "semen germinatum")
  
  data_spica <- x %>%
    filter(grepl("Spica ", biologicalsource, fixed = TRUE))
  data_spica$newbiologicalsource <-
    gsub("Spica", "\\1", data_spica$biologicalsource)
  data_spica$newbiologicalsource <-
    trimws(data_spica$newbiologicalsource)
  data_spica$newbiologicalsource <-
    capitalize(data_spica$newbiologicalsource)
  if (nrow(data_spica) != 0)
    data_spica$newbiologicalsource <-
    paste(data_spica$newbiologicalsource, "spica")
  
  data_stamen <- x %>%
    filter(grepl("^Stamen", biologicalsource))
  data_stamen$newbiologicalsource <-
    gsub("Stamen", "\\1", data_stamen$biologicalsource)
  data_stamen$newbiologicalsource <-
    trimws(data_stamen$newbiologicalsource)
  data_stamen$newbiologicalsource <-
    capitalize(data_stamen$newbiologicalsource)
  if (nrow(data_stamen) != 0)
    data_stamen$newbiologicalsource <-
    paste(data_stamen$newbiologicalsource, "stamen")
  
  data_stigma <- x %>%
    filter(grepl("Stigma ", biologicalsource, fixed = TRUE))
  data_stigma$newbiologicalsource <-
    gsub("Stigma", "\\1", data_stigma$biologicalsource)
  data_stigma$newbiologicalsource <-
    trimws(data_stigma$newbiologicalsource)
  data_stigma$newbiologicalsource <-
    capitalize(data_stigma$newbiologicalsource)
  if (nrow(data_stigma) != 0)
    data_stigma$newbiologicalsource <-
    paste(data_stigma$newbiologicalsource, "stigma")
  
  data_storax <- x %>%
    filter(grepl("^Storax", biologicalsource))
  data_storax$newbiologicalsource <-
    gsub("Storax", "\\1", data_storax$biologicalsource)
  data_storax$newbiologicalsource <-
    trimws(data_storax$newbiologicalsource)
  data_storax$newbiologicalsource <-
    capitalize(data_storax$newbiologicalsource)
  if (nrow(data_storax) != 0)
    data_storax$newbiologicalsource <-
    paste(data_storax$newbiologicalsource, "storax")
  
  data_thallus <- x %>%
    filter(grepl("^Thallus", biologicalsource, fixed = TRUE))
  data_thallus$newbiologicalsource <-
    gsub("Thallus", "\\1", data_thallus$biologicalsource)
  data_thallus$newbiologicalsource <-
    trimws(data_thallus$newbiologicalsource)
  data_thallus$newbiologicalsource <-
    capitalize(data_thallus$newbiologicalsource)
  if (nrow(data_thallus) != 0)
    data_thallus$newbiologicalsource <-
    paste(data_thallus$newbiologicalsource, "thallus")
  #not tuber
  
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
    apply(x_3, 1, function(x)
      tail(na.omit(x), 1))
  
  x_3 <- x_3 %>%
    select(biologicalsource,
           newnewbiologicalsource)
  
  x_4 <-
    left_join(x_3, food_names_list, by = c("newnewbiologicalsource" = "name"))
  
  x_4 <-
    left_join(x_4, tcm_names_list, by = c("newnewbiologicalsource" = "latin"))
  
  x_4 <-
    left_join(x_4,
              tcm_names_list,
              by = c("newnewbiologicalsource" = "common"))
  
  x_4$newnewnewbiologicalsource <-
    apply(x_4, 1, function(x)
      tail(na.omit(x), 1))
  
  x_4 <- x_4 %>%
    select(biologicalsource = biologicalsource.x,
           newnewnewbiologicalsource)
  
  data_standard_2 <- left_join(data_standard, x_4) %>%
    select(-biologicalsource) %>%
    select(biologicalsource = newnewnewbiologicalsource,
           everything())
  
  data_standard_2
}

#######################################################
#######################################################

name2inchi <- function(i)
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
{
  tryCatch({
    cpd <- dataTranslatedNominal[i, "nameSanitized"]
    url <-
      paste("https://cactus.nci.nih.gov/chemical/structure/",
            cpd,
            "/stdinchi",
            sep = "")
    url <- gsub(pattern = "\\s",
                replacement = "%20",
                x = url)
    read_html(url) %>%
      html_text()
  }
  , error = function(e) {
    NA
  })
}

#######################################################
#######################################################

tcm_inverting <- function(x)
{
  data_bulbus <- x %>%
    filter(grepl("^Bulbus", biologicalsource))
  data_bulbus$newbiologicalsource <-
    gsub("Bulbus", "\\1", data_bulbus$biologicalsource)
  data_bulbus$newbiologicalsource <-
    trimws(data_bulbus$newbiologicalsource)
  data_bulbus$newbiologicalsource <-
    capitalize(data_bulbus$newbiologicalsource)
  if (nrow(data_bulbus) != 0)
    data_bulbus$newbiologicalsource <-
    paste(data_bulbus$newbiologicalsource, "bulbus")
  
  data_caulis <- x %>%
    filter(grepl("^Caulis", biologicalsource))
  data_caulis$newbiologicalsource <-
    gsub("Caulis", "\\1", data_caulis$biologicalsource)
  data_caulis$newbiologicalsource <-
    trimws(data_caulis$newbiologicalsource)
  data_caulis$newbiologicalsource <-
    capitalize(data_caulis$newbiologicalsource)
  if (nrow(data_caulis) != 0)
    data_caulis$newbiologicalsource <-
    paste(data_caulis$newbiologicalsource, "caulis")
  
  data_caluis_et_folium <- x %>%
    filter(grepl("^Caluis_et_folium", biologicalsource))
  data_caluis_et_folium$newbiologicalsource <-
    gsub("Caluis et folium",
         "\\1",
         data_caluis_et_folium$biologicalsource)
  data_caluis_et_folium$newbiologicalsource <-
    trimws(data_caluis_et_folium$newbiologicalsource)
  data_caluis_et_folium$newbiologicalsource <-
    capitalize(data_caluis_et_folium$newbiologicalsource)
  if (nrow(data_caluis_et_folium) != 0)
    data_caluis_et_folium$newbiologicalsource <-
    paste(data_caluis_et_folium$newbiologicalsource,
          "caluis et folium")
  
  data_corolla <- x %>%
    filter(grepl("^Corolla", biologicalsource))
  data_corolla$newbiologicalsource <-
    gsub("Corolla", "\\1", data_corolla$biologicalsource)
  data_corolla$newbiologicalsource <-
    trimws(data_corolla$newbiologicalsource)
  data_corolla$newbiologicalsource <-
    capitalize(data_corolla$newbiologicalsource)
  if (nrow(data_corolla) != 0)
    data_corolla$newbiologicalsource <-
    paste(data_corolla$newbiologicalsource, "corolla")
  
  data_cortex <- x %>%
    filter(grepl("^Cortex", biologicalsource))
  data_cortex$newbiologicalsource <-
    gsub("Cortex", "\\1", data_cortex$biologicalsource)
  data_cortex$newbiologicalsource <-
    trimws(data_cortex$newbiologicalsource)
  data_cortex$newbiologicalsource <-
    capitalize(data_cortex$newbiologicalsource)
  if (nrow(data_cortex) != 0)
    data_cortex$newbiologicalsource <-
    paste(data_cortex$newbiologicalsource, "cortex")
  
  data_exocarpium <- x %>%
    filter(grepl("^Exocarpium", biologicalsource))
  data_exocarpium$newbiologicalsource <-
    gsub("Exocarpium", "\\1", data_exocarpium$biologicalsource)
  data_exocarpium$newbiologicalsource <-
    trimws(data_exocarpium$newbiologicalsource)
  data_exocarpium$newbiologicalsource <-
    capitalize(data_exocarpium$newbiologicalsource)
  if (nrow(data_exocarpium) != 0)
    data_exocarpium$newbiologicalsource <-
    paste(data_exocarpium$newbiologicalsource, "exocarpium")
  
  data_exocarpium_rubrum <- x %>%
    filter(grepl("^Exocarpium rubrum", biologicalsource))
  data_exocarpium_rubrum$newbiologicalsource <-
    gsub("Exocarpium rubrum",
         "\\1",
         data_exocarpium_rubrum$biologicalsource)
  data_exocarpium_rubrum$newbiologicalsource <-
    trimws(data_exocarpium_rubrum$newbiologicalsource)
  data_exocarpium_rubrum$newbiologicalsource <-
    capitalize(data_exocarpium_rubrum$newbiologicalsource)
  if (nrow(data_exocarpium_rubrum) != 0)
    data_exocarpium_rubrum$newbiologicalsource <-
    paste(data_exocarpium_rubrum$newbiologicalsource,
          "exocarpium rubrum")
  
  data_flos <- x %>%
    filter(grepl("^Flos", biologicalsource))
  data_flos$newbiologicalsource <-
    gsub("Flos", "\\1", data_flos$biologicalsource)
  data_flos$newbiologicalsource <-
    trimws(data_flos$newbiologicalsource)
  data_flos$newbiologicalsource <-
    capitalize(data_flos$newbiologicalsource)
  if (nrow(data_flos) != 0)
    data_flos$newbiologicalsource <-
    paste(data_flos$newbiologicalsource, "flos")
  
  data_folium <- x %>%
    filter(grepl("^Folium", biologicalsource))
  data_folium$newbiologicalsource <-
    gsub("Folium", "\\1", data_folium$biologicalsource)
  data_folium$newbiologicalsource <-
    trimws(data_folium$newbiologicalsource)
  data_folium$newbiologicalsource <-
    capitalize(data_folium$newbiologicalsource)
  if (nrow(data_folium) != 0)
    data_folium$newbiologicalsource <-
    paste(data_folium$newbiologicalsource, "folium")
  
  data_folium_et_cacumen <- x %>%
    filter(grepl("^Folium et cacumen", biologicalsource))
  data_folium_et_cacumen$newbiologicalsource <-
    gsub("Folium et cacumen",
         "\\1",
         data_folium_et_cacumen$biologicalsource)
  data_folium_et_cacumen$newbiologicalsource <-
    trimws(data_folium_et_cacumen$newbiologicalsource)
  data_folium_et_cacumen$newbiologicalsource <-
    capitalize(data_folium_et_cacumen$newbiologicalsource)
  if (nrow(data_folium_et_cacumen) != 0)
    data_folium_et_cacumen$newbiologicalsource <-
    paste(data_folium_et_cacumen$newbiologicalsource,
          "folium et cacumen")
  
  data_folium_et_caulis <- x %>%
    filter(grepl("^Folium et caulis", biologicalsource))
  data_folium_et_caulis$newbiologicalsource <-
    gsub("Folium et caulis",
         "\\1",
         data_folium_et_caulis$biologicalsource)
  data_folium_et_caulis$newbiologicalsource <-
    trimws(data_folium_et_caulis$newbiologicalsource)
  data_folium_et_caulis$newbiologicalsource <-
    capitalize(data_folium_et_caulis$newbiologicalsource)
  if (nrow(data_folium_et_caulis) != 0)
    data_folium_et_caulis$newbiologicalsource <-
    paste(data_folium_et_caulis$newbiologicalsource,
          "folium et caulis")
  
  data_fructus <- x %>%
    filter(grepl("^Fructus", biologicalsource))
  data_fructus$newbiologicalsource <-
    gsub("Fructus", "\\1", data_fructus$biologicalsource)
  data_fructus$newbiologicalsource <-
    trimws(data_fructus$newbiologicalsource)
  data_fructus$newbiologicalsource <-
    capitalize(data_fructus$newbiologicalsource)
  if (nrow(data_fructus) != 0)
    data_fructus$newbiologicalsource <-
    paste(data_fructus$newbiologicalsource, "fructus")
  
  data_fructus_germinatus <- x %>%
    filter(grepl("^Fructus germinatus", biologicalsource))
  data_fructus_germinatus$newbiologicalsource <-
    gsub("Fructus germinatus",
         "\\1",
         data_fructus_germinatus$biologicalsource)
  data_fructus_germinatus$newbiologicalsource <-
    trimws(data_fructus_germinatus$newbiologicalsource)
  data_fructus_germinatus$newbiologicalsource <-
    capitalize(data_fructus_germinatus$newbiologicalsource)
  if (nrow(data_fructus_germinatus) != 0)
    data_fructus_germinatus$newbiologicalsource <-
    paste(data_fructus_germinatus$newbiologicalsource,
          "fructus germinatus")
  
  data_fructus_immaturus <- x %>%
    filter(grepl("^Fructus immaturus", biologicalsource))
  data_fructus_immaturus$newbiologicalsource <-
    gsub("Fructus immaturus",
         "\\1",
         data_fructus_immaturus$biologicalsource)
  data_fructus_immaturus$newbiologicalsource <-
    trimws(data_fructus_immaturus$newbiologicalsource)
  data_fructus_immaturus$newbiologicalsource <-
    capitalize(data_fructus_immaturus$newbiologicalsource)
  if (nrow(data_fructus_immaturus) != 0)
    data_fructus_immaturus$newbiologicalsource <-
    paste(data_fructus_immaturus$newbiologicalsource,
          "fructus immaturus")
  
  data_fructus_retinervus <- x %>%
    filter(grepl("^Fructus retinervus", biologicalsource))
  data_fructus_retinervus$newbiologicalsource <-
    gsub("Fructus retinervus",
         "\\1",
         data_fructus_retinervus$biologicalsource)
  data_fructus_retinervus$newbiologicalsource <-
    trimws(data_fructus_retinervus$newbiologicalsource)
  data_fructus_retinervus$newbiologicalsource <-
    capitalize(data_fructus_retinervus$newbiologicalsource)
  if (nrow(data_fructus_retinervus) != 0)
    data_fructus_retinervus$newbiologicalsource <-
    paste(data_fructus_retinervus$newbiologicalsource,
          "fructus retinervus")
  
  data_fructus_rotundus <- x %>%
    filter(grepl("^Fructus rotundus", biologicalsource))
  data_fructus_rotundus$newbiologicalsource <-
    gsub("Fructus rotundus",
         "\\1",
         data_fructus_rotundus$biologicalsource)
  data_fructus_rotundus$newbiologicalsource <-
    trimws(data_fructus_rotundus$newbiologicalsource)
  data_fructus_rotundus$newbiologicalsource <-
    capitalize(data_fructus_rotundus$newbiologicalsource)
  if (nrow(data_fructus_rotundus) != 0)
    data_fructus_rotundus$newbiologicalsource <-
    paste(data_fructus_rotundus$newbiologicalsource,
          "fructus rotundus")
  
  data_herba <- x %>%
    filter(grepl("^Herba", biologicalsource))
  data_herba$newbiologicalsource <-
    gsub("Herba", "\\1", data_herba$biologicalsource)
  data_herba$newbiologicalsource <-
    trimws(data_herba$newbiologicalsource)
  data_herba$newbiologicalsource <-
    capitalize(data_herba$newbiologicalsource)
  if (nrow(data_herba) != 0)
    data_herba$newbiologicalsource <-
    paste(data_herba$newbiologicalsource, "herba")
  
  data_lignum <- x %>%
    filter(grepl("^Lignum", biologicalsource))
  data_lignum$newbiologicalsource <-
    gsub("Lignum", "\\1", data_lignum$biologicalsource)
  data_lignum$newbiologicalsource <-
    trimws(data_lignum$newbiologicalsource)
  data_lignum$newbiologicalsource <-
    capitalize(data_lignum$newbiologicalsource)
  if (nrow(data_lignum) != 0)
    data_lignum$newbiologicalsource <-
    paste(data_lignum$newbiologicalsource, "lignum")
  
  data_medulla <- x %>%
    filter(grepl("^Medulla", biologicalsource))
  data_medulla$newbiologicalsource <-
    gsub("Medulla", "\\1", data_medulla$biologicalsource)
  data_medulla$newbiologicalsource <-
    trimws(data_medulla$newbiologicalsource)
  data_medulla$newbiologicalsource <-
    capitalize(data_medulla$newbiologicalsource)
  if (nrow(data_medulla) != 0)
    data_medulla$newbiologicalsource <-
    paste(data_medulla$newbiologicalsource, "medulla")
  
  data_pericarpum <- x %>%
    filter(grepl("^Pericarpum", biologicalsource))
  data_pericarpum$newbiologicalsource <-
    gsub("Pericarpum", "\\1", data_pericarpum$biologicalsource)
  data_pericarpum$newbiologicalsource <-
    trimws(data_pericarpum$newbiologicalsource)
  data_pericarpum$newbiologicalsource <-
    capitalize(data_pericarpum$newbiologicalsource)
  if (nrow(data_pericarpum) != 0)
    data_pericarpum$newbiologicalsource <-
    paste(data_pericarpum$newbiologicalsource, "pericarpum")
  
  data_petiolus <- x %>%
    filter(grepl("^Petiolus", biologicalsource))
  data_petiolus$newbiologicalsource <-
    gsub("Petiolus", "\\1", data_petiolus$biologicalsource)
  data_petiolus$newbiologicalsource <-
    trimws(data_petiolus$newbiologicalsource)
  data_petiolus$newbiologicalsource <-
    capitalize(data_petiolus$newbiologicalsource)
  if (nrow(data_petiolus) != 0)
    data_petiolus$newbiologicalsource <-
    paste(data_petiolus$newbiologicalsource, "petiolus")
  
  data_pollen <- x %>%
    filter(grepl("^Pollen", biologicalsource))
  data_pollen$newbiologicalsource <-
    gsub("Pollen", "\\1", data_pollen$biologicalsource)
  data_pollen$newbiologicalsource <-
    trimws(data_pollen$newbiologicalsource)
  data_pollen$newbiologicalsource <-
    capitalize(data_pollen$newbiologicalsource)
  if (nrow(data_pollen) != 0)
    data_pollen$newbiologicalsource <-
    paste(data_pollen$newbiologicalsource, "pollen")
  
  data_radicis_cortex <- x %>%
    filter(grepl("^Radicis cortex", biologicalsource))
  data_radicis_cortex$newbiologicalsource <-
    gsub("Radicis cortex",
         "\\1",
         data_radicis_cortex$biologicalsource)
  data_radicis_cortex$newbiologicalsource <-
    trimws(data_radicis_cortex$newbiologicalsource)
  data_radicis_cortex$newbiologicalsource <-
    capitalize(data_radicis_cortex$newbiologicalsource)
  if (nrow(data_radicis_cortex) != 0)
    data_radicis_cortex$newbiologicalsource <-
    paste(data_radicis_cortex$newbiologicalsource, "radicis cortex")
  
  data_radix <- x %>%
    filter(grepl("^Radix", biologicalsource))
  data_radix$newbiologicalsource <-
    gsub("Radix", "\\1", data_radix$biologicalsource)
  data_radix$newbiologicalsource <-
    trimws(data_radix$newbiologicalsource)
  data_radix$newbiologicalsource <-
    capitalize(data_radix$newbiologicalsource)
  if (nrow(data_radix) != 0)
    data_radix$newbiologicalsource <-
    paste(data_radix$newbiologicalsource, "radix")
  
  data_radix_et_rhizoma <- x %>%
    filter(grepl("^Radix et rhizoma", biologicalsource))
  data_radix_et_rhizoma$newbiologicalsource <-
    gsub("Radix et rhizoma",
         "\\1",
         data_radix_et_rhizoma$biologicalsource)
  data_radix_et_rhizoma$newbiologicalsource <-
    trimws(data_radix_et_rhizoma$newbiologicalsource)
  data_radix_et_rhizoma$newbiologicalsource <-
    capitalize(data_radix_et_rhizoma$newbiologicalsource)
  if (nrow(data_radix_et_rhizoma) != 0)
    data_radix_et_rhizoma$newbiologicalsource <-
    paste(data_radix_et_rhizoma$newbiologicalsource,
          "radix et rhizoma")
  
  data_radix_preparata <- x %>%
    filter(grepl("^Radix preparata", biologicalsource))
  data_radix_preparata$newbiologicalsource <-
    gsub("Radix preparata",
         "\\1",
         data_radix_preparata$biologicalsource)
  data_radix_preparata$newbiologicalsource <-
    trimws(data_radix_preparata$newbiologicalsource)
  data_radix_preparata$newbiologicalsource <-
    capitalize(data_radix_preparata$newbiologicalsource)
  if (nrow(data_radix_preparata) != 0)
    data_radix_preparata$newbiologicalsource <-
    paste(data_radix_preparata$newbiologicalsource,
          "radix preparata")
  
  data_ramulus <- x %>%
    filter(grepl("^Ramulus", biologicalsource))
  data_ramulus$newbiologicalsource <-
    gsub("Ramulus", "\\1", data_ramulus$biologicalsource)
  data_ramulus$newbiologicalsource <-
    trimws(data_ramulus$newbiologicalsource)
  data_ramulus$newbiologicalsource <-
    capitalize(data_ramulus$newbiologicalsource)
  if (nrow(data_ramulus) != 0)
    data_ramulus$newbiologicalsource <-
    paste(data_ramulus$newbiologicalsource, "ramulus")
  
  data_ramulus_cum_uncus <- x %>%
    filter(grepl("^Ramulus cum uncus", biologicalsource))
  data_ramulus_cum_uncus$newbiologicalsource <-
    gsub("Ramulus cum uncus",
         "\\1",
         data_ramulus_cum_uncus$biologicalsource)
  data_ramulus_cum_uncus$newbiologicalsource <-
    trimws(data_ramulus_cum_uncus$newbiologicalsource)
  data_ramulus_cum_uncus$newbiologicalsource <-
    capitalize(data_ramulus_cum_uncus$newbiologicalsource)
  if (nrow(data_ramulus_cum_uncus) != 0)
    data_ramulus_cum_uncus$newbiologicalsource <-
    paste(data_ramulus_cum_uncus$newbiologicalsource,
          "ramulus cum uncus")
  
  data_ramulus_et_folium <- x %>%
    filter(grepl("^Ramulus et folium", biologicalsource))
  data_ramulus_et_folium$newbiologicalsource <-
    gsub("Ramulus et folium",
         "\\1",
         data_ramulus_et_folium$biologicalsource)
  data_ramulus_et_folium$newbiologicalsource <-
    trimws(data_ramulus_et_folium$newbiologicalsource)
  data_ramulus_et_folium$newbiologicalsource <-
    capitalize(data_ramulus_et_folium$newbiologicalsource)
  if (nrow(data_ramulus_et_folium) != 0)
    data_ramulus_et_folium$newbiologicalsource <-
    paste(data_ramulus_et_folium$newbiologicalsource,
          "ramulus et folium")
  
  data_rhizoma <- x %>%
    filter(grepl("^Rhizoma", biologicalsource))
  data_rhizoma$newbiologicalsource <-
    gsub("Rhizoma", "\\1", data_rhizoma$biologicalsource)
  data_rhizoma$newbiologicalsource <-
    trimws(data_rhizoma$newbiologicalsource)
  data_rhizoma$newbiologicalsource <-
    capitalize(data_rhizoma$newbiologicalsource)
  if (nrow(data_rhizoma) != 0)
    data_rhizoma$newbiologicalsource <-
    paste(data_rhizoma$newbiologicalsource, "rhizoma")
  
  data_rhizoma_alba <- x %>%
    filter(grepl("^Rhizoma alba", biologicalsource))
  data_rhizoma_alba$newbiologicalsource <-
    gsub("Rhizoma alba", "\\1", data_rhizoma_alba$biologicalsource)
  data_rhizoma_alba$newbiologicalsource <-
    trimws(data_rhizoma_alba$newbiologicalsource)
  data_rhizoma_alba$newbiologicalsource <-
    capitalize(data_rhizoma_alba$newbiologicalsource)
  if (nrow(data_rhizoma_alba) != 0)
    data_rhizoma_alba$newbiologicalsource <-
    paste(data_rhizoma_alba$newbiologicalsource, "rhizoma alba")
  
  data_rhizoma_et_radix <- x %>%
    filter(grepl("^Rhizoma et radix", biologicalsource))
  data_rhizoma_et_radix$newbiologicalsource <-
    gsub("Rhizoma et radix",
         "\\1",
         data_rhizoma_et_radix$biologicalsource)
  data_rhizoma_et_radix$newbiologicalsource <-
    trimws(data_rhizoma_et_radix$newbiologicalsource)
  data_rhizoma_et_radix$newbiologicalsource <-
    capitalize(data_rhizoma_et_radix$newbiologicalsource)
  if (nrow(data_rhizoma_et_radix) != 0)
    data_rhizoma_et_radix$newbiologicalsource <-
    paste(data_rhizoma_et_radix$newbiologicalsource,
          "rhizoma et radix")
  
  data_semen <- x %>%
    filter(grepl("^Semen", biologicalsource))
  data_semen$newbiologicalsource <-
    gsub("Semen", "\\1", data_semen$biologicalsource)
  data_semen$newbiologicalsource <-
    trimws(data_semen$newbiologicalsource)
  data_semen$newbiologicalsource <-
    capitalize(data_semen$newbiologicalsource)
  if (nrow(data_semen) != 0)
    data_semen$newbiologicalsource <-
    paste(data_semen$newbiologicalsource, "semen")
  
  data_semen_germinatum <- x %>%
    filter(grepl("^Semen germinatum", biologicalsource))
  data_semen_germinatum$newbiologicalsource <-
    gsub("Semen germinatum",
         "\\1",
         data_semen_germinatum$biologicalsource)
  data_semen_germinatum$newbiologicalsource <-
    trimws(data_semen_germinatum$newbiologicalsource)
  data_semen_germinatum$newbiologicalsource <-
    capitalize(data_semen_germinatum$newbiologicalsource)
  if (nrow(data_semen_germinatum) != 0)
    data_semen_germinatum$newbiologicalsource <-
    paste(data_semen_germinatum$newbiologicalsource,
          "semen germinatum")
  
  data_spica <- x %>%
    filter(grepl("Spica ", biologicalsource, fixed = TRUE))
  data_spica$newbiologicalsource <-
    gsub("Spica", "\\1", data_spica$biologicalsource)
  data_spica$newbiologicalsource <-
    trimws(data_spica$newbiologicalsource)
  data_spica$newbiologicalsource <-
    capitalize(data_spica$newbiologicalsource)
  if (nrow(data_spica) != 0)
    data_spica$newbiologicalsource <-
    paste(data_spica$newbiologicalsource, "spica")
  
  data_stamen <- x %>%
    filter(grepl("^Stamen", biologicalsource))
  data_stamen$newbiologicalsource <-
    gsub("Stamen", "\\1", data_stamen$biologicalsource)
  data_stamen$newbiologicalsource <-
    trimws(data_stamen$newbiologicalsource)
  data_stamen$newbiologicalsource <-
    capitalize(data_stamen$newbiologicalsource)
  if (nrow(data_stamen) != 0)
    data_stamen$newbiologicalsource <-
    paste(data_stamen$newbiologicalsource, "stamen")
  
  data_stigma <- x %>%
    filter(grepl("Stigma ", biologicalsource, fixed = TRUE))
  data_stigma$newbiologicalsource <-
    gsub("Stigma", "\\1", data_stigma$biologicalsource)
  data_stigma$newbiologicalsource <-
    trimws(data_stigma$newbiologicalsource)
  data_stigma$newbiologicalsource <-
    capitalize(data_stigma$newbiologicalsource)
  if (nrow(data_stigma) != 0)
    data_stigma$newbiologicalsource <-
    paste(data_stigma$newbiologicalsource, "stigma")
  
  data_storax <- x %>%
    filter(grepl("^Storax", biologicalsource))
  data_storax$newbiologicalsource <-
    gsub("Storax", "\\1", data_storax$biologicalsource)
  data_storax$newbiologicalsource <-
    trimws(data_storax$newbiologicalsource)
  data_storax$newbiologicalsource <-
    capitalize(data_storax$newbiologicalsource)
  if (nrow(data_storax) != 0)
    data_storax$newbiologicalsource <-
    paste(data_storax$newbiologicalsource, "storax")
  
  data_thallus <- x %>%
    filter(grepl("^Thallus", biologicalsource, fixed = TRUE))
  data_thallus$newbiologicalsource <-
    gsub("Thallus", "\\1", data_thallus$biologicalsource)
  data_thallus$newbiologicalsource <-
    trimws(data_thallus$newbiologicalsource)
  data_thallus$newbiologicalsource <-
    capitalize(data_thallus$newbiologicalsource)
  if (nrow(data_thallus) != 0)
    data_thallus$newbiologicalsource <-
    paste(data_thallus$newbiologicalsource, "thallus")
  #not tuber
  
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
    apply(x_3, 1, function(x)
      tail(na.omit(x), 1))
  
  x_4 <- x_3 %>%
    select(-biologicalsource,-newbiologicalsource) %>%
    mutate(biologicalsource = newnewbiologicalsource) %>%
    select(-newnewbiologicalsource)
  
  return(x_4)
}

#######################################################
#######################################################

tcm_cleaning <- function(x)
{
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
  
  #not tuber
  data_final <- data %>%
    mutate(latin = newbiologicalsource) %>%
    select(latin, common, biologicalsource)
  
  return(data_final)
}


#######################################################
#######################################################

standardizing_original <- function(data_selected,
                                   db,
                                   structure_field)
{
  data_selected[setdiff(c("name",
                          "biologicalsource",
                          "reference"),
                        names(data_selected))] <- NA
  
  data_standard <- data.frame(data_selected) %>%
    mutate(database = db) %>%
    select(database,
           name,
           all_of(structure_field),
           biologicalsource,
           reference) %>%
    distinct_at(vars(all_of(structure_field),
                     biologicalsource),
                .keep_all = TRUE)
  
  data_standard[] <-
    lapply(data_standard, function(x)
      gsub("\r\n", " ", x))
  data_standard[] <-
    lapply(data_standard, function(x)
      gsub("\r", " ", x))
  data_standard[] <-
    lapply(data_standard, function(x)
      gsub("\n", " ", x))
  data_standard[] <-
    lapply(data_standard, function(x)
      gsub("\t", " ", x))
  
  return(data_standard)
}

#######################################################
#######################################################

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
#######################################################

split_data_table <-
  function(x, no_rows_per_frame, text, path_to_store) {
    split_vec <- seq(1, nrow(x), no_rows_per_frame)
    
    for (split_cut in split_vec) {
      sample <- x[split_cut:(split_cut + (no_rows_per_frame - 1))]
      write.table(
        sample,
        paste(path_to_store,
              text,
              as.integer(split_cut + (
                no_rows_per_frame - 1
              )),
              ".tsv",
              sep = ""),
        row.names = FALSE,
        quote = FALSE,
        sep = "\t",
        fileEncoding = "UTF-8"
      )
    }
  }

#######################################################
#######################################################

tcm_pharmacopoeia_cleaning <- function(x)
{
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
  
  #not tuber
  data_final <- data %>%
    select(-organismTranslated) %>%
    mutate(organismTranslated = newbiologicalsource) %>%
    select(-newbiologicalsource)
  
  return(data_final)
}

#######################################################
#######################################################

gnfinder_cleaning <- function(num, organismCol) {
  
  if (organismCol == "organismOriginal")
  inpath_organism_f <- paste(
    pathOriginalOrganismDistinct,
    "originalOrganismGnfinderUntil_",
    num,
    ".tsv",
    sep = ""
  )
  
  if (organismCol == "organismInterim")
    inpath_organism_f <- paste(
      pathTranslatedOrganismDistinct,
      "translatedOrganismGnfinderUntil_",
      num,
      ".tsv",
      sep = ""
    )
  
  if (organismCol == "organismOriginal")
  inpath_gnfinder_f <-
    paste(pathSanitizedOrganismOriginalDirJson,
          "originalOrganismGnfinderUntil_",
          num,
          ".json",
          sep = "")
  
  if (organismCol == "organismInterim")
    inpath_gnfinder_f <-
      paste(pathSanitizedOrganismTranslatedDirJson,
            "sanitizedOrganismGnfinderUntil_",
            num,
            ".json",
            sep = "")
  
  if (organismCol == "organismOriginal")
  outpath_f <-
    paste(pathSanitizedOrganismOriginalDirTsv,
          "originalOrganismGnfinderUntil_",
          num,
          ".tsv.zip",
          sep = "")
  
  if (organismCol == "organismInterim")
    outpath_f <-
      paste(pathSanitizedOrganismTranslatedDirTsv,
            "translatedOrganismGnfinderUntil_",
            num,
            ".tsv.zip",
            sep = "")
  
  gnfound <- data.frame(fromJSON(txt = inpath_gnfinder_f,
                                 simplifyDataFrame = TRUE))
  
  data_bio <- read_delim(
    file = inpath_organism_f,
    delim = "\t",
    escape_double = FALSE,
    trim_ws = FALSE
  ) %>%
    mutate_all(as.character)
  
  data_bio <- data_bio[!is.na(data_bio[,organismCol]),]
  
  data_bio_clean <- biocleaning(x = gnfound,
                                y = data_bio,
                                organismCol = organismCol)
  
  return(data_bio_clean)
}

#######################################################
#######################################################

taxo_cleaning_auto <- function(dfsel) {
  test <- dfsel %>%
    filter(!is.na(organismSanitized)) %>%
    distinct(organismSanitized,
             organism_database,
             .keep_all = TRUE)
  
  test$organism_1_kingdom <-
    y_as_na(x = test$organism_1_kingdom,
            y = "Not assigned")
  
  test$organism_2_phylum <-
    y_as_na(x = test$organism_2_phylum,
            y = "Not assigned")
  
  test$organism_3_class <-
    y_as_na(x = test$organism_3_class,
            y = "Not assigned")
  
  test$organism_4_order <-
    y_as_na(x = test$organism_4_order,
            y = "Not assigned")
  
  test$organism_5_family <-
    y_as_na(x = test$organism_5_family,
            y = "Not assigned")
  
  test$organism_6_genus <-
    y_as_na(x = test$organism_6_genus,
            y = "Not assigned")
  
  test$organism_7_species <-
    y_as_na(x = test$organism_7_species,
            y = "Not assigned")
  
  species <- test %>%
    filter(!is.na(organism_7_species)) %>%
    arrange(match(x = organism_database,
                  table =  c("Catalogue of Life",
                             "NCBI"))) %>%
    group_by(organism_1_kingdom,
             organism_7_species) %>%
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
             .keep_all = TRUE) %>%
    select(-organismOriginal, -organismTranslated, -organismSanitized)
  
  species_fill <- test %>%
    filter(!is.na(organism_7_species)) %>%
    select(organismOriginal,
           organismTranslated,
           organismSanitized,
           organism_7_species)
  
  species_full <- left_join(species_fill, species)
  
  unspecies <- test %>%
    filter(is.na(organism_7_species))
  
  genus_1 <- rbind(species_full, unspecies)
  
  genus <- genus_1 %>%
    filter(!is.na(organism_6_genus)) %>%
    arrange(match(x = organism_database,
                  table = c("Catalogue of Life",
                            "NCBI"))) %>%
    group_by(organism_1_kingdom,
             organism_6_genus) %>%
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
             .keep_all = TRUE) %>%
    select(-organismOriginal,
           -organismTranslated,
           -organismSanitized,
           -organism_7_species)
  
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
    arrange(match(x = organism_database,
                  table = c("Catalogue of Life",
                            "NCBI"))) %>%
    group_by(organism_1_kingdom,
             organism_5_family) %>%
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
             .keep_all = TRUE) %>%
    select(
      -organismOriginal,-organismTranslated,-organismSanitized,-organism_7_species,-organism_6_genus
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
    arrange(match(x = organism_database,
                  table = c("Catalogue of Life",
                            "NCBI"))) %>%
    group_by(organism_1_kingdom,
             organism_4_order) %>%
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
             .keep_all = TRUE) %>%
    select(
      -organismOriginal,-organismTranslated,-organismSanitized,-organism_7_species,-organism_6_genus,-organism_5_family
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
    arrange(match(x = organism_database,
                  table = c("Catalogue of Life",
                            "NCBI"))) %>%
    group_by(organism_1_kingdom,
             organism_3_class) %>%
    mutate(organism_2_phylum = na.locf(
      object = organism_2_phylum,
      na.rm = FALSE,
      fromLast = TRUE
    )) %>%
    ungroup() %>%
    distinct(organism_3_class,
             .keep_all = TRUE) %>%
    select(
      -organismOriginal,-organismTranslated,-organismSanitized,-organism_7_species,-organism_6_genus,-organism_5_family,-organism_4_order,
      
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
    arrange(match(x = organism_database,
                  table = c("Catalogue of Life",
                            "NCBI"))) %>%
    group_by(organism_2_phylum) %>%
    mutate(organism_1_kingdom = na.locf(
      object = organism_1_kingdom,
      na.rm = FALSE,
      fromLast = TRUE
    )) %>%
    ungroup() %>%
    distinct(organism_2_phylum,
             .keep_all = TRUE) %>%
    select(
      -organismOriginal,-organismTranslated,-organismSanitized,-organism_7_species,-organism_6_genus,-organism_5_family,-organism_4_order,-organism_3_class
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
    select(-organismOriginal,-organismTranslated)
  
  tojoin <- dfsel
  
  newdf = left_join(tojoin,
                    kingdom_tojoin,
                    by = c("organismSanitized" = "organismSanitized")) %>%
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
             .keep_all = TRUE)
  
  newdf$organism_modified_taxonomy_auto <-
    y_as_na(x = newdf$organism_modified_taxonomy_auto,
            y = "")
  
  return(newdf)
}

#######################################################
#######################################################

getref <- function(X) {
  tryCatch({
    cr_works(
      flq = c(query.bibliographic = X),
      sort = 'score',
      order = "desc",
      limit = 1
    )
  },
  error = function(e) {
    NA
  })
}

#######################################################
#######################################################

getrefPubmed <- function(X) {
  tryCatch({
    df <- entrez_summary(db = "pubmed", id = X)
    
    translatedDoi <-
      ifelse(test = "doi" %in% df[["articleids"]][, 1],
             yes = trimws(df[["articleids"]][["value"]][[which(df[["articleids"]] == "doi")]]),
             no = NA)
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
  })
}

#######################################################
#######################################################

getrefDoi <- function(X) {
  tryCatch({
    cr_works(dois = X)
  },
  error = function(e) {
    NA
  })
}

#######################################################
#######################################################
