#######################################################
####################   Functions   ####################
#######################################################

library(Hmisc)
library(data.table)
library(jsonlite)
library(parallel)
library(pbmcapply)
library(splitstackshape)
library(stringi)
# library(zoo)
library(tidyverse)

#######################################################
#######################################################

cut <- 10000

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
#######################################################

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
#######################################################

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
#######################################################

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
#######################################################