# title: "Reference cleaner"

# loading paths
source("paths.R")

# loading functions
source("functions/reference.R")

# loading files
print(x = "loading crossref translations file, this may take a while")
dataTranslated <- read_delim(
  file = gzfile(pathDataInterimTablesTranslatedReferenceFile),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
)

print(x = "loading pmcid file, this may take a while")

PMC_ids <- read_delim(
  file = gzfile(pathDataExternalTranslationSourcePubmedFile),
  delim = ",",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  filter(!is.na(DOI) | !is.na(PMID)) %>%
  select(DOI,
         PMCID,
         PMID) %>%
  mutate_all(as.character)

### Find something appropriate
# test <- dataTranslated %>%
#   filter(!is.na(referenceTranslatedTitle)) %>%
#   group_by(referenceTranslatedTitle) %>%
#   distinct(
#     referenceOriginal_doi,
#     referenceOriginal_original,
#     referenceOriginal_publishingDetails,
#     referenceOriginal_pubmed,
#     referenceOriginal_split,
#     referenceOriginal_title,
#     .keep_all = TRUE
#   ) %>%
#   add_count() %>%
#   arrange(desc(n))

# COMMENT: JUST AS FOR BIOLOGICAL SOURCES THOSE LINES SHOULD COME THEN
## BEFORE TRANSLATION OF THE REFERENCE. WE DO IT AFTER BECAUSE YOU CAN
### NOT GUESS IT BUT THEN BUILDING A REPLACEMENT DIC SEEMS GOOD FOR ME

# test2 <- test %>%
#   filter(grepl(
#     "Harborne, The Handbook of Natural Flavonoids",
#     referenceOriginal_original
#   )) %>%
#   mutate(
#     referenceSplitNew = gsub(
#       pattern = "Harborne, The Handbook of Natural Flavonoids, 2, (1999), 1,Anthocyanins",
#       replacement = "",
#       x = referenceOriginal_original,
#       fixed = TRUE
#     ),
#     referenceSplitNew = gsub(
#       pattern = "Harborne, The Handbook of Natural Flavonoids, 1, (1999), 2.Flavones, John Wiley & Son",
#       replacement = "",
#       x = referenceSplitNew,
#       fixed = TRUE
#     ),
#     referenceSplitNew = gsub(
#       pattern = "Harborne, The Handbook of Natural Flavonoids, 1, (1999), 3.Flavone O-glycosides, John Wiley & Son",
#       replacement = "",
#       x = referenceSplitNew,
#       fixed = TRUE
#     ),
#     referenceSplitNew = gsub(
#       pattern = "Harborne, The Handbook of Natural Flavonoids, 1, (1999), 3.Flavone O-glycosides,",
#       replacement = "",
#       x = referenceSplitNew,
#       fixed = TRUE
#     ),
#     referenceSplitNew = gsub(
#       pattern = "Harborne, The Handbook of Natural Flavonoids, 1, (1999), 8.Flavone O-glycosides,",
#       replacement = "",
#       x = referenceSplitNew,
#       fixed = TRUE
#     ),
#     referenceSplitNew = gsub(
#       pattern = "Harborne, The Handbook of Natural Flavonoids, 1, (1999), 9.Flavone O-glycosides,",
#       replacement = "",
#       x = referenceSplitNew,
#       fixed = TRUE
#     ),
#     referenceSplitNew = gsub(
#       pattern = "Harborne, The Handbook of Natural Flavonoids, 1, (1999), 11.Flavone O-glycosides,",
#       replacement = "",
#       x = referenceSplitNew,
#       fixed = TRUE
#     ),
#     referenceSplitNew = gsub(
#       pattern = "Harborne, The Handbook of Natural Flavonoids, 1, (1999), 12.Flavone O-glycosides,",
#       replacement = "",
#       x = referenceSplitNew,
#       fixed = TRUE
#     ),
#     referenceSplitNew = gsub(
#       pattern = "Harborne, The Handbook of Natural Flavonoids, 1, (1999), 13.Flavone O-glycosides,",
#       replacement = "",
#       x = referenceSplitNew,
#       fixed = TRUE
#     ),
#     referenceSplitNew = gsub(
#       pattern = "Harborne, The Handbook of Natural Flavonoids, 1, (1999), 181.Flavonols",
#       replacement = "",
#       x = referenceSplitNew,
#       fixed = TRUE
#     ),
#     referenceSplitNew = gsub(
#       pattern = "Harborne, The Handbook of Natural Flavonoids, 1, (1999), 297.Flavonol O-glycosides",
#       replacement = "",
#       x = referenceSplitNew,
#       fixed = TRUE
#     ),
#     referenceSplitNew = gsub(
#       pattern = "Harborne, The Handbook of Natural Flavonoids, 2, (1999), 355,Flavans and proanthocyanidins",
#       replacement = "",
#       x = referenceSplitNew,
#       fixed = TRUE
#     ),
#     referenceSplitNew = gsub(
#       pattern = "Harborne, The Handbook of Natural Flavonoids, 2, (1999), 517,Isoflavonoids and neoflavonoids",
#       replacement = "",
#       x = referenceSplitNew,
#       fixed = TRUE
#     ),
#     referenceSplitNew = gsub(
#       pattern = "Harborne, The Handbook of Natural Flavonoids, 1, (1999), 549,C-glycosylflavones",
#       replacement = "",
#       x = referenceSplitNew,
#       fixed = TRUE
#     ),
#     referenceSplitNew = gsub(
#       pattern = "Harborne, The Handbook of Natural Flavonoids, 1, (1999), 645,Biflavonyls",
#       replacement = "",
#       x = referenceSplitNew,
#       fixed = TRUE
#     ),
#     referenceSplitNew = gsub(
#       pattern = "Harbone, Comparative Biochemisty of the Flavonoids,(1965),39, Academic Press",
#       replacement = "",
#       x = referenceSplitNew,
#       fixed = TRUE
#     ),
#     referenceSplitNew = gsub(
#       pattern = "Harborne, Comparative Biochemstry of the Flavonoides,(1967),41, Academic Press",
#       replacement = "",
#       x = referenceSplitNew,
#       fixed = TRUE
#     ),
#     referenceSplitNew = gsub(
#       pattern = "Harborne,Comparative Biochemistry of the Flavonoids, (1967) Academic Press",
#       replacement = "",
#       x = referenceSplitNew,
#       fixed = TRUE
#     ),
#     referenceSplitNew = gsub(
#       pattern = "Chalcones,dihydrochalcones and aurones",
#       replacement = "",
#       x = referenceSplitNew,
#       fixed = FALSE
#     ),
#     referenceSplitNew = gsub(
#       pattern = "Harborne, The Handbook of Natural Flavonoids, 2, (1999)",
#       replacement = "",
#       x = referenceSplitNew,
#       fixed = TRUE
#     ),
#     referenceSplitNew = gsub(
#       pattern = ", [0-9]{3},",
#       replacement = "",
#       x = referenceSplitNew,
#       fixed = FALSE
#     ),
#     referenceSplitNew = gsub(
#       pattern = "Harborne, The Handbook of Natural Flavonoids, 2, (1999)Flavanones and dihydroflavonols",
#       replacement = "",
#       x = referenceSplitNew,
#       fixed = TRUE
#     )
#   )

# test3 <- test %>%
#   filter(n == 631) %>%
#   mutate(
#     referenceSplitNew = gsub(
#       pattern = "; Khimiya .*",
#       replacement = "",
#       x = referenceOriginal_original,
#       fixed = FALSE
#     )
#   )

# test4 <- test %>%
#   filter(
#     n != 631 &
#       !grepl(
#         "Harborne, The Handbook of Natural Flavonoids",
#         referenceOriginal_original
#       ) &
#       n > 10
#   )

# test5 <- test4 %>%
#   filter(n == 274) # problematic issue from PlantaMed
#
# test6 <- test4 %>%
#   filter(n == 194) # problematic too... DNP

## example
### Luesch, Hendrik; Yoshida, Wesley Y.; Moore, Richard E.; Paul, Valerie J.; Journal of Natural Products; vol. 63; 10; (2000); p. 1437 - 1439.
#### returning erroneous result
##### Journal of Natural Products; vol. 63; 10; (2000); p. 1437 - 1439
###### test <- cr_works(query = "Journal of Natural Products; vol. 63; 10; (2000); p. 1437 - 1439") #returning right result (but very low score)

print(x = "cleaning, this may take a while if running full mode")

dataCleaned <- dataTranslated %>%
  mutate(referenceCleanedValue = referenceTranslatedValue,
         referenceCleanedType = referenceTranslatedType)

dataCleanedScore <- dataCleaned %>%
  filter(referenceCleanedType == "title") %>%
  filter(!is.na(organismOriginal) &
           !is.na(referenceCleanedValue)) %>%
  distinct(organismOriginal,
           organismCleaned,
           referenceCleanedValue,
           level) %>%
  rowwise() %>%
  mutate(referenceCleaned_scoreTitleOrganism = ifelse(
    test = str_detect(
      string = fixed(referenceCleanedValue),
      pattern = paste(word(organismOriginal, 1),
                      word(organismCleaned, 1),
                      sep = "|")
    ),
    yes = 1,
    no = 0
  )) %>%
  filter(referenceCleaned_scoreTitleOrganism == 1)

dataCleanedJoined <- left_join(dataCleaned, dataCleanedScore)

dataCleanedJoinedWide <- dataCleanedJoined %>%
  pivot_wider(names_from = referenceCleanedType,
              names_prefix = "referenceCleaned_",
              values_from = referenceCleanedValue) %>%
  select(-referenceCleaned_NA) %>%
  distinct(
    organismOriginal,
    level,
    referenceValue,
    referenceCleaned_doi,
    referenceCleaned_journal,
    referenceCleaned_title,
    referenceCleaned_date,
    referenceCleaned_author,
    referenceCleaned_scoreCrossref,
    referenceCleaned_scoreDistance,
    referenceCleaned_scoreTitleOrganism,
    .keep_all = TRUE
  ) %>%
  group_by(organismOriginal, referenceValue) %>%
  arrange(desc(referenceCleaned_scoreCrossref)) %>%
  arrange(referenceCleaned_scoreDistance) %>%
  arrange(desc(referenceCleaned_scoreTitleOrganism)) %>%
  ungroup() %>%
  distinct(organismOriginal,
           referenceValue,
           referenceTranslatedType,
           .keep_all = TRUE) %>%
  select(-organismOriginal, -organismCleaned) %>%
  mutate_all(as.character)

dataCleanedJoinedLong <- dataCleanedJoinedWide %>%
  pivot_longer(
    cols = 6:ncol(.),
    names_to = c("drop", "referenceCleanedType"),
    names_sep = "_",
    values_to = "referenceCleanedValue",
    values_drop_na = TRUE
  ) %>%
  distinct(referenceType,
           referenceValue,
           level,
           referenceCleanedType,
           referenceCleanedValue)

dataCleanedJoinedLongWide <- dataCleanedJoinedLong %>%
  pivot_wider(names_from = referenceCleanedType,
              names_prefix = "referenceCleaned_",
              values_from = referenceCleanedValue)

subDataClean_doi <- dataCleanedJoinedLongWide %>%
  filter(!is.na(referenceCleaned_doi)) %>%
  distinct(referenceCleaned_doi)

subDataClean_pmid <- dataCleaned %>%
  filter(referenceType == "pubmed") %>%
  distinct(referenceValue)

df_doi <- left_join(subDataClean_doi,
                    PMC_ids,
                    by = c("referenceCleaned_doi" = "DOI")) %>%
  filter(!is.na(PMID) | !is.na(PMCID)) %>%
  select(
    referenceCleaned_doi,
    referenceCleaned_pmid = PMID,
    referenceCleaned_pmcid = PMCID
  )

df_pubmed <- left_join(subDataClean_pmid,
                       PMC_ids,
                       by = c("referenceValue" = "PMID")) %>%
  filter(!is.na(referenceValue) | !is.na(PMCID)) %>%
  select(referenceValue,
         referenceCleaned_pmcid = PMCID) %>%
  mutate(referenceCleaned_pmid = referenceValue)

tableJoined <- left_join(dataCleanedJoinedLongWide, df_doi)

referenceTable <-
  left_join(tableJoined,
            df_pubmed,
            by = c("referenceValue" = "referenceValue")) %>%
  mutate(
    referenceCleaned_pmid = ifelse(
      test = !is.na(referenceCleaned_pmid.x),
      yes = referenceCleaned_pmid.x,
      no = referenceCleaned_pmid.y
    ),
    referenceCleaned_pmcid = ifelse(
      test = !is.na(referenceCleaned_pmcid.x),
      yes = referenceCleaned_pmcid.x,
      no = referenceCleaned_pmcid.y
    )
  ) %>%
  select(
    referenceType,
    referenceValue,
    referenceCleanedDoi = referenceCleaned_doi,
    referenceCleanedPmcid = referenceCleaned_pmcid,
    referenceCleanedPmid = referenceCleaned_pmid,
    referenceCleaned_title,
    referenceCleaned_journal,
    referenceCleaned_date,
    referenceCleaned_author,
    referenceCleaned_score_crossref = referenceCleaned_scoreCrossref,
    referenceCleaned_score_distance = referenceCleaned_scoreDistance,
    referenceCleaned_score_titleOrganism = referenceCleaned_scoreTitleOrganism,
  ) %>%
  distinct(
    referenceType,
    referenceValue,
    referenceCleanedDoi,
    referenceCleanedPmcid,
    referenceCleanedPmid,
    referenceCleaned_title,
    referenceCleaned_journal,
    referenceCleaned_date,
    referenceCleaned_author,
    referenceCleaned_score_crossref,
    referenceCleaned_score_distance,
    referenceCleaned_score_titleOrganism,
  )

# exporting
## creating directories if they do not exist
ifelse(
  !dir.exists(pathDataInterimTablesCleaned),
  dir.create(pathDataInterimTablesCleaned),
  FALSE
)

ifelse(
  !dir.exists(pathDataInterimTablesCleanedReference),
  dir.create(pathDataInterimTablesCleanedReference),
  FALSE
)

write.table(
  x = referenceTable,
  file = gzfile(
    description = pathDataInterimTablesCleanedReferenceFile,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
