cat("This script checks for organism presence in cleaned reference title \n")

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")

cat("... functions \n")
source("functions/reference.R")
source("functions/helpers.R")

cat("loading crossref translations file, this may take a while \n")
dataTranslated <- read_delim(
  file = gzfile(pathDataInterimTablesTranslatedReferenceFile),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  tibble()

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

cat("cleaning \n")
dataCleaned <- dataTranslated %>%
  filter(!is.na(referenceTranslatedType)) %>%
  filter(!is.na(referenceTranslatedValue)) %>%
  mutate(
    referenceCleanedValue = referenceTranslatedValue,
    referenceCleanedType = referenceTranslatedType
  ) %>%
  select(
    -referenceTranslatedValue,
    -referenceTranslatedType
  )

rm(dataTranslated)

cat("checking for organism in title, may take a while if running full mode \n")
dataCleanedScore <- dataCleaned %>%
  filter(referenceCleanedType == "title") %>%
  filter(!is.na(organismOriginal) &
    !is.na(referenceCleanedValue)) %>%
  mutate(referenceCleaned_scoreTitleOrganism = ifelse(
    test = str_detect(
      string = fixed(tolower(
        referenceCleanedValue
      )),
      pattern = fixed(tolower(
        word(organismCleaned, 1)
      ))
    ),
    yes = 1,
    no = 0
  )) %>%
  filter(referenceCleaned_scoreTitleOrganism == 1) %>%
  pivot_wider(
    names_from = referenceCleanedType,
    names_prefix = "referenceCleaned_",
    values_from = referenceCleanedValue
  ) %>%
  mutate_all(as.character) %>%
  pivot_longer(
    cols = (ncol(.) - 1):ncol(.),
    names_to = c("drop", "referenceCleanedType"),
    names_sep = "_",
    values_to = "referenceCleanedValue",
    values_drop_na = TRUE
  ) %>%
  filter(referenceCleanedType == "scoreTitleOrganism") %>%
  select(-drop) %>%
  tibble() %>%
  distinct()

dataCleanedJoined <- bind_rows(dataCleaned, dataCleanedScore) %>%
  tibble() %>%
  filter(!is.na(organismOriginal))

rm(dataCleaned)

# this is because sadly crossref does not always give the same score, therefore
## we do not have unique values ...
subDataCleanedJoined_1 <- dataCleanedJoined %>%
  filter(referenceCleanedType == "scoreCrossref") %>%
  distinct(organismOriginal,
    referenceType,
    referenceValue,
    organismCleaned,
    level,
    .keep_all = TRUE
  )

# this is because sadly crossref does not always give the same DOI, therefore
## we do not have unique values ...
subDataCleanedJoined_2 <- dataCleanedJoined %>%
  filter(referenceCleanedType != "scoreCrossref") %>%
  distinct(
    organismOriginal,
    referenceType,
    referenceValue,
    organismCleaned,
    level,
    referenceCleanedType,
    .keep_all = TRUE
  )

dataCleanedJoinedUnique <-
  bind_rows(
    subDataCleanedJoined_1,
    subDataCleanedJoined_2
  )

rm(
  dataCleanedJoined,
  subDataCleanedJoined_1,
  subDataCleanedJoined_2
)

gc()

cat("manipulating and keeping best result only (long step) \n")
dataCleanedJoinedWide <- dataCleanedJoinedUnique %>%
  pivot_wider(
    names_from = referenceCleanedType,
    names_prefix = "referenceCleaned_",
    values_from = referenceCleanedValue
  )

rm(dataCleanedJoinedUnique)

dataCleanedJoinedWide_1 <- dataCleanedJoinedWide %>%
  filter(referenceType == "doi" |
    referenceType == "pubmed" |
    referenceType == "title") %>%
  group_by(organismOriginal, organismCleaned, referenceValue) %>%
  arrange(desc(as.numeric(referenceCleaned_scoreCrossref))) %>%
  arrange(desc(as.numeric(referenceCleaned_scoreTitleOrganism))) %>%
  arrange(as.numeric(referenceCleaned_scoreDistance)) %>%
  ungroup() %>%
  distinct(organismOriginal,
    organismCleaned,
    referenceValue,
    .keep_all = TRUE
  ) %>%
  select(-level) %>%
  mutate(
    referenceCleaned_scoreComplement_date = NA,
    referenceCleaned_scoreComplement_author = NA,
    referenceCleaned_scoreComplement_journal = NA,
    referenceCleaned_scoreComplement_total = NA
  )

dataCleanedJoinedWide_2 <- dataCleanedJoinedWide %>%
  filter(
    referenceType == "original" |
      referenceType == "publishingDetails" |
      referenceType == "split"
  ) %>%
  mutate(
    referenceCleaned_scoreComplement_date = ifelse(
      str_detect(
        string = referenceValue,
        pattern = substr(x = referenceCleaned_date, start = 1, stop = 4)
      ),
      yes = 1,
      no = 0
    ),
    referenceCleaned_scoreComplement_author = ifelse(
      str_detect(
        string = tolower(referenceValue),
        pattern = fixed(tolower(word(
          referenceCleaned_author, 1
        )))
      ),
      yes = 1,
      no = 0
    ),
    referenceCleaned_scoreComplement_journal = ifelse(
      str_detect(
        string = tolower(referenceValue),
        pattern = fixed(tolower(referenceCleaned_journal))
      ),
      yes = 1,
      no = 0
    )
  ) %>%
  mutate(
    referenceCleaned_scoreComplement_total =
      referenceCleaned_scoreComplement_date +
        referenceCleaned_scoreComplement_author +
        referenceCleaned_scoreComplement_journal
  ) %>%
  group_by(organismOriginal, organismCleaned, referenceValue) %>%
  arrange(desc(as.numeric(referenceCleaned_scoreCrossref))) %>%
  arrange(desc(as.numeric(referenceCleaned_scoreComplement_total))) %>%
  arrange(desc(as.numeric(referenceCleaned_scoreTitleOrganism))) %>%
  ungroup() %>%
  distinct(organismOriginal,
    organismCleaned,
    referenceValue,
    .keep_all = TRUE
  ) %>%
  select(-level)

rm(dataCleanedJoinedWide)

dataCleanedJoinedWideScore <- bind_rows(
  dataCleanedJoinedWide_1,
  dataCleanedJoinedWide_2
) %>%
  mutate(referenceCleaned_doi = toupper(referenceCleaned_doi))

rm(
  dataCleanedJoinedWide_1,
  dataCleanedJoinedWide_2
)

subDataClean_doi <- dataCleanedJoinedWideScore %>%
  filter(!is.na(referenceCleaned_doi)) %>%
  distinct(referenceCleaned_doi) %>%
  mutate_all(as.character)

subDataClean_pmid <- dataCleanedJoinedWideScore %>%
  filter(referenceType == "pubmed") %>%
  distinct(referenceValue) %>%
  mutate_all(as.character)

cat("loading pmcid file, this may take a while \n")
# here because of memory
PMC_ids <- read_delim(
  file = gzfile(pathDataExternalTranslationSourcePubmedFile),
  delim = ",",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  filter(!is.na(DOI) | !is.na(PMID)) %>%
  select(
    DOI,
    PMCID,
    PMID
  ) %>%
  mutate(DOI = toupper(DOI)) %>%
  mutate_all(as.character) %>%
  tibble()

cat("adding PMID and PMCID \n")
df_doi <- left_join(subDataClean_doi,
  PMC_ids,
  by = c("referenceCleaned_doi" = "DOI")
) %>%
  filter(!is.na(PMID) | !is.na(PMCID)) %>%
  select(
    referenceCleaned_doi,
    referenceCleaned_pmid = PMID,
    referenceCleaned_pmcid = PMCID
  )

df_pubmed <- left_join(subDataClean_pmid,
  PMC_ids,
  by = c("referenceValue" = "PMID")
) %>%
  filter(!is.na(referenceValue) | !is.na(PMCID)) %>%
  select(referenceValue,
    referenceCleaned_pmcid = PMCID
  ) %>%
  mutate(referenceCleaned_pmid = referenceValue)

tableJoined <- left_join(dataCleanedJoinedWideScore, df_doi)

referenceTable <-
  left_join(tableJoined,
    df_pubmed,
    by = c("referenceValue" = "referenceValue")
  ) %>%
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
    organismOriginal,
    organismCleaned,
    referenceType,
    referenceValue,
    referenceCleanedDoi = referenceCleaned_doi,
    referenceCleanedPmcid = referenceCleaned_pmcid,
    referenceCleanedPmid = referenceCleaned_pmid,
    referenceCleanedTitle = referenceCleaned_title,
    referenceCleaned_journal,
    referenceCleaned_date,
    referenceCleaned_author,
    referenceCleaned_score_crossref = referenceCleaned_scoreCrossref,
    referenceCleaned_score_distance = referenceCleaned_scoreDistance,
    referenceCleaned_score_titleOrganism = referenceCleaned_scoreTitleOrganism,
    referenceCleaned_score_complementDate = referenceCleaned_scoreComplement_date,
    referenceCleaned_score_complementAuthor = referenceCleaned_scoreComplement_author,
    referenceCleaned_score_complementJournal =  referenceCleaned_scoreComplement_journal,
    referenceCleaned_score_complementTotal = referenceCleaned_scoreComplement_total
  ) %>%
  distinct() %>%
  mutate(across(everything(), ~ y_as_na(.x, "NULL")))

cat("ensuring directories exist \n")
ifelse(
  test = !dir.exists(pathDataInterimTablesCleaned),
  yes = dir.create(pathDataInterimTablesCleaned),
  no = paste(pathDataInterimTablesCleaned, "exists")
)

ifelse(
  test = !dir.exists(pathDataInterimTablesCleanedReference),
  yes = dir.create(pathDataInterimTablesCleanedReference),
  no = paste(pathDataInterimTablesCleanedReference, "exists")
)

cat("exporting ... \n")
cat(pathDataInterimTablesCleanedReferenceFile, "\n")
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

end <- Sys.time()

cat("Script finished in", format(end - start), "\n")