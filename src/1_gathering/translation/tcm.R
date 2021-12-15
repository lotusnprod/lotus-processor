# title: "TCM to scientific name translatoR"

# loading
## paths
source("paths.R")

library(dplyr)
library(readr)
library(readxl)
library(splitstackshape)
library(stringi)
library(stringr)

## functions
source("r/capitalize.R")
source("r/tcm_cleaning.R")
source("r/tcm_inverting.R")
source("r/y_as_na.R")

if (mode == "full") {
  ## files
  ### dictionary from TMMC
  database <- databases$get("tmmc")
  tcmNamesDic_1 <- read_excel(database$sourceFiles$tsv,
    sheet = 1
  ) %>%
    mutate_all(as.character) %>%
    select(
      latin = LATIN,
      common = COMMON,
      biologicalsource = SCIENCE
    ) %>%
    distinct(biologicalsource, .keep_all = TRUE)

  ### dictionary from TCMID
  if (file.exists(pathDataExternalTranslationSourceTcmTcmid)) {
    tcmNamesDic_2 <-
      read_delim(file = pathDataExternalTranslationSourceTcmTcmid) %>%
      select(
        latin = `Latin Name`,
        common = `English Name`
      ) %>%
      filter(!is.na(common) &
        !is.na(latin)) %>%
      mutate(biologicalsource = latin)
  }

  ## dictionary from Chinese Medicine Board of Australia
  tcmNamesDic_3 <-
    read_excel(
      path = pathDataExternalTranslationSourceTcmCmba,
      sheet = 1
    ) %>%
    mutate_all(as.character) %>%
    select(
      latin = 6,
      common = 9
    ) %>%
    mutate(biologicalsource = latin) %>%
    filter(row.names(.) != 1) %>%
    filter(common != "N/A" &
      latin != "N/A")

  # latin genitive dictionaries
  ## i
  latinGenitiveIDic <-
    read_delim(
      file = pathDataInterimDictionariesLatinGenitiveI,
      delim = "\t",
    ) %>%
    mutate(n = str_count(string = genitive)) %>%
    arrange(desc(n)) %>%
    select(-n)

  ## is
  latinGenitiveIsDic <-
    read_delim(
      file = pathDataInterimDictionariesLatinGenitiveIs,
      delim = "\t",
    ) %>%
    mutate(n = str_count(string = genitive)) %>%
    arrange(desc(n)) %>%
    select(-n)

  ## parts
  latinGenitivePartsDic <-
    read_delim(
      file = pathDataInterimDictionariesLatinPlantParts,
      delim = "\t",
    )

  ## manually subtracted entries
  manualSubtraction <-
    read_delim(
      file = pathDataInterimDictionariesTcmManualSubtraction,
      delim = "\t",
    )

  ## forcing organisms in test file to be in sampled dic
  test_organisms <- read_delim(
    file = "../tests/tests_min.tsv",
    delim = "\t"
  )

  # cleaning
  ## tcm_1
  tcmNamesDic_1$biologicalsource <- gsub(
    pattern = "^*<i>\\s*|</i>.*",
    replacement = "",
    x = tcmNamesDic_1$biologicalsource
  )

  tcmNamesDic_1$biologicalsource <- gsub(
    pattern = "].*",
    replacement = "",
    x = tcmNamesDic_1$biologicalsource
  )

  tcmNamesDic_1$biologicalsource <- gsub(
    pattern = ".*\\[<i>",
    replacement = "",
    x = tcmNamesDic_1$biologicalsource
  )

  ## tcm_3
  tcmNamesDic_3 <- tcmNamesDic_3 %>%
    cSplit(
      splitCols = "biologicalsource",
      sep = "1. ",
      fixed = TRUE,
      stripWhite = FALSE
    ) %>%
    select(common,
      biologicalsource = biologicalsource_2
    ) %>%
    filter(!is.na(biologicalsource)) %>%
    cSplit(
      splitCols = "biologicalsource",
      sep = "2. ",
      fixed = TRUE,
      stripWhite = FALSE
    ) %>%
    cSplit(
      splitCols = "biologicalsource_2",
      sep = "3. ",
      fixed = TRUE,
      stripWhite = FALSE
    ) %>%
    cSplit(
      splitCols = "biologicalsource_2_2",
      sep = "4. ",
      fixed = TRUE,
      stripWhite = FALSE
    ) %>%
    cSplit(
      splitCols = "biologicalsource_2_2_2",
      sep = "5. ",
      fixed = TRUE,
      stripWhite = FALSE
    ) %>%
    cSplit(
      splitCols = "biologicalsource_2_2_2_2",
      sep = "6. ",
      fixed = TRUE,
      stripWhite = FALSE
    ) %>%
    cSplit(
      splitCols = "biologicalsource_2_2_2_2_2",
      sep = "7. ",
      fixed = TRUE,
      stripWhite = FALSE
    ) %>%
    group_by(rownames(.)) %>%
    summarise(
      common = first(common),
      biologicalsource = paste(
        biologicalsource_1[!is.na(biologicalsource_1)],
        biologicalsource_2_1[!is.na(biologicalsource_2_1)],
        biologicalsource_2_2_1[!is.na(biologicalsource_2_2_1)],
        biologicalsource_2_2_2_1[!is.na(biologicalsource_2_2_2_1)],
        biologicalsource_2_2_2_2_1[!is.na(biologicalsource_2_2_2_2_1)],
        biologicalsource_2_2_2_2_2_1[!is.na(biologicalsource_2_2_2_2_2_1)],
        biologicalsource_2_2_2_2_2_2[!is.na(biologicalsource_2_2_2_2_2_2)],
        sep = " "
      )
    ) %>%
    select(common, biologicalsource)

  # capitalizing first letter only
  tcmNamesDic_1$latin <- tolower(x = tcmNamesDic_1$latin)
  tcmNamesDic_1$latin <- capitalize(string = tcmNamesDic_1$latin)

  tcmNamesDic_1$common <- tolower(x = tcmNamesDic_1$common)
  tcmNamesDic_1$common <- capitalize(string = tcmNamesDic_1$common)

  if (file.exists(pathDataExternalTranslationSourceTcmTcmid)) {
    tcmNamesDic_2$latin <- tolower(x = tcmNamesDic_2$latin)
    tcmNamesDic_2$latin <- capitalize(string = tcmNamesDic_2$latin)

    tcmNamesDic_2$common <- tolower(x = tcmNamesDic_2$common)
    tcmNamesDic_2$common <-
      capitalize(string = tcmNamesDic_2$common)

    # reordering tcm names (eg. not Radix gentianae but Gentianae radix)
    tcmNamesDic_2 <-
      tcm_inverting(x = tcmNamesDic_2) %>%
      mutate(latin = biologicalsource) %>%
      select(
        latin,
        common,
        biologicalsource
      ) %>%
      filter(!is.na(biologicalsource)) %>%
      distinct(biologicalsource, .keep_all = TRUE)

    tcmNamesDic_2 <-
      tcm_cleaning(x = tcmNamesDic_2) %>%
      filter(!is.na(biologicalsource)) %>%
      distinct(biologicalsource, .keep_all = TRUE) %>%
      mutate(biologicalsource = latin) %>%
      select(common, biologicalsource)
  }

  # cleaning most occurring biological parts from names (eg. radix, flos, folium, cortrex, ...)
  tcmNamesDic_1 <-
    tcm_cleaning(x = tcmNamesDic_1) %>%
    filter(!is.na(biologicalsource)) %>%
    distinct(biologicalsource, .keep_all = TRUE)

  tcmNamesDic_1_1 <- tcmNamesDic_1 %>%
    select(common = latin, biologicalsource) %>%
    filter(!is.na(common))

  tcmNamesDic_1_2 <- tcmNamesDic_1 %>%
    select(common, biologicalsource) %>%
    filter(!is.na(common))

  # joining
  tcmNamesDic <-
    rbind(
      tcmNamesDic_1_1,
      tcmNamesDic_1_2,
      if (file.exists(pathDataExternalTranslationSourceTcmTcmid)) {
        tcmNamesDic_2
      },
      tcmNamesDic_3
    ) %>%
    distinct(common, .keep_all = TRUE) %>%
    cSplit(
      splitCols = "biologicalsource",
      sep = " "
    ) %>%
    select(common,
      biologicalsource = biologicalsource_01
    )

  # filtering empty entries
  tcmNamesDic <- tcmNamesDic %>%
    filter(biologicalsource != " NA ")

  # manual curation
  ## ae
  tcmNamesDic_ae <- tcmNamesDic %>%
    distinct(biologicalsource, .keep_all = TRUE) %>%
    arrange(biologicalsource) %>%
    filter(grepl(
      pattern = "\\w+ae\\b",
      x = biologicalsource
    )) %>%
    mutate(newbiologicalsource = biologicalsource)

  tcmNamesDic_ae$newbiologicalsource <- gsub(
    pattern = ".{2}$",
    replacement = "a",
    x = tcmNamesDic_ae$newbiologicalsource
  )

  ## i
  tcmNamesDic_i <- tcmNamesDic %>%
    distinct(biologicalsource, .keep_all = TRUE) %>%
    arrange(biologicalsource) %>%
    filter(grepl(
      pattern = "\\w+i\\b",
      x = biologicalsource
    )) %>%
    mutate(newbiologicalsource = biologicalsource)

  a <- paste0("\\b", latinGenitiveIDic$genitive, "\\b")
  b <- latinGenitiveIDic$nominative

  tcmNamesDic_i$newbiologicalsource <-
    stri_replace_all_regex(
      str = tcmNamesDic_i$newbiologicalsource,
      pattern = a,
      replacement = b,
      case_insensitive = FALSE,
      vectorize_all = FALSE
    )
  # example Alhagi
  # example Litchi
  # example Muscari
  # example rice bean -> dangerous
  # example Seseli
  # maybe example Strychni Strychnos
  # example Rose apple beautiful Syzygii Syzygium
  # example Tadehagi or Thlaspi

  ## is
  tcmNamesDic_is <- tcmNamesDic %>%
    distinct(biologicalsource, .keep_all = TRUE) %>%
    arrange(biologicalsource) %>%
    filter(grepl(
      pattern = "\\w+is\\b",
      x = biologicalsource
    )) %>%
    mutate(newbiologicalsource = biologicalsource)

  c <- paste0("\\b", latinGenitiveIsDic$genitive, "\\b")
  d <- latinGenitiveIsDic$nominative

  tcmNamesDic_is$newbiologicalsource <-
    stri_replace_all_regex(
      str = tcmNamesDic_is$newbiologicalsource,
      pattern = c,
      replacement = d,
      case_insensitive = FALSE,
      vectorize_all = FALSE
    )

  # example Adonis
  # example Agathis, Aleuritopteris, Amaryllis, ...

  # joining
  tcmNamesDicCurated <- full_join(tcmNamesDic, tcmNamesDic_ae)

  # removing some organism parts
  e <- paste0("\\b", latinGenitivePartsDic$part, "\\b")
  f <- as.character(latinGenitivePartsDic$replacement)

  tcmNamesDicCurated$newbiologicalsource <-
    stri_replace_all_regex(
      str = tcmNamesDicCurated$newbiologicalsource,
      pattern = e,
      replacement = f,
      case_insensitive = FALSE,
      vectorize_all = FALSE
    )

  tcmNamesDicCurated$newbiologicalsource <-
    y_as_na(
      x = tcmNamesDicCurated$newbiologicalsource,
      y = "FALSE"
    )

  # joining
  tcmNamesDicCurated <-
    rbind(tcmNamesDicCurated, tcmNamesDic_i, tcmNamesDic_is) %>%
    arrange(newbiologicalsource) %>%
    distinct(common, .keep_all = TRUE)

  # deleting ambiguous entries
  tcmNamesDicCurated <- tcmNamesDicCurated %>%
    filter(!common %in% manualSubtraction$common)

  # sorting in appropriate order
  tcmNamesDicCurated <- tcmNamesDicCurated %>%
    mutate(n = str_count(string = common)) %>%
    arrange(desc(n)) %>%
    # sorting for replacements like "sea cucumber" and so on...
    select(
      vernacularName = common,
      canonicalName = biologicalsource,
      newCanonicalName = newbiologicalsource
    ) %>%
    distinct() %>%
    filter(vernacularName != canonicalName |
      is.na(newCanonicalName))

  ## sampling rows for test mode
  "%ni%" <- Negate("%in%")
  set.seed(
    seed = 42,
    kind = "Mersenne-Twister",
    normal.kind = "Inversion"
  )
  tcmNamesDicCurated_sampled_1 <- tcmNamesDicCurated %>%
    sample_n(500)

  tcmNamesDicCurated_sampled_2 <- tcmNamesDicCurated %>%
    filter(grepl(
      pattern = paste(test_organisms$organismValue, collapse = "|"),
      x = vernacularName,
      ignore.case = TRUE
    ))

  tcmNamesDicCurated_sampled <-
    bind_rows(
      tcmNamesDicCurated_sampled_1,
      tcmNamesDicCurated_sampled_2
    ) %>%
    mutate(n = str_count(string = vernacularName)) %>%
    arrange(desc(n)) %>%
    select(-n)

  # exporting
  write_delim(
    x = tcmNamesDicCurated,
    delim = "\t",
    file = gzfile(
      description = pathDataInterimDictionariesTcmNames,
      compression = 9,
      encoding = "UTF-8"
    ),
    na = "",
    quote = "none",
    escape = "double"
  )

  write_delim(
    x = tcmNamesDicCurated_sampled,
    delim = "\t",
    file = gzfile(
      description = pathTestsDicTcmFile,
      compression = 9,
      encoding = "UTF-8"
    ),
    na = "",
    quote = "none",
    escape = "double"
  )
  ## because of univocity parser settings
}
