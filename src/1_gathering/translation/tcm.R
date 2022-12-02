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
  tcmNamesDic_1 <-
    readxl::read_excel(
      path = database$sourceFiles$tsv,
      sheet = 1
    ) |>
    dplyr::mutate_all(as.character) |>
    dplyr::select(
      latin = LATIN,
      common = COMMON,
      biologicalsource = SCIENCE
    ) |>
    dplyr::distinct(biologicalsource, .keep_all = TRUE)

  ### dictionary from TCMID
  if (file.exists(pathDataExternalTranslationSourceTcmTcmid)) {
    tcmNamesDic_2 <-
      readr::read_delim(file = pathDataExternalTranslationSourceTcmTcmid) |>
      dplyr::select(
        latin = `Latin Name`,
        common = `English Name`
      ) |>
      dplyr::filter(!is.na(common) &
        !is.na(latin)) |>
      dplyr::mutate(biologicalsource = latin)
  }

  ## dictionary from Chinese Medicine Board of Australia
  tcmNamesDic_3 <-
    readxl::read_excel(
      path = pathDataExternalTranslationSourceTcmCmba,
      sheet = 1
    ) %>%
    dplyr::mutate_all(as.character) %>%
    dplyr::select(
      latin = 6,
      common = 9
    ) %>%
    dplyr::mutate(biologicalsource = latin) %>%
    dplyr::filter(row.names(.) != 1) %>%
    dplyr::filter(common != "N/A" &
      latin != "N/A")

  # latin genitive dictionaries
  ## i
  latinGenitiveIDic <-
    readr::read_delim(
      file = pathDataInterimDictionariesLatinGenitiveI,
      delim = "\t"
    ) |>
    dplyr::mutate(n = stringr::str_count(string = genitive)) |>
    dplyr::arrange(dplyr::desc(n)) |>
    dplyr::select(-n)

  ## is
  latinGenitiveIsDic <-
    readr::read_delim(
      file = pathDataInterimDictionariesLatinGenitiveIs,
      delim = "\t"
    ) |>
    dplyr::mutate(n = stringr::str_count(string = genitive)) |>
    dplyr::arrange(dplyr::desc(n)) %>%
    dplyr::select(-n)

  ## parts
  latinGenitivePartsDic <-
    readr::read_delim(
      file = pathDataInterimDictionariesLatinPlantParts,
      delim = "\t"
    )

  ## manually subtracted entries
  manualSubtraction <-
    readr::read_delim(
      file = pathDataInterimDictionariesTcmManualSubtraction,
      delim = "\t"
    )

  ## forcing organisms in test file to be in sampled dic
  test_organisms <-
    readr::read_delim(
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
    splitstackshape::cSplit(
      splitCols = "biologicalsource",
      sep = "1. ",
      fixed = TRUE,
      stripWhite = FALSE
    ) %>%
    dplyr::select(common,
      biologicalsource = biologicalsource_2
    ) %>%
    dplyr::filter(!is.na(biologicalsource)) %>%
    splitstackshape::cSplit(
      splitCols = "biologicalsource",
      sep = "2. ",
      fixed = TRUE,
      stripWhite = FALSE
    ) %>%
    splitstackshape::cSplit(
      splitCols = "biologicalsource_2",
      sep = "3. ",
      fixed = TRUE,
      stripWhite = FALSE
    ) %>%
    splitstackshape::cSplit(
      splitCols = "biologicalsource_2_2",
      sep = "4. ",
      fixed = TRUE,
      stripWhite = FALSE
    ) %>%
    splitstackshape::cSplit(
      splitCols = "biologicalsource_2_2_2",
      sep = "5. ",
      fixed = TRUE,
      stripWhite = FALSE
    ) %>%
    splitstackshape::cSplit(
      splitCols = "biologicalsource_2_2_2_2",
      sep = "6. ",
      fixed = TRUE,
      stripWhite = FALSE
    ) %>%
    splitstackshape::cSplit(
      splitCols = "biologicalsource_2_2_2_2_2",
      sep = "7. ",
      fixed = TRUE,
      stripWhite = FALSE
    ) %>%
    dplyr::group_by(rownames(.)) %>%
    dplyr::summarise(
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
    ) |>
    dplyr::select(common, biologicalsource)

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
      tcm_inverting(x = tcmNamesDic_2) |>
      dplyr::mutate(latin = biologicalsource) |>
      dplyr::select(
        latin,
        common,
        biologicalsource
      ) |>
      dplyr::filter(!is.na(biologicalsource)) |>
      dplyr::distinct(biologicalsource, .keep_all = TRUE)

    tcmNamesDic_2 <-
      tcm_cleaning(x = tcmNamesDic_2) |>
      dplyr::filter(!is.na(biologicalsource)) |>
      dplyr::distinct(biologicalsource, .keep_all = TRUE) |>
      dplyr::mutate(biologicalsource = latin) |>
      dplyr::select(common, biologicalsource)
  }

  # cleaning most occurring biological parts from names (eg. radix, flos, folium, cortrex, ...)
  tcmNamesDic_1 <-
    tcm_cleaning(x = tcmNamesDic_1) |>
    dplyr::filter(!is.na(biologicalsource)) |>
    dplyr::distinct(biologicalsource, .keep_all = TRUE)

  tcmNamesDic_1_1 <- tcmNamesDic_1 |>
    dplyr::select(common = latin, biologicalsource) |>
    dplyr::filter(!is.na(common))

  tcmNamesDic_1_2 <- tcmNamesDic_1 |>
    dplyr::select(common, biologicalsource) |>
    dplyr::filter(!is.na(common))

  # joining
  tcmNamesDic <-
    rbind(
      tcmNamesDic_1_1,
      tcmNamesDic_1_2,
      if (file.exists(pathDataExternalTranslationSourceTcmTcmid)) {
        tcmNamesDic_2
      },
      tcmNamesDic_3
    ) |>
    dplyr::distinct(common, .keep_all = TRUE) |>
    splitstackshape::cSplit(
      splitCols = "biologicalsource",
      sep = " "
    ) |>
    dplyr::select(common,
      biologicalsource = biologicalsource_01
    )

  # filtering empty entries
  tcmNamesDic <- tcmNamesDic |>
    dplyr::filter(biologicalsource != " NA ")

  # manual curation
  ## ae
  tcmNamesDic_ae <- tcmNamesDic |>
    dplyr::distinct(biologicalsource, .keep_all = TRUE) |>
    dplyr::arrange(biologicalsource) |>
    dplyr::filter(grepl(
      pattern = "\\w+ae\\b",
      x = biologicalsource
    )) |>
    dplyr::mutate(newbiologicalsource = biologicalsource)

  tcmNamesDic_ae$newbiologicalsource <- gsub(
    pattern = ".{2}$",
    replacement = "a",
    x = tcmNamesDic_ae$newbiologicalsource
  )

  ## i
  tcmNamesDic_i <- tcmNamesDic |>
    dplyr::distinct(biologicalsource, .keep_all = TRUE) |>
    dplyr::arrange(biologicalsource) |>
    dplyr::filter(grepl(
      pattern = "\\w+i\\b",
      x = biologicalsource
    )) |>
    dplyr::mutate(newbiologicalsource = biologicalsource)

  a <- paste0("\\b", latinGenitiveIDic$genitive, "\\b")
  b <- latinGenitiveIDic$nominative

  tcmNamesDic_i$newbiologicalsource <-
    stringi::stri_replace_all_regex(
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
  tcmNamesDic_is <- tcmNamesDic |>
    dplyr::distinct(biologicalsource, .keep_all = TRUE) |>
    dplyr::arrange(biologicalsource) |>
    dplyr::filter(grepl(
      pattern = "\\w+is\\b",
      x = biologicalsource
    )) |>
    dplyr::mutate(newbiologicalsource = biologicalsource)

  c <- paste0("\\b", latinGenitiveIsDic$genitive, "\\b")
  d <- latinGenitiveIsDic$nominative

  tcmNamesDic_is$newbiologicalsource <-
    stringi::stri_replace_all_regex(
      str = tcmNamesDic_is$newbiologicalsource,
      pattern = c,
      replacement = d,
      case_insensitive = FALSE,
      vectorize_all = FALSE
    )

  # example Adonis
  # example Agathis, Aleuritopteris, Amaryllis, ...

  # joining
  tcmNamesDicCurated <-
    dplyr::full_join(tcmNamesDic, tcmNamesDic_ae)

  # removing some organism parts
  e <- paste0("\\b", latinGenitivePartsDic$part, "\\b")
  f <- as.character(latinGenitivePartsDic$replacement)

  tcmNamesDicCurated$newbiologicalsource <-
    stringi::stri_replace_all_regex(
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
    rbind(tcmNamesDicCurated, tcmNamesDic_i, tcmNamesDic_is) |>
    dplyr::arrange(newbiologicalsource) |>
    dplyr::distinct(common, .keep_all = TRUE)

  # deleting ambiguous entries
  tcmNamesDicCurated <- tcmNamesDicCurated |>
    dplyr::filter(!common %in% manualSubtraction$common)

  # sorting in appropriate order
  tcmNamesDicCurated <- tcmNamesDicCurated |>
    dplyr::mutate(n = stringr::str_count(string = common)) |>
    dplyr::arrange(dplyr::desc(n)) |>
    # sorting for replacements like "sea cucumber" and so on...
    dplyr::select(
      vernacularName = common,
      canonicalName = biologicalsource,
      newCanonicalName = newbiologicalsource
    ) |>
    dplyr::distinct() |>
    dplyr::filter(vernacularName != canonicalName |
      is.na(newCanonicalName))

  tcmNamesDicCurated$vernacularName <-
    iconv(
      x = tcmNamesDicCurated$vernacularName,
      from = "UTF-8",
      to = "UTF-8",
      sub = ""
    )

  tcmNamesDicCurated$canonicalName <-
    iconv(
      x = tcmNamesDicCurated$canonicalName,
      from = "UTF-8",
      to = "UTF-8",
      sub = ""
    )

  tcmNamesDicCurated$newCanonicalName <-
    iconv(
      x = tcmNamesDicCurated$newCanonicalName,
      from = "UTF-8",
      to = "UTF-8",
      sub = ""
    )

  ## sampling rows for test mode
  "%ni%" <- Negate("%in%")
  set.seed(
    seed = 42,
    kind = "Mersenne-Twister",
    normal.kind = "Inversion"
  )
  tcmNamesDicCurated_sampled_1 <- tcmNamesDicCurated |>
    dplyr::sample_n(500)

  tcmNamesDicCurated_sampled_2 <- tcmNamesDicCurated |>
    dplyr::filter(grepl(
      pattern = paste(test_organisms$organismValue, collapse = "|"),
      x = vernacularName,
      ignore.case = TRUE
    ))

  tcmNamesDicCurated_sampled <-
    dplyr::bind_rows(
      tcmNamesDicCurated_sampled_1,
      tcmNamesDicCurated_sampled_2
    ) |>
    dplyr::mutate(n = stringr::str_count(string = vernacularName)) |>
    dplyr::arrange(dplyr::desc(n)) |>
    dplyr::select(-n)

  # exporting
  readr::write_delim(
    x = tcmNamesDicCurated,
    delim = "\t",
    file = gzfile(
      description = pathDataInterimDictionariesTcmNames,
      compression = 9,
      encoding = "UTF-8"
    ),
    quote = "none",
    escape = "double"
  )

  readr::write_delim(
    x = tcmNamesDicCurated_sampled,
    delim = "\t",
    file = gzfile(
      description = pathTestsDicTcmFile,
      compression = 9,
      encoding = "UTF-8"
    ),
    quote = "none",
    escape = "double"
  )
  ## because of univocity parser settings
}
