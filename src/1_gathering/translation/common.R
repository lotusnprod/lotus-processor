# title: "Common to scientific name translatoR"

# loading
## paths
source("paths.R")

library(dplyr)
library(readr)
library(splitstackshape)
library(stringr)

source("r/capitalize.R")
source("r/y_as_na.R")

if (mode == "full") {
  ##  files
  ### common names from PhenolExplorer
  commonSciPhe <-
    readr::read_delim(file = pathDataExternalTranslationSourceCommonPhenolexplorer) |>
    dplyr::select(
      vernacularName = name,
      canonicalName = food_source_scientific_name
    ) |>
    dplyr::filter(!is.na(vernacularName))

  ### common names from FooDB
  if (file.exists(pathDataExternalTranslationSourceCommonFoodb)) {
    commonSciFoo <-
      readr::read_delim(file = pathDataExternalTranslationSourceCommonFoodb) |>
      dplyr::select(
        vernacularName = name,
        canonicalName = name_scientific
      ) |>
      dplyr::filter(!is.na(vernacularName))
  }

  ### common names from DrDuke
  commonDuk <-
    readr::read_delim(file = pathDataExternalTranslationSourceCommonDrdukeCommon) |>
    dplyr::select(
      vernacularName = CNNAM,
      FNFNUM
    )

  ### scientific names from DrDuke
  sciDuk <-
    readr::read_delim(file = pathDataExternalTranslationSourceCommonDrdukeScientific) |>
    dplyr::select(FNFNUM,
      canonicalName = TAXON
    )

  ## forcing organisms in test file to be in sampled dic
  test_organisms <- readr::read_delim(
    file = "../tests/tests.tsv",
    delim = "\t"
  )

  commonSciDuk <- dplyr::left_join(sciDuk, commonDuk) |>
    dplyr::select(-FNFNUM) |>
    dplyr::filter(!is.na(vernacularName))

  ### GBIF
  #### taxa
  taxa <-
    readr::read_delim(
      file = unz(
        description = pathDataExternalTranslationSourceCommonGbif,
        filename = "backbone/Taxon.tsv"
      )
    ) |>
    dplyr::filter(!is.na(canonicalName)) |>
    dplyr::distinct(
      taxonID,
      canonicalName,
      genericName,
      specificEpithet
    ) |>
    dplyr::filter(!grepl(
      pattern = "\\?",
      x = canonicalName
    ))

  #### taxa
  vernacular <- readr::read_delim(
    file = unz(
      description = pathDataExternalTranslationSourceCommonGbif,
      filename = "backbone/VernacularName.tsv"
    )
  ) |>
    dplyr::filter(language == "en") |>
    dplyr::distinct(
      taxonID,
      vernacularName
    ) |>
    dplyr::filter(!grepl(
      pattern = "\\?",
      x = vernacularName
    ))

  ### manually subtracted entries
  manualSubtraction <-
    readr::read_delim(
      file = pathDataInterimDictionariesCommonManualSubtraction,
      delim = "\t"
    )

  # removing parts of PhenolExplorer names
  commonSciPhe$vernacularName <-
    sub(
      pattern = "(\\[.*)",
      replacement = "",
      x = commonSciPhe$vernacularName
    )
  commonSciPhe$vernacularName <-
    sub(
      pattern = "(\\,.*)",
      replacement = "",
      x = commonSciPhe$vernacularName
    )

  # joining taxa and vernacular names from GBIF
  taxaVernacular <- dplyr::left_join(taxa, vernacular) |>
    dplyr::filter(!is.na(vernacularName)) |>
    dplyr::distinct(canonicalName, vernacularName) |>
    dplyr::group_by(vernacularName) |>
    dplyr::arrange(dplyr::desc(stringr::str_count(canonicalName))) |>
    dplyr::ungroup() |>
    dplyr::distinct(vernacularName, .keep_all = TRUE) |>
    dplyr::arrange(dplyr::desc(stringr::str_count(vernacularName))) |>
    dplyr::filter(canonicalName != "Boa constrictor")

  # deleting vernacular names corresponding to generic epithets for safety reasons
  ## they are almost safe (see Cacao) but just to be on the safe side...
  # list <- commonSciSub |>
  #   dplyr::filter(vernacularName %in% taxa$genericName)
  taxaVernacular <- taxaVernacular |>
    dplyr::filter(!tolower(vernacularName) %in% tolower(taxa$genericName))

  # joining common names from PhenolExplorer and FooDB
  if (file.exists(pathDataExternalTranslationSourceCommonFoodb)) {
    commonSciPheFoo <- dplyr::full_join(commonSciPhe, commonSciFoo)
  } else {
    commonSciPheFoo <- commonSciPhe
  }

  # joining common names from DrDukes
  commonSciPheFooDuk <-
    dplyr::full_join(commonSciPheFoo, commonSciDuk) |>
    dplyr::filter(!tolower(vernacularName) %in% tolower(taxa$genericName))

  # adding
  ## plurals
  # library(textclean) DOES NOT GIVE GOOD RESULTS (ex. cherry tomatoes)
  # library(SemNetCleaner) also does not work (ex. Zebrafish)
  # food2sci$plural <- make_plural(food2sci$name)

  # deleting vernacular names corresponding to generic epithets for safety reasons
  ## they are almost safe (see Cacao) but just to be on the safe side...
  # list <- commonSciSub %>%
  #   filter(vernacularName %in% taxa$genericName)

  commonSciPheFooDuk <- commonSciPheFooDuk

  ### normal
  commonSciPlural_1 <- commonSciPheFooDuk |>
    dplyr::filter(
      !grepl(
        pattern = "(.+[^aeiou])y$",
        x = vernacularName
      ) &
        !grepl(
          pattern = "(.+[^aeiou])o$",
          x = vernacularName
        ) &
        !grepl(
          pattern = ".+s$|.+sh$|.+ch$.+x$|.+z$",
          x = vernacularName
        ) &
        !grepl(
          pattern = ".+f$|.+fe$",
          x = vernacularName
        )
    )

  commonSciPlural_1$vernacularName <-
    paste0(commonSciPlural_1$vernacularName, "s")

  ### "o" and "y" plurals (mango -> mangoes, berry -> berries)
  commonSciPlural_2 <- commonSciPheFooDuk |>
    dplyr::filter(
      grepl(
        pattern = "(.+[^aeiou])y$",
        x = vernacularName
      ) |
        grepl(
          pattern = "(.+[^aeiou])o$",
          x = vernacularName
        )
    )

  commonSciPlural_2$vernacularName <- gsub(
    pattern = "(.+[^aeiou])y$",
    replacement = "\\1ies",
    x = commonSciPlural_2$vernacularName
  )

  commonSciPlural_2$vernacularName <- gsub(
    pattern = "(.+[^aeiuy])o$",
    replacement = "\\1oes",
    x = commonSciPlural_2$vernacularName
  )

  ### s, sh, ch, x, z
  commonSciPlural_3 <- commonSciPheFooDuk |>
    dplyr::filter(grepl(
      pattern = ".+s$|.+sh$|.+ch$.+x$|.+z$",
      x = vernacularName
    ))

  commonSciPlural_3$vernacularName <-
    paste0(commonSciPlural_3$vernacularName, "es")

  ### f, fe
  commonSciPlural_4 <- commonSciPheFooDuk |>
    dplyr::filter(grepl(
      pattern = "(.+)f$|(.+)fe$",
      x = vernacularName
    ))

  commonSciPlural_4$vernacularName <- gsub(
    pattern = "(.+)f$|(.+)fe$",
    replacement = "\\1ves",
    x = commonSciPlural_4$vernacularName
  )

  # joining
  commonSciPluralized <-
    rbind(
      commonSciPheFooDuk,
      commonSciPlural_1,
      commonSciPlural_2,
      commonSciPlural_3,
      commonSciPlural_4
    ) |>
    dplyr::distinct(vernacularName, canonicalName)

  # joining common names from GBIF
  commonSci <- dplyr::full_join(commonSciPluralized, taxaVernacular)

  # capitalizing
  commonSci$vernacularName <-
    capitalize(string = commonSci$vernacularName)

  # removing some useless characters
  commonSci$vernacularName <- trimws(x = commonSci$vernacularName)

  commonSci <- commonSci %>%
    dplyr::mutate_all(~ iconv(x = ., from = "utf-8", to = "utf-8//ignore"))

  commonSci$vernacularName <- gsub(
    pattern = "/",
    replacement = " ",
    x = commonSci$vernacularName
  )

  commonSci$vernacularName <- gsub(
    pattern = "\\s*\\([^\\)]+\\)",
    replacement = "",
    x = commonSci$vernacularName
  )

  commonSci$vernacularName <- gsub(
    pattern = "\\s*\\([^\\)]",
    replacement = "",
    x = commonSci$vernacularName
  )

  # removing approximative additions of specific names
  commonSci <- commonSci |>
    dplyr::filter(vernacularName != stringr::word(
      string = canonicalName,
      start = 1
    ))

  ## explanation
  explanation <- commonSci |>
    dplyr::filter(vernacularName == stringr::word(
      string = canonicalName,
      start = 1
    ))

  # filtering only results with canonical name
  commonSci$canonicalName <- y_as_na(
    x = commonSci$canonicalName,
    y = "N/A"
  )

  commonSci$canonicalName <- y_as_na(
    x = commonSci$canonicalName,
    y = ""
  )

  commonSci$canonicalName <- trimws(commonSci$canonicalName)

  commonSci$vernacularName <- trimws(commonSci$vernacularName)

  commonSci <- commonSci |>
    dplyr::filter(!is.na(canonicalName)) |>
    dplyr::distinct(vernacularName, canonicalName)

  # filtering common names with only one translation
  commonSci_1 <- commonSci |>
    dplyr::arrange(canonicalName) |>
    dplyr::group_by(vernacularName) |>
    dplyr::add_count() |>
    dplyr::filter(n == 1) |>
    dplyr::select(-n) |>
    dplyr::ungroup() |>
    dplyr::arrange(vernacularName)

  # filtering common names with more than one translation
  commonSci_2 <- commonSci |>
    dplyr::arrange(canonicalName) |>
    dplyr::group_by(vernacularName) |>
    dplyr::add_count(name = "vernacularCount") |>
    dplyr::filter(vernacularCount != 1) |>
    dplyr::arrange(vernacularName) |>
    splitstackshape::cSplit(
      splitCols = "canonicalName",
      sep = " ",
      drop = FALSE
    ) |>
    dplyr::group_by(vernacularName, canonicalName_1, canonicalName_2) |>
    dplyr::add_count(name = "specificCount") |>
    dplyr::group_by(vernacularName, canonicalName_1) |>
    dplyr::add_count(name = "genericCount") |>
    dplyr::ungroup() |>
    dplyr::mutate(
      specificRatio = specificCount / vernacularCount,
      genericRatio = genericCount / vernacularCount
    )

  # filtering ambiguous entries and outliers
  commonSci_3 <- commonSci_2 |>
    dplyr::filter(specificRatio > 0.5 | genericRatio > 0.5) |>
    dplyr::group_by(vernacularName) |>
    dplyr::add_count(name = "vernacularCount") |>
    dplyr::group_by(vernacularName, canonicalName_1, canonicalName_2) |>
    dplyr::add_count(name = "specificCount") |>
    dplyr::group_by(vernacularName, canonicalName_1) |>
    dplyr::add_count(name = "genericCount") |>
    dplyr::ungroup()

  commonSciAmbiguous <- commonSci_2 |>
    dplyr::filter(specificRatio == 0.5 & genericRatio == 0.5) |>
    dplyr::group_by(vernacularName) |>
    dplyr::add_count(name = "vernacularCount") |>
    dplyr::group_by(vernacularName, canonicalName_1, canonicalName_2) |>
    dplyr::add_count(name = "canonicalCount") |>
    dplyr::ungroup()

  ## specific names matching
  commonSci_4 <- commonSci_3 |>
    dplyr::filter(specificRatio > 0.5) |>
    dplyr::arrange(dplyr::desc(specificRatio)) |>
    dplyr::distinct(vernacularName, .keep_all = TRUE) |>
    dplyr::mutate(newCanonicalName = paste(canonicalName_1,
      canonicalName_2,
      sep = " "
    )) |>
    dplyr::select(vernacularName,
      canonicalName = newCanonicalName
    ) |>
    dplyr::arrange(vernacularName)

  ## generic name matching
  commonSci_5 <- dplyr::anti_join(commonSci_3,
    commonSci_4,
    by = "vernacularName"
  ) |>
    dplyr::filter(specificRatio <= 0.5 &
      genericRatio > 0.5) |>
    dplyr::arrange(dplyr::desc(genericRatio)) |>
    dplyr::distinct(vernacularName, .keep_all = TRUE) |>
    dplyr::select(vernacularName,
      canonicalName = canonicalName_1
    ) |>
    dplyr::arrange(vernacularName)

  # joining again cleaned results
  commonSciJoined <- rbind(commonSci_1, commonSci_4, commonSci_5)

  # deleting ambiguous entries
  commonSciSub <- commonSciJoined |>
    dplyr::filter(!tolower(vernacularName) %in% tolower(manualSubtraction$name))

  commonSciSub$canonicalName <-
    y_as_na(commonSciSub$canonicalName, "\"\"")

  ## sorting in appropriate order
  common2Sci <- commonSciSub |>
    dplyr::mutate(n = stringr::str_count(string = vernacularName)) |>
    dplyr::arrange(dplyr::desc(n)) %>%
    # sorting for replacements like "sea cucumber" and so on...
    dplyr::filter(n >= 4) %>%
    # names with 3 char are not enough
    dplyr::select(
      vernacularName,
      canonicalName
    ) |>
    dplyr::filter(!grepl(pattern = "\\?", x = canonicalName)) |>
    dplyr::filter(!grepl(pattern = "\\)", x = vernacularName)) |>
    splitstackshape::cSplit(
      "vernacularName",
      sep = "   ",
      fixed = TRUE,
      direction = "long"
    )

  common2Sci$canonicalName <-
    iconv(
      x = common2Sci$canonicalName,
      from = "UTF-8",
      to = "UTF-8",
      sub = ""
    )

  common2Sci$vernacularName <-
    iconv(
      x = common2Sci$vernacularName,
      from = "UTF-8",
      to = "UTF-8",
      sub = ""
    )

  common2Sci <- common2Sci |>
    dplyr::mutate(vernacularName = gsub(
      pattern = "\\[.*\\]",
      replacement = "",
      x = vernacularName
    )) |>
    dplyr::mutate(vernacularName = gsub(
      pattern = "\\[",
      replacement = "",
      x = vernacularName
    )) |>
    dplyr::mutate(vernacularName = gsub(
      pattern = "\\]",
      replacement = "",
      x = vernacularName
    ))

  ## New in GBIF backbone
  common2Sci <- common2Sci |>
    dplyr::filter(!grepl(
      pattern = "| ",
      x = vernacularName,
      fixed = TRUE
    ))

  ## sampling rows for test mode
  "%ni%" <- Negate("%in%")
  set.seed(
    seed = 42,
    kind = "Mersenne-Twister",
    normal.kind = "Inversion"
  )
  common2Sci_sampled_1 <- common2Sci |>
    dplyr::sample_n(500)

  common2Sci_sampled_2 <- common2Sci |>
    dplyr::filter(grepl(
      pattern = paste(test_organisms$organismValue, collapse = "|"),
      x = vernacularName,
      ignore.case = TRUE
    ))

  common2Sci_sampled <-
    dplyr::bind_rows(
      common2Sci_sampled_1,
      common2Sci_sampled_2
    ) |>
    dplyr::mutate(n = stringr::str_count(string = vernacularName)) |>
    dplyr::arrange(dplyr::desc(n)) |>
    dplyr::select(-n)

  # exporting
  readr::write_delim(
    x = common2Sci,
    delim = "\t",
    file = gzfile(
      description = pathDataInterimDictionariesCommonNames,
      compression = 9,
      encoding = "UTF-8"
    ),
    quote = "none",
    escape = "double"
  )

  readr::write_delim(
    x = common2Sci_sampled,
    delim = "\t",
    file = gzfile(
      description = pathTestsDicCommonFile,
      compression = 9,
      encoding = "UTF-8"
    ),
    quote = "none",
    escape = "double"
  )
  ## because of univocity parser settings
}
