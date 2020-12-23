# title: "TCM to scientific name translatoR"

# loading
## paths
source("paths.R")

library(Hmisc)
library(readxl)
library(splitstackshape)
library(stringi)
library(tidyverse)
library(vroom)

## functions
source("r/database.R")
source("r/y_as_na.R")
source("r/tcm_cleaning.R")
source("r/tcm_inverting.R")
source("r/vroom_safe.R")

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
tcmNamesDic_2 <- vroom(
  file = pathDataExternalTranslationSourceTcmTcmid,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE,
  col_names = TRUE,
  id = NULL,
  progress = TRUE,
  quote = ""
) %>%
  select(
    latin = `Latin Name`,
    common = `English Name`
  ) %>%
  filter(!is.na(common) &
    !is.na(latin)) %>%
  mutate(biologicalsource = latin)

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
latinGenitiveIDic <- read_delim(
  file = pathDataInterimDictionariesLatinGenitiveI,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate(n = str_count(string = genitive)) %>%
  arrange(desc(n)) %>%
  select(-n)

## is
latinGenitiveIsDic <- read_delim(
  file = pathDataInterimDictionariesLatinGenitiveIs,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate(n = str_count(string = genitive)) %>%
  arrange(desc(n)) %>%
  select(-n)

## parts
latinGenitivePartsDic <- read_delim(
  file = pathDataInterimDictionariesLatinPlantParts,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

## manually subtracted entries
manualSubtraction <- read_delim(
  file = pathDataInterimDictionariesTcmManualSubtraction,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
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

tcmNamesDic_2$latin <- tolower(x = tcmNamesDic_2$latin)
tcmNamesDic_2$latin <- capitalize(string = tcmNamesDic_2$latin)

tcmNamesDic_2$common <- tolower(x = tcmNamesDic_2$common)
tcmNamesDic_2$common <- capitalize(string = tcmNamesDic_2$common)

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

tcmNamesDic_2 <-
  tcm_cleaning(x = tcmNamesDic_2) %>%
  filter(!is.na(biologicalsource)) %>%
  distinct(biologicalsource, .keep_all = TRUE) %>%
  mutate(biologicalsource = latin) %>%
  select(common, biologicalsource)

# joining
tcmNamesDic <-
  rbind(
    tcmNamesDic_1_1,
    tcmNamesDic_1_2,
    tcmNamesDic_2,
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

a <- paste("\\b", latinGenitiveIDic$genitive, "\\b", sep = "")
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

c <- paste("\\b", latinGenitiveIsDic$genitive, "\\b", sep = "")
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
e <- paste("\\b", latinGenitivePartsDic$part, "\\b", sep = "")
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
  arrange(desc(n)) %>% # sorting for replacements like "sea cucumber" and so on...
  select(
    vernacularName = common,
    canonicalName = biologicalsource,
    newCanonicalName = newbiologicalsource
  ) %>%
  distinct() %>%
  filter(vernacularName != canonicalName | is.na(newCanonicalName))

# exporting
vroom_write_safe(
  x = tcmNamesDicCurated,
  path = pathDataInterimDictionariesTcmNames
)
