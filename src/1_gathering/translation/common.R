# title: "Common to scientific name translatoR"

# loading paths
source("paths.R")

# loading functions
source("functions.R")

# loading files
## common names from PhenolExplorer
commonSciPhe <- read_delim(
  file = pathDataExternalTranslationSourceCommonPhenolexplorer,
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  select(vernacularName = name,
         canonicalName = food_source_scientific_name) %>%
  filter(!is.na(vernacularName))

## common names from FooDB
commonSciFoo <- read_delim(
  file = pathDataExternalTranslationSourceCommonFoodb,
  delim = ",",
  escape_double = TRUE,
  trim_ws = TRUE
) %>%
  select(vernacularName = name,
         canonicalName = name_scientific) %>%
  filter(!is.na(vernacularName))

## common names from DrDuke
commonDuk <- read_delim(
  file = pathDataExternalTranslationSourceCommonDrdukeCommon,
  delim = ",",
  escape_double = TRUE,
  trim_ws = TRUE
) %>%
  select(vernacularName = CNNAM,
         FNFNUM)

sciDuk <- read_delim(
  file = pathDataExternalTranslationSourceCommonDrdukeScientific,
  delim = ",",
  escape_double = TRUE,
  trim_ws = TRUE
) %>%
  select(FNFNUM,
         canonicalName = TAXON)

commonSciDuk <- left_join(sciDuk, commonDuk) %>%
  select(-FNFNUM) %>%
  filter(!is.na(vernacularName))

# GBIF
taxa <- read_delim(
  file = pathDataExternalTranslationSourceCommonGbifScientific,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = FALSE,
  progress = TRUE
) %>%
  filter(!is.na(canonicalName)) %>%
  select(taxonID,
         canonicalName,
         genericName,
         specificEpithet) %>%
  filter(!grepl(pattern = "\\?",
                x = canonicalName))

vernacular <- read_delim(
  file = pathDataExternalTranslationSourceCommonGbifVernacular,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = FALSE
) %>%
  filter(language == "en") %>%
  select(taxonID,
         vernacularName) %>%
  filter(!grepl(pattern = "\\?",
                x = vernacularName))

# manually subtracted entries
manualSubtraction <- read_delim(
  file = pathDataInterimDictionariesCommonManualSubtraction,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
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
taxaVernacular <- left_join(taxa, vernacular) %>%
  filter(!is.na(vernacularName)) %>%
  distinct(canonicalName, vernacularName) %>%
  group_by(vernacularName) %>%
  arrange(desc(str_count(canonicalName))) %>%
  ungroup() %>%
  distinct(vernacularName, .keep_all = TRUE) %>%
  arrange(desc(str_count(vernacularName)))

# joining common names from PhenolExplorer and FooDB
commonSciPheFoo <- full_join(commonSciPhe, commonSciFoo)

# joining common names from DrDukes
commonSciPheFooDuk <- full_join(commonSciPheFoo, commonSciDuk)

# adding
## plurals
# library(textclean) DOES NOT GIVE GOOD RESULTS (ex. cherry tomatoes)
# library(SemNetCleaner) also does not work (ex. Zebrafish)
# food2sci$plural <- make_plural(food2sci$name)

### normal
commonSciPlural_1 <- commonSciPheFooDuk %>%
  filter(
    !grepl(pattern = "(.+[^aeiou])y$",
           x = vernacularName) &
      !grepl(pattern = "(.+[^aeiou])o$",
             x = vernacularName) &
      !grepl(pattern = ".+s$|.+sh$|.+ch$.+x$|.+z$",
             x =  vernacularName) &
      !grepl(pattern = ".+f$|.+fe$",
             x = vernacularName)
  )

commonSciPlural_1$vernacularName <-
  paste(commonSciPlural_1$vernacularName, "s", sep = "")

### "o" and "y" plurals (mango -> mangoes, berry -> berries)
commonSciPlural_2 <- commonSciPheFooDuk %>%
  filter(
    grepl(pattern = "(.+[^aeiou])y$",
          x = vernacularName) |
      grepl(pattern = "(.+[^aeiou])o$",
            x = vernacularName)
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
commonSciPlural_3 <- commonSciPheFooDuk %>%
  filter(grepl(pattern = ".+s$|.+sh$|.+ch$.+x$|.+z$",
               x = vernacularName))

commonSciPlural_3$vernacularName <-
  paste(commonSciPlural_3$vernacularName, "es", sep = "")

### f, fe
commonSciPlural_4 <- commonSciPheFooDuk %>%
  filter(grepl(pattern = "(.+)f$|(.+)fe$",
               x = vernacularName))

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
  ) %>%
  distinct(vernacularName, canonicalName)

# joining common names from GBIF
commonSci <- full_join(commonSciPluralized, taxaVernacular)

# capitalizing
commonSci$vernacularName <-
  capitalize(string = commonSci$vernacularName)

# removing some useless characters
commonSci$vernacularName <- trimws(x = commonSci$vernacularName)

commonSci <- commonSci %>%
  mutate_if(is.character, ~ gsub(
    pattern = '[^ -~]',
    replacement = '',
    x = .
  ))

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
  replacement =  "",
  x =  commonSci$vernacularName
)

# removing approximative additions of specific names
commonSci <- commonSci %>%
  filter(vernacularName != word(string = canonicalName,
                                start =  1))

# filtering only results with canonical name
commonSci$canonicalName <- y_as_na(x = commonSci$canonicalName,
                                   y =  "N/A")

commonSci$canonicalName <- y_as_na(x = commonSci$canonicalName,
                                   y = "")

commonSci$canonicalName <- trimws(commonSci$canonicalName)

commonSci$vernacularName <- trimws(commonSci$vernacularName)

commonSci <- commonSci %>%
  filter(!is.na(canonicalName)) %>%
  distinct(vernacularName, canonicalName)

# filtering common names with only one translation
commonSci_1 <- commonSci %>%
  arrange(canonicalName) %>%
  group_by(vernacularName) %>%
  add_count() %>%
  filter(n == 1) %>%
  select(-n) %>%
  ungroup() %>%
  arrange(vernacularName)

# filtering common names with more than one translation
commonSci_2 <- commonSci %>%
  arrange(canonicalName) %>%
  group_by(vernacularName) %>%
  add_count(name = "vernacularCount") %>%
  filter(vernacularCount != 1) %>%
  arrange(vernacularName) %>%
  cSplit(splitCols = "canonicalName",
         sep = " ",
         drop = FALSE) %>%
  group_by(vernacularName, canonicalName_1, canonicalName_2) %>%
  add_count(name = "specificCount") %>%
  group_by(vernacularName, canonicalName_1) %>%
  add_count(name = "genericCount") %>%
  ungroup() %>%
  mutate(
    specificRatio = specificCount / vernacularCount,
    genericRatio = genericCount / vernacularCount
  )

# filtering ambiguous entries and outliers
commonSci_3 <- commonSci_2 %>%
  filter(specificRatio > 0.5 | genericRatio > 0.5) %>%
  group_by(vernacularName) %>%
  add_count(name = "vernacularCount") %>%
  group_by(vernacularName, canonicalName_1, canonicalName_2) %>%
  add_count(name = "specificCount") %>%
  group_by(vernacularName, canonicalName_1) %>%
  add_count(name = "genericCount") %>%
  ungroup()

commonSciAmbiguous <- commonSci_2 %>%
  filter(specificRatio == 0.5 & genericRatio == 0.5) %>%
  group_by(vernacularName) %>%
  add_count(name = "vernacularCount") %>%
  group_by(vernacularName, canonicalName_1, canonicalName_2) %>%
  add_count(name = "canonicalCount") %>%
  ungroup()

## specific names matching
commonSci_4 <- commonSci_3 %>%
  filter(specificRatio > 0.5) %>%
  arrange(desc(specificRatio)) %>%
  distinct(vernacularName, .keep_all = TRUE) %>%
  mutate(newCanonicalName = paste(canonicalName_1,
                                  canonicalName_2,
                                  sep = " ")) %>%
  select(vernacularName,
         canonicalName = newCanonicalName) %>%
  arrange(vernacularName)

## generic name matching
commonSci_5 <- anti_join(commonSci_3,
                         commonSci_4,
                         by = "vernacularName") %>%
  filter(specificRatio <= 0.5 &
           genericRatio > 0.5) %>%
  arrange(desc(genericRatio)) %>%
  distinct(vernacularName, .keep_all = TRUE) %>%
  select(vernacularName,
         canonicalName = canonicalName_1) %>%
  arrange(vernacularName)

# joining again cleaned results
commonSciJoined <- rbind(commonSci_1, commonSci_4, commonSci_5)

# deleting ambiguous entries
commonSciSub <- commonSciJoined %>%
  filter(!vernacularName %in% manualSubtraction$name)

## sorting in appropriate order
common2Sci <- commonSciSub %>%
  mutate(n = str_count(string = vernacularName)) %>%
  arrange(desc(n)) %>%  #sorting for replacements like "sea cucumber" and so on...
  filter(n >= 4) %>% #names with 3 char are not enough
  select(vernacularName,
         canonicalName) %>%
  filter(!grepl("\\?", canonicalName)) %>%
  filter(!grepl("\\)", vernacularName))

# exporting
write.table(
  x = common2Sci,
  file = gzfile(
    description = pathDataInterimDictionariesCommonNames,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)