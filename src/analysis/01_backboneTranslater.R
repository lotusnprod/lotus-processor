# title: "Backbone translatoR"

# setting working directory
setwd("~/GitLab/opennaturalproductsdb/src/")

# loading paths
source("paths.R")

# loading functions
source("functions.R")

# loading files
## tcm names
tcmNamesDic <- read_delim(
  file = gzfile(pathInterimTcmNamesDic),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = FALSE
)

## common names
commonNamesDic <- read_delim(
  file = gzfile(pathInterimCommonNamesDic),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = FALSE
)

## taxa
taxa <- read_delim(
  pathTranslationSourceCommonGbifScientific,
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
                x = canonicalName)) #problems otherwise


## manual substraction
taxaRemovalDic <- read_delim(
  pathInterimTaxaManualSubtraction,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = FALSE
)

taxa <- taxa %>%
  mutate(specificEpithetCapitalized = capitalize(string = specificEpithet))

# deleting some (luckily) found erroneous entries
taxaSub <- taxa %>%
  filter(!genericName %in% taxaRemovalDic$genericName)

# quick replacement of the common names
## vs generic
problematicVernacularCommonGeneric <-
  left_join(taxaSub,
            commonNamesDic,
            by = c("genericName" = "vernacularName")) %>%
  filter(!is.na(canonicalName.y)) %>%
  distinct(canonicalName.x, .keep_all = TRUE) %>%
  select(name = canonicalName.x)

## vs specific
problematicVernacularCommonSpecific <-
  left_join(taxaSub,
            commonNamesDic,
            by = c("specificEpithetCapitalized" = "vernacularName")) %>%
  filter(!is.na(canonicalName.y)) %>%
  distinct(canonicalName.x, .keep_all = TRUE) %>%
  select(name = canonicalName.x)

# quick replacement of the common names
## vs generic
problematicVernacularTcmGeneric <-
  left_join(taxaSub,
            tcmNamesDic,
            by = c("genericName" = "vernacularName")) %>%
  filter(!is.na(canonicalName.y)) %>%
  distinct(canonicalName.x, .keep_all = TRUE) %>%
  select(name = canonicalName.x)

## vs specific
problematicVernacularTcmSpecific <-
  left_join(taxaSub,
            tcmNamesDic,
            by = c("specificEpithetCapitalized" = "vernacularName")) %>%
  filter(!is.na(canonicalName.y)) %>%
  distinct(canonicalName.x, .keep_all = TRUE) %>%
  select(name = canonicalName.x)

# joining
problematicVernacular <-
  rbind(
    problematicVernacularCommonGeneric,
    problematicVernacularCommonSpecific,
    problematicVernacularTcmGeneric,
    problematicVernacularTcmSpecific
  ) %>%
  distinct(name) %>%
  arrange(name)

problematicVernacularMin <- problematicVernacular %>%
  filter(str_count(string = name,
                   pattern =  "\\w+") == 1) %>%
  mutate(name = tolower(name))

problematicVernacularFull <-
  rbind(problematicVernacular,
        problematicVernacularMin) %>%
  mutate(nameMutated = paste("MUTATED", name, "mutated", sep = "")) %>%
  distinct(name, nameMutated) %>%
  arrange(name)

# exporting
write.table(
  x = problematicVernacularFull,
  file = gzfile(
    description = pathInterimProblematicNamesDic,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
