# title: "Organisms (sanitized) compileR"

# setting working directory
setwd("~/GitLab/opennaturalproductsdb/src/")

# loading paths
source("paths.R")

# loading functions
source("functions.R")

#writing path
## dictionaries
### taxa levels
taxaRanksDictionary <- read_delim(
  file = pathDataInterimDictionariesTaxaRanks,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

### black listed strings
blacklistDictionary <- read_delim(
  file = pathDataInterimDictionariesCommonBlackDic,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate(n = str_count(string = blackName)) %>%
  arrange(desc(n)) %>%
  select(-n)

## translated organisms
dataTranslatedOrganism <- read_delim(
  file = gzfile(pathTranslatedOrganism),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

dir <- pathTranslatedOrganismDistinct

length <- length(list.files(path = dir,
                            pattern = 'tsv'))

cut <- 10000

num <- as.integer(seq(
  from = 1 * cut,
  to = length * cut,
  by = cut
))

dataOrganismClean <- list()

# cleaning GNFinder output
for (i in num) {
  j <- i / cut
  tryCatch({
    dataOrganismClean[[j]] <-
      gnfinder_cleaning(num = i)
  }, error = function(e) {
    cat("ERROR :", conditionMessage(e), "\n")
  })
}

# filling missing organisms and taxonomies
dataOrganismFilled <- list()

for (i in 1:length(dataOrganismClean))
{
  tryCatch({
    dataOrganismFilled[[i]] <-
      biofilling(x = dataOrganismClean[[i]])
  }, error = function(e) {
    cat("ERROR :", conditionMessage(e), "\n")
  })
}

# explanation on taxize filling step: run:
# before <- dataOrganismClean[[1]] %>%
# filter(!is.na(canonicalname)) %>%
#   distinct(organismTranslated,
#            canonicalname,
#            .keep_all = TRUE)
# after <- dataOrganismFilled[[1]] %>%
# filter(!is.na(canonicalname)) %>%
#   distinct(organismTranslated,
#            canonicalname,
#            .keep_all = TRUE)
# diff <- anti_join(after, before)

# selecting and reordering
dataOrganismSanitized <- bind_rows(dataOrganismFilled) %>%
  select(
    organismTranslated,
    organismSanitized = canonicalname,
    organism_database = dbTaxo,
    everything()
  )

# adding original organism
dataOrganismSanitizedFilled <-
  left_join(dataTranslatedOrganism, dataOrganismSanitized)

# exporting
write.table(
  x = dataOrganismSanitizedFilled,
  file = gzfile(
    description = pathSanitizedOrganism,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
