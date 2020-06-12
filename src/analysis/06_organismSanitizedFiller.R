# title: "Organisms (sanitized) filleR"

# setting working directory
setwd("~/GitLab/opennaturalproductsdb/src/")

# loading paths
source("paths.R")

# loading functions
source("functions.R")

# writing paths
## input
### files
filenames <- list.files(pathSanitizedOrganismDirTsv,
                        pattern = "*.tsv.zip",
                        full.names = TRUE)

# loading
## files
dataOrganismClean <- lapply(filenames,
                            function(x)
                            {
                              read_delim(
                                file = gzfile(x),
                                delim = "\t",
                                escape_double = FALSE,
                                trim_ws = TRUE
                              ) %>%
                                filter(!is.na(organismTranslated)) %>%
                                mutate_all(as.character)
                            })

## dictionaries
### taxa levels
taxaLevelsDictionary <- read_delim(
  file = pathInterimTaxaLevelsDic,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

### black listed strings
blackListDictionary <- read_delim(
  file = pathInterimBlackDic,
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

# selecting and reordering
dataOrganismSanitized <- bind_rows(dataOrganismFilled) %>%
  select(
    organismTranslated,
    organismSanitized = canonicalname,
    organism_database = db_taxo,
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
