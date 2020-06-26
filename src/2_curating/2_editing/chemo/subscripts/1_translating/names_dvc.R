# title: "Structure (nominal) translatoR"

# loading paths
source("paths.R")

# loading functions
source("functions/chemo.R")

# loading files
dataOriginal <- read_delim(
  file = gzfile(pathOriginalStructureNominal),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

# preparing names
dataTranslatedNominal <- preparing_name(x = dataOriginal)

# preparing structures
dataTranslatedNominal$inchiNominal <- invisible(
  pbmclapply(
    FUN = name2inchi,
    X = 1:nrow(dataTranslatedNominal),
    mc.preschedule = TRUE,
    mc.set.seed = TRUE,
    mc.silent = TRUE,
    mc.cores = (parallel::detectCores() - 2),
    mc.cleanup = TRUE,
    mc.allow.recursive = TRUE
  )
)

dataTranslatedNominal$inchiNominal <-
  as.character(dataTranslatedNominal$inchiNominal)
dataTranslatedNominal$inchiNominal <-
  y_as_na(x = dataTranslatedNominal$inchiNominal,
          y = "character(0)")
dataTranslatedNominal$inchiNominal <-
  y_as_na(x = dataTranslatedNominal$inchiNominal,
          y = "NA")
dataTranslatedNominal$inchiNominal <- gsub(pattern = "^http.*",
                                           replacement = NA,
                                           x = dataTranslatedNominal$inchiNominal)
dataTranslatedNominal$inchiNominal <- gsub(pattern = "^NCI.*",
                                           replacement = NA,
                                           x = dataTranslatedNominal$inchiNominal)

dataTranslated <- left_join(dataOriginal, dataTranslatedNominal) %>%
  mutate(structureTranslatedNominal = ifelse(
    test =  grepl(pattern = "^InChI=.*",
                  x = inchiNominal),
    yes = inchiNominal,
    no = NA
  )) %>%
  select(-inchiNominal)

# exporting
write.table(
  x = dataTranslated,
  file = gzfile(
    description = pathTranslatedStructureNominal_dvc,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
