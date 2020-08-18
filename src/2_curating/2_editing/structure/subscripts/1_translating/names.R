# title: "Structure (nominal) translatoR"

# loading paths
source("paths.R")

# loading functions
source("functions/helpers.R")
source("functions/chemo.R")

# loading files
dataOriginal <- read_delim(
  file = gzfile(pathDataInterimTablesOriginalStructureNominal),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

# preparing names
dataTranslatedNominal <- preparing_name(x = dataOriginal)

# translating structures (cactus) (fast)
print(x = "translating structures with cactus (fast)")
dataTranslatedNominal_cactus <- dataTranslatedNominal %>%
  mutate(inchiNominal_cactus = invisible(
    pbmclapply(
      FUN = name2inchi_cactus,
      X = 1:nrow(dataTranslatedNominal),
      mc.preschedule = TRUE,
      mc.set.seed = TRUE,
      mc.silent = TRUE,
      mc.cores = (parallel::detectCores() - 2),
      mc.cleanup = TRUE,
      mc.allow.recursive = TRUE,
      ignore.interactive = TRUE
    )
  )) %>%
  mutate(inchiNominal_cactus = as.character(inchiNominal_cactus)) %>%
  mutate(inchiNominal_cactus = y_as_na(x = inchiNominal_cactus,
                                       y = "character(0)")) %>%
  mutate(inchiNominal_cactus = y_as_na(x = inchiNominal_cactus,
                                       y = "NA")) %>%
  mutate(inchiNominal_cactus = gsub(
    pattern = "^http.*",
    replacement = NA,
    x = inchiNominal_cactus
  )) %>%
  mutate(inchiNominal_cactus = gsub(
    pattern = "^NCI.*",
    replacement = NA,
    x = inchiNominal_cactus
  ))

dataTranslated <- left_join(dataOriginal,
                            dataTranslatedNominal_cactus) %>%
  mutate(structureTranslatedNominal_cactus = ifelse(
    test =  grepl(pattern = "^InChI=.*",
                  x = inchiNominal_cactus),
    yes = inchiNominal_cactus,
    no = NA
  )) %>%
  select(-inchiNominal_cactus)

# translating structures (cts) (slow but more results)
print(x = "translating structures with CTS (slow but more results)")
dataTranslatedNominal_cts <- dataTranslated %>%
  filter(is.na(structureTranslatedNominal_cactus))

dataTranslatedNominal_cts <- dataTranslatedNominal_cts %>%
  mutate(inchiNominal_cts = invisible(
    pbmclapply(
      FUN = name2inchi_cts,
      X = 1:nrow(dataTranslatedNominal_cts),
      mc.preschedule = TRUE,
      mc.set.seed = TRUE,
      mc.silent = TRUE,
      mc.cores = (parallel::detectCores() - 2),
      mc.cleanup = TRUE,
      mc.allow.recursive = TRUE,
      ignore.interactive = TRUE
    )
  )) %>%
  mutate(inchiNominal_cts = as.character(inchiNominal_cts)) %>%
  mutate(inchiNominal_cts = y_as_na(x = inchiNominal_cts,
                                    y = "character(0)")) %>%
  mutate(inchiNominal_cts = y_as_na(x = inchiNominal_cts,
                                    y = "NA")) %>%
  mutate(inchiNominal_cts = gsub(
    pattern = "^http.*",
    replacement = NA,
    x = inchiNominal_cts
  )) %>%
  mutate(inchiNominal_cts = gsub(
    pattern = "^NCI.*",
    replacement = NA,
    x = inchiNominal_cts
  ))

dataTranslated <- left_join(dataTranslated,
                            dataTranslatedNominal_cts) %>%
  mutate(structureTranslatedNominal_cts = ifelse(
    test =  grepl(pattern = "^InChI=.*",
                  x = inchiNominal_cts),
    yes = inchiNominal_cts,
    no = NA
  )) %>%
  select(-inchiNominal_cts) %>%
  mutate(
    structureTranslated_nominal = ifelse(
      test = !is.na(structureTranslatedNominal_cts),
      yes = structureTranslatedNominal_cts,
      no = structureTranslatedNominal_cactus
    )
  )

# exporting
## creating directories if they do not exist
ifelse(
  !dir.exists(pathDataInterimTablesTranslated),
  dir.create(pathDataInterimTablesTranslated),
  FALSE
)

ifelse(
  !dir.exists(pathDataInterimTablesTranslatedStructure),
  dir.create(pathDataInterimTablesTranslatedStructure),
  FALSE
)

write.table(
  x = dataTranslated,
  file = gzfile(
    description = pathDataInterimTablesTranslatedStructureNominal,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
