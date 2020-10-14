cat("This script performs chemical name translation. \n")

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")

cat("... functions \n")
source("functions/helpers.R")
source("functions/chemo.R")

cat("loading chemical names lists \n")
dataOriginal <- read_delim(
  file = gzfile(pathDataInterimTablesOriginalStructureNominal),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

cat("preparing names \n")
dataTranslatedNominal <- preparing_name(x = dataOriginal)

## abandoned because cactus returns too much very incorrect spotted...
### see https://cactus.nci.nih.gov/chemical/structure?identifier=terpinen-4-ol&representation=smiles
### or https://cactus.nci.nih.gov/chemical/structure?identifier=PONGAMOL&representation=smiles

# cat("translating structures with cactus (fast) \n")
# dataTranslatedNominal_cactus <- dataTranslatedNominal %>%
#   mutate(inchiNominal_cactus = invisible(
#     pbmclapply(
#       FUN = name2inchi_cactus,
#       X = 1:nrow(dataTranslatedNominal),
#       mc.preschedule = TRUE,
#       mc.set.seed = TRUE,
#       mc.silent = TRUE,
#       mc.cores = (parallel::detectCores() - 2),
#       mc.cleanup = TRUE,
#       mc.allow.recursive = TRUE,
#       ignore.interactive = TRUE
#     )
#   )) %>%
#   mutate(inchiNominal_cactus = as.character(inchiNominal_cactus)) %>%
#   mutate(inchiNominal_cactus = y_as_na(x = inchiNominal_cactus,
#                                        y = "character(0)")) %>%
#   mutate(inchiNominal_cactus = y_as_na(x = inchiNominal_cactus,
#                                        y = "NA")) %>%
#   mutate(inchiNominal_cactus = gsub(
#     pattern = "^http.*",
#     replacement = NA,
#     x = inchiNominal_cactus
#   )) %>%
#   mutate(inchiNominal_cactus = gsub(
#     pattern = "^NCI.*",
#     replacement = NA,
#     x = inchiNominal_cactus
#   ))
#
# dataTranslated <- left_join(dataOriginal,
#                             dataTranslatedNominal_cactus) %>%
#   mutate(structureTranslatedNominal_cactus = ifelse(
#     test =  grepl(pattern = "^InChI=.*",
#                   x = inchiNominal_cactus),
#     yes = inchiNominal_cactus,
#     no = NA
#   )) %>%
#   select(-inchiNominal_cactus)

# cat("translating structures with CTS (slow but more results) \n")
# dataTranslatedNominal_cts <- dataTranslated %>%
#   filter(is.na(structureTranslatedNominal_cactus))

cat("translating structures with CTS (slow but more results) \n")
dataTranslatedNominal_cts <- dataTranslatedNominal

if (nrow(dataTranslatedNominal_cts) == 0)
  dataTranslatedNominal_cts[1,] <- NA

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

# dataTranslated <- left_join(dataTranslated,
#                             dataTranslatedNominal_cts) %>%
#   mutate(structureTranslatedNominal_cts = ifelse(
#     test =  grepl(pattern = "^InChI=.*",
#                   x = inchiNominal_cts),
#     yes = inchiNominal_cts,
#     no = NA
#   )) %>%
#   select(-inchiNominal_cts) %>%
#   mutate(
#     structureTranslated_nominal = ifelse(
#       test = !is.na(structureTranslatedNominal_cts),
#       yes = structureTranslatedNominal_cts,
#       no = structureTranslatedNominal_cactus
#     )
#   )

dataTranslated <-  dataTranslatedNominal_cts %>%
  mutate(structureTranslatedNominal_cts = ifelse(
    test =  grepl(pattern = "^InChI=.*",
                  x = inchiNominal_cts),
    yes = inchiNominal_cts,
    no = NA
  )) %>%
  select(-inchiNominal_cts) %>%
  mutate(structureTranslated_nominal = structureTranslatedNominal_cts)

cat("ensuring directories exist \n")
ifelse(
  test = !dir.exists(pathDataInterimTablesTranslated),
  yes = dir.create(pathDataInterimTablesTranslated),
  no = paste(pathDataInterimTablesTranslated, "exists")
)

ifelse(
  test = !dir.exists(pathDataInterimTablesTranslatedStructure),
  yes = dir.create(pathDataInterimTablesTranslatedStructure),
  no = paste(pathDataInterimTablesTranslatedStructure, "exists")
)

cat("exporting ... \n")
cat(pathDataInterimTablesTranslatedStructureNominal, "\n")
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

end <- Sys.time()

cat("Script finished in", format(end - start), "\n")
