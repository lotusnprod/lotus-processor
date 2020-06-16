# title: "MIBIG scrapeR"

# loading paths
source("paths.R")

# loading functions
source("functions.R")

df <- rbind(lapply(pathDataExternalDbSourceMibigOriginal, fromJSON))

x <- 1:length(df)

getid <- function(x) {
  j <- df[[x]][["cluster"]][["mibig_accession"]]
  j
}

id <- pbmclapply(
  FUN = getid,
  X = x,
  mc.preschedule = TRUE,
  mc.set.seed = TRUE,
  mc.silent = TRUE,
  mc.cores = (parallel::detectCores() - 2),
  mc.cleanup = TRUE,
  mc.allow.recursive = TRUE
)

getsmiles <- function(x) {
  j <- df[[x]][["cluster"]][["compounds"]][["chem_struct"]]
  j
}

smiles <- pbmclapply(
  FUN = getsmiles,
  X = x,
  mc.preschedule = TRUE,
  mc.set.seed = TRUE,
  mc.silent = TRUE,
  mc.cores = (parallel::detectCores() - 2),
  mc.cleanup = TRUE,
  mc.allow.recursive = TRUE
)

getorganism <- function(x) {
  j <- df[[x]][["cluster"]][["organism_name"]]
  j
}

organism <- pbmclapply(
  FUN = getorganism,
  X = x,
  mc.preschedule = TRUE,
  mc.set.seed = TRUE,
  mc.silent = TRUE,
  mc.cores = (parallel::detectCores() - 2),
  mc.cleanup = TRUE,
  mc.allow.recursive = TRUE
)

getreference <- function(x) {
  j <- df[[x]][["cluster"]][["publications"]]
  j
}

reference <- pbmclapply(
  FUN = getreference,
  X = x,
  mc.preschedule = TRUE,
  mc.set.seed = TRUE,
  mc.silent = TRUE,
  mc.cores = (parallel::detectCores() - 2),
  mc.cleanup = TRUE,
  mc.allow.recursive = TRUE
)

data <- tibble(id, smiles, organism, reference) %>%
  filter(smiles != "NULL") %>%
  unnest(id) %>%
  unnest(smiles) %>%
  unnest(organism) %>%
  unnest(reference) %>%
  distinct(id, .keep_all = TRUE) %>%
  filter(!is.na(smiles)) %>%
  select(id,
         smiles,
         biologicalsource = organism,
         reference)

data$name <- NA

# standardizing
data_standard <- standardizing_original(
  data_selected = data,
  db = "mib_1",
  structure_field = c("name", "smiles")
)

# exporting
write.table(
  x = data,
  file = gzfile(
    description = pathDataInterimDbMibig,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
