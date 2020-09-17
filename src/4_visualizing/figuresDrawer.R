# title: "Figures drawer"

### NOT FINISHED


cat("sourcing paths and functions \n")
source("functions.R")
source("paths.R")

cat("loading files, if running fullmode, this may take a while ... \n")

cat("... open DB triplets \n")
openDbTriplets <- read_delim(
  file = gzfile(pathDataInterimTablesAnalysedOpenDbTriplets),
  col_types = cols(.default = "c"),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  distinct(database,
           organismCleaned,
           structureCleanedInchikey3D) %>%
  data.frame()

cat("... GOLD open DB triplets \n")
openDbTripletsGold <- read_delim(
  file = gzfile(pathDataInterimTablesAnalysedGold),
  col_types = cols(.default = "c"),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)  %>%
  distinct(database,
           organismCleaned,
           structureCleanedInchikey3D) %>%
  data.frame()

cat("... DNP DB triplets \n")
dnpDbTriplets <- read_delim(
  file = gzfile(pathDataInterimTablesAnalysedDnpDbTriplets),
  col_types = cols(.default = "c"),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  data.frame() %>%
  distinct(database,
           organismCleaned,
           structureCleanedInchikey3D) %>%
  data.frame()

cat("... metadata ... \n")
cat("... organisms \n")
organismMetadata <- read_delim(
  file = gzfile(pathDataInterimDictionariesOrganismMetadata),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  data.frame() %>%
  mutate(
    organismCleaned_dbTaxo_1kingdom = ifelse(
      test = organismCleaned_dbTaxo_1kingdom == "Viridiplantae",
      yes = "Plantae",
      no = organismCleaned_dbTaxo_1kingdom
    ),
    organismCleaned_dbTaxo_1kingdom = ifelse(
      test = organismCleaned_dbTaxo_1kingdom == "Metazoa",
      yes = "Animalia",
      no = organismCleaned_dbTaxo_1kingdom
    ),
    organismCleaned_dbTaxo_1kingdom = ifelse(
      test = organismCleaned_dbTaxo_1kingdom == "Proteobacteria",
      yes = "Bacteria",
      no = organismCleaned_dbTaxo_1kingdom
    ),
    organismCleaned_dbTaxo_1kingdom = ifelse(
      test = organismCleaned_dbTaxo_1kingdom == "Protozoa",
      yes = "Protista",
      no = organismCleaned_dbTaxo_1kingdom
    ),
    organismCleaned_dbTaxo_2phylum = ifelse(
      test = organismCleaned_dbTaxo_2phylum == "Streptophyta" |
        organismCleaned_dbTaxo_2phylum == "Magnoliophyta" |
        organismCleaned_dbTaxo_2phylum == "Gymnospermophyta",
      yes = "Tracheophyta",
      no = organismCleaned_dbTaxo_2phylum
    ),
    organismCleaned_dbTaxo_1kingdom = ifelse(
      test = organismCleaned_dbTaxo_1kingdom == "Hepaticae",
      yes = "Plantae",
      no = organismCleaned_dbTaxo_1kingdom
    ),
    organismCleaned_dbTaxo_2phylum = ifelse(
      test = organismCleaned_dbTaxo_1kingdom == "Hepaticae",
      yes = "Marchantiophyta",
      no = organismCleaned_dbTaxo_2phylum
    ),
    organismCleaned_dbTaxo_1kingdom = ifelse(
      test = organismCleaned_dbTaxo_1kingdom == "Cyanobacteria",
      yes = "Bacteria",
      no = organismCleaned_dbTaxo_1kingdom
    ),
    organismCleaned_dbTaxo_2phylum = ifelse(
      test = organismCleaned_dbTaxo_1kingdom == "Cyanobacteria",
      yes = "Cyanobacteria",
      no = organismCleaned_dbTaxo_2phylum
    )
  )

cat("... structures (1) \n")
structureMetadata_1 <- read_delim(
  file = gzfile(pathDataInterimDictionariesStructureMetadata),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  data.frame()

cat("... structures (2) \n")
structureMetadata_2 <- read_delim(
  file = gzfile(
    "../../data/interim/09_structureContextualizer/classy.tsv.gz"
  ),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  data.frame()
