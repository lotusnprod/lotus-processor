cat("This script performs taxonomy alignment. \n")

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")

cat("... functions \n")
source("functions/log.R")
source("functions/helpers.R")
source("functions/bio.R")
source("2_curating/2_editing/organism/functions/manipulating_taxo.R")

log_debug(" Step 4")
cat("defining function \n")
taxo_cleaning_auto <- function(dfsel) {
  dfsel$count <-
    rowSums(!is.na(dfsel[(grepl(
      pattern = "organism_",
      x = colnames(dfsel),
      fixed = TRUE
    ))]))
  
  df1 <- dfsel %>%
    group_by(organismOriginal) %>%
    mutate(count_max = max(count)) %>%
    ungroup() %>% 
    filter(count == count_max) %>%
    select(-count, -count_max)
  
  return(df1)
}

cat("Step 4 \n")
cat("loading cleaned organisms \n")
dataCleanedOrganismManipulated <- read_delim(
  file = gzfile(pathDataInterimTablesCleanedOrganismTranslatedTable),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = FALSE
) %>%
  relocate(organismCleaned, .after = organismTaxonomy)

# curating taxonomy
## auto
dataCuratedOrganismAuto <-
  taxo_cleaning_auto(dfsel = dataCleanedOrganismManipulated)

## manual
# dataCuratedOrganism <-
#   taxo_cleaning_manual(dfsel = dataCleanedOrganismManipulated)

cat("keeping lowest taxon \n")
dataCuratedOrganism <- dataCuratedOrganismAuto %>%
  mutate(organismCleaned = as.character(apply(dataCuratedOrganismAuto[7:15], 1, function(x) {
    tail(na.omit(x), 1)
  })))

dataCuratedOrganism$organismCleaned <-
  y_as_na(x = dataCuratedOrganism$organismCleaned,
          y = "character(0)")

dataCuratedOrganism$organismCleaned <-
  y_as_na(x = dataCuratedOrganism$organismCleaned,
          y = "NA")

cat("selecting \n")
dataCuratedOrganism[setdiff(
  x = c(
    "organismOriginal",
    "organismCleaned",
    "organismDbTaxo",
    "organismDbTaxoQuality",
    "organismTaxonIds",
    "organismTaxonRanks",
    "organismTaxonomy",
    "organism_1_kingdom",
    "organism_2_phylum",
    "organism_3_class",
    "organism_4_order",
    "organism_5_family",
    "organism_6_genus",
    "organism_7_species",
    "organism_8_quality"
  ),
  y = names(dataCuratedOrganism)
)] <- NA

dataCuratedOrganism <- dataCuratedOrganism %>%
  select(
    organismOriginal,
    organismCleaned,
    organismDbTaxo,
    organismDbTaxoQuality,
    # organismModifiedTaxonomyAuto = organism_modified_taxonomy_auto,
    # organismModifiedTaxonomyManual = organism_modified_taxonomy_manual,
    organismTaxonIds,
    organismTaxonRanks,
    organismTaxonomy,
    organism_1_kingdom,
    organism_2_phylum,
    organism_3_class,
    organism_4_order,
    organism_5_family,
    organism_6_genus,
    organism_7_species,
    organism_8_variety
  )

cat("exporting ... \n")
cat(pathDataInterimTablesCleanedOrganismFinal, "\n")
write.table(
  x = dataCuratedOrganism,
  file = gzfile(
    description = pathDataInterimTablesCleanedOrganismFinal,
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