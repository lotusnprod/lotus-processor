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
  df1 <- dfsel %>%
    group_by(organismOriginal, organism_7_species) %>%
    fill(organism_8_variety, .direction = "downup") %>%
    group_by(organismOriginal, organism_6_genus) %>%
    fill(organism_7_species, .direction = "downup") %>%
    group_by(organismOriginal, organism_5_family) %>%
    fill(organism_6_genus, .direction = "downup") %>%
    group_by(organismOriginal, organism_4_order) %>%
    fill(organism_5_family, .direction = "downup") %>%
    group_by(organismOriginal, organism_3_class) %>%
    fill(organism_4_order, .direction = "downup") %>%
    group_by(organismOriginal, organism_2_phylum) %>%
    fill(organism_3_class, .direction = "downup") %>%
    group_by(organismOriginal, organism_1_kingdom) %>%
    fill(organism_2_phylum, .direction = "downup") %>%
    group_by(organismOriginal) %>%
    fill(organism_1_kingdom, .direction = "downup") %>%
    ungroup() %>%
    mutate(organismCleanedBis = as.character(apply(.[7:15], 1, function(x) {
      tail(na.omit(x), 1)
    })))

  df1$organismCleanedBis <-
    y_as_na(
      x = df1$organismCleanedBis,
      y = "character(0)"
    )

  df1$organismCleanedBis <-
    y_as_na(
      x = df1$organismCleanedBis,
      y = "NA"
    )

  df2 <- df1 %>%
    filter(organismCleaned == organismCleanedBis) %>%
    distinct(
      organismOriginal,
      organismDetected,
      organismDbTaxo,
      organismDbTaxoQuality,
      organismTaxonIds,
      organismTaxonRanks,
      organismTaxonomy,
      organismCleaned
    )

  df3 <- left_join(df2, dfsel)

  return(df3)
}

cat("Step 4 \n")
cat("loading cleaned organisms \n")
dataCleanedOrganismManipulated <- read_delim(
  file = gzfile(pathDataInterimTablesCleanedOrganismTranslatedTable),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = FALSE
) %>%
  select(everything(), organismDetected = organismCleaned) %>%
  relocate(organismDetected, .after = organismTaxonomy) %>%
  distinct()

dataCuratedOrganism <- dataCleanedOrganismManipulated %>%
  mutate(organismCleaned = as.character(apply(dataCleanedOrganismManipulated[7:15], 1, function(x) {
    tail(na.omit(x), 1)
  })))

dataCuratedOrganism$organismCleaned <-
  y_as_na(
    x = dataCuratedOrganism$organismCleaned,
    y = "character(0)"
  )

dataCuratedOrganism$organismCleaned <-
  y_as_na(
    x = dataCuratedOrganism$organismCleaned,
    y = "NA"
  )

cat("cleaning duplicate upstream taxa \n")
dataCuratedOrganismAuto <-
  taxo_cleaning_auto(dfsel = dataCuratedOrganism)

cat("selecting \n")
dataCuratedOrganismAuto[setdiff(
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
  y = names(dataCuratedOrganismAuto)
)] <- NA

dataCuratedOrganismAuto <- dataCuratedOrganismAuto %>%
  select(
    organismOriginal,
    organismDetected,
    organismCleaned,
    organismDbTaxo,
    organismDbTaxoQuality,
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
  x = dataCuratedOrganismAuto,
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