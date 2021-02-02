cat("This script outputs some metrics related to the DB \n")

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")

cat("... libraries \n")
library(tidyverse)
library(data.table)
source("r/vroom_safe.R")

cat("loading ... \n")
cat("databases list ... \n")
dataset <- read_delim(
  file = "../docs/dataset.tsv",
  delim = "\t"
) %>%
  select(
    -`initial retrieved unique observations`,
    -`cleaned documented structure-organism pairs`,
    -`pairs validated for wikidata export`
  )

cat("initial table ... \n")
dbTable <- lapply(pathDataInterimDbDir, vroom_read_safe) %>%
  rbindlist(l = ., fill = TRUE) %>%
  select(
    database,
    organismOriginal_clean = organism_clean,
    organismOriginal_dirty = organism_dirty,
    structureOriginal_inchi = structure_inchi,
    structureOriginal_nominal = structure_name,
    structureOriginal_smiles = structure_smiles,
    referenceOriginal_authors = reference_authors,
    referenceOriginal_doi = reference_doi,
    referenceOriginal_external = reference_external,
    referenceOriginal_isbn = reference_isbn,
    referenceOriginal_journal = reference_journal,
    referenceOriginal_original = reference_original,
    referenceOriginal_pubmed = reference_pubmed,
    referenceOriginal_publishingDetails = reference_publishingDetails,
    referenceOriginal_split = reference_split,
    referenceOriginal_title = reference_title,
  ) %>%
  tibble()

cat("final table ... \n")
inhouseDbMinimal <-
  vroom_read_safe(path = pathDataInterimTablesCuratedTable)

cat("validated for export ... \n")
openDb <-
  vroom_read_safe(path = pathDataInterimTablesAnalysedPlatinum) %>%
  tibble()

cat("... dnp db \n")
dnpDb <-
  vroom_read_safe(path = file.path(pathDataInterimTablesAnalysed, "dnp.tsv.gz")) %>%
  data.frame()

inhouseDb <- bind_rows(dnpDb, openDb)

initial_stats <- dbTable %>%
  group_by(database) %>%
  count(name = "initial retrieved unique observations")

cleaned_stats <- inhouseDbMinimal %>%
  group_by(database) %>%
  distinct(
    organismCleaned,
    structureCleanedInchikey,
    referenceCleanedTitle
  ) %>% ## we should decide if 2 or 3D
  count(name = "cleaned documented structure-organism pairs")

final_stats <- openDb %>%
  group_by(database) %>%
  distinct(
    organismCleaned,
    structureCleanedInchikey,
    referenceCleanedTitle
  ) %>% ## we should decide if 2 or 3D
  count(name = "pairs validated for wikidata export")

stats_table <- left_join(initial_stats, cleaned_stats) %>%
  left_join(., final_stats)

dataset <- dataset %>%
  left_join(., stats_table) %>%
  select(
    database,
    `type`,
    `initial retrieved unique observations`,
    `cleaned documented structure-organism pairs`,
    `pairs validated for wikidata export`,
    everything()
  )

pairsOpenDb <- openDb %>%
  filter(!is.na(organismCleaned) &
    !is.na(structureCleaned_inchikey2D)) %>%
  distinct(structureCleaned_inchikey2D,
    organismCleaned,
    .keep_all = TRUE
  )

pairsOutsideDnp <- inhouseDb %>%
  filter(!is.na(organismCleaned) &
    !is.na(structureCleaned_inchikey2D)) %>%
  distinct(structureCleaned_inchikey2D,
    organismCleaned,
    .keep_all = TRUE
  ) %>%
  filter(database != "dnp_1")

pairsFull <- bind_rows(openDb, dnpDb) %>%
  filter(!is.na(organismCleaned) &
    !is.na(structureCleaned_inchikey2D)) %>%
  distinct(structureCleaned_inchikey2D,
    organismCleaned,
    .keep_all = TRUE
  )

pairsDNP <- dnpDb %>%
  distinct(structureCleaned_inchikey2D,
    organismCleaned,
    .keep_all = TRUE
  )

stats <- pairsOutsideDnp %>%
  group_by(database) %>%
  count() %>%
  arrange(desc(n))

# unique
cat("analysing unique organisms per db \n")
## biological taxa
### open NP DB
openDbOrganism <- openDb %>%
  filter(!is.na(organismCleaned)) %>%
  distinct(organismCleaned)

### inhouseDB
inhouseDbOrganism <- inhouseDb %>%
  filter(!is.na(organismCleaned)) %>%
  distinct(organismCleaned)

cat(paste(
  "inhouse:",
  nrow(inhouseDbOrganism),
  "distinct organisms \n",
  sep = " "
))

### DNP
dnpDbOrganism <- dnpDb %>%
  filter(!is.na(organismCleaned)) %>%
  distinct(organismCleaned)

cat(paste("dnp:", nrow(dnpDbOrganism), "distinct organisms \n", sep = " "))

## structures
cat("analysing unique structures per db \n")
### open NP DB
openDbStructure <- openDb %>%
  filter(!is.na(structureCleaned_inchikey2D)) %>%
  distinct(structureCleaned_inchikey2D, .keep_all = TRUE)

### inhouseDB
inhouseDbStructure <- inhouseDb %>%
  filter(!is.na(structureCleaned_inchikey2D)) %>%
  distinct(structureCleaned_inchikey2D, .keep_all = TRUE)

cat(paste(
  "inhouse:",
  nrow(inhouseDbStructure),
  "distinct structures \n",
  sep = " "
))

### DNP
dnpDbStructure <- dnpDb %>%
  filter(!is.na(structureCleaned_inchikey2D)) %>%
  distinct(structureCleaned_inchikey2D, .keep_all = TRUE)

cat(paste("dnp:", nrow(dnpDbStructure), "distinct structures \n", sep = " "))

structuresPerOrganism <- pairsOpenDb %>%
  filter(grepl(pattern = "species", x = organismCleaned_dbTaxoTaxonRanks)) %>%
  distinct(organismCleaned, structureCleaned_inchikey2D) %>%
  group_by(organismCleaned) %>%
  count()

tableStructures <-
  c(
    "only001_structure" = sum(structuresPerOrganism$n == 1),
    "between001and010_structures" = sum(structuresPerOrganism$n >= 1 &
      structuresPerOrganism$n <= 9),
    "between010and100_structures" = sum(structuresPerOrganism$n >= 10 &
      structuresPerOrganism$n <= 99),
    "above100_structures" = sum(structuresPerOrganism$n >= 100)
  ) %>%
  data.frame()
colnames(tableStructures)[1] <- "organisms"

organismsPerStructure <- pairsOpenDb %>%
  filter(grepl(pattern = "species", x = organismCleaned_dbTaxoTaxonRanks)) %>%
  distinct(organismCleaned, structureCleaned_inchikey2D) %>%
  group_by(structureCleaned_inchikey2D) %>%
  count()

tableOrganisms <-
  c(
    "only001_organism" = sum(organismsPerStructure$n == 1),
    "between001and010_organisms" = sum(organismsPerStructure$n >= 1 &
      organismsPerStructure$n <= 9),
    "between010and100_organisms" = sum(organismsPerStructure$n >= 10 &
      organismsPerStructure$n <= 99),
    "above100_organisms" = sum(organismsPerStructure$n >= 100)
  ) %>%
  data.frame()
colnames(tableOrganisms)[1] <- "structures"

# writing tabular stats
## species by kingdom
# cat("analysing species by kingdom \n")
# inhouseSpeciesByKingdom <- openDbTripletsGoldMeta %>%
#   group_by(organismCleaned_dbTaxo_1kingdom) %>%
#   distinct(organismCleaned_dbTaxo_7species, .keep_all = TRUE) %>%
#   count(organismCleaned_dbTaxo_1kingdom) %>%
#   ungroup() %>%
#   mutate(speciesPercent = 100 * n / sum(n)) %>%
#   select(kingdom = organismCleaned_dbTaxo_1kingdom,
#          species = n,
#          speciesPercent) %>%
#   arrange(desc(speciesPercent)) %>%
#   head(10)

## structures by class
# inhouseStructuresByClass <- openDbTripletsGoldMeta %>%
#   group_by(structureCleaned_3class) %>%
#   distinct(structureCleanedInchikey, .keep_all = TRUE) %>%
#   count(structureCleaned_3class) %>%
#   ungroup() %>%
#   mutate(structuresPercent = 100 * n / sum(n)) %>%
#   select(class = structureCleaned_3class,
#          structures = n,
#          structuresPercent) %>%
#   arrange(desc(structuresPercent)) %>%
#   head(10)

## structures by kingdom
# cat("analysing structures by kingdom \n")
# inhouseStructuresByOrganismKingdom <- openDbTripletsGoldMeta %>%
#   group_by(organismCleaned_dbTaxo_1kingdom) %>%
#   distinct(structureCleanedInchikey, organismCleaned_dbTaxo_1kingdom, .keep_all = TRUE) %>%
#   count(organismCleaned_dbTaxo_1kingdom) %>%
#   ungroup() %>%
#   mutate(structuresPercent = 100 * n / sum(n)) %>%
#   select(kingdom = organismCleaned_dbTaxo_1kingdom,
#          structures = n,
#          structuresPercent) %>%
#   arrange(desc(structuresPercent)) %>%
#   head(10)

## unique structures per kingdom
# cat("analysing unique structures by kingdom \n")
# inhouseUniqueStructuresPerKingdom <- openDbTripletsGoldMeta %>%
#   group_by(structureCleanedInchikey) %>%
#   add_count(structureCleanedInchikey) %>%
#   ungroup() %>%
#   filter(n == 1) %>%
#   select(-n) %>%
#   group_by(organismCleaned_dbTaxo_1kingdom) %>%
#   distinct(structureCleanedInchikey, organismCleaned_dbTaxo_1kingdom, .keep_all = TRUE) %>%
#   count(organismCleaned_dbTaxo_1kingdom) %>%
#   ungroup() %>%
#   select(kingdom = organismCleaned_dbTaxo_1kingdom,
#          specificStructures = n) %>%
#   arrange(desc(specificStructures)) %>%
#   head(10)

## structures by kingdom
# cat("joining \n")
# inhouseStructuresByKingdom <-
#   full_join(inhouseSpeciesByKingdom,
#             inhouseStructuresByOrganismKingdom) %>%
#   mutate(strucuresPerSpecies = structures / species)
#
# inhouseStructuresByKingdom <-
#   full_join(inhouseStructuresByKingdom,
#             inhouseUniqueStructuresPerKingdom) %>%
#   mutate(structuresSpecificity = 100 * specificStructures / structures) %>%
#   select(1, 2, 3, 4, 7, 8, 5, 6)

## unique structures per species
# cat("analysing unique structures by species \n")
# inhouseUniqueStructuresPerSpecies <- openDbTripletsGoldMeta %>%
#   group_by(structureCleanedInchikey) %>%
#   add_count(structureCleanedInchikey) %>%
#   ungroup() %>%
#   filter(n == 1) %>%
#   select(-n) %>%
#   group_by(organismCleaned_dbTaxo_7species) %>%
#   distinct(structureCleanedInchikey, organismCleaned_dbTaxo_7species, .keep_all = TRUE) %>%
#   count(organismCleaned_dbTaxo_7species) %>%
#   ungroup() %>%
#   select(species = organismCleaned_dbTaxo_7species,
#          specificStructures = n) %>%
#   arrange(desc(specificStructures)) %>%
#   filter(!is.na(species)) %>%
#   head(10)

## widespread metabolites
# cat("analysing widespread metabolites \n")
# openDbWidespread <- inhouseDbTripletsGoldMeta %>%
#   group_by(structureCleanedInchikey) %>%
#   filter(!is.na(organismCleaned_dbTaxo_1kingdom)) %>%
#   distinct(organismCleaned_dbTaxo_1kingdom, .keep_all = TRUE) %>%
#   add_count() %>%
#   ungroup() %>%
#   filter(n >= 6) %>%
#   arrange(structureCleanedInchikey)

## word(species,1) != genus
# cat("analysing mismatched genera \n")
# mismatchedGenera <- inhouseDbOrganism %>%
#   filter(word(organism_7_species, 1) != organism_6_genus)

## redundancy table
# heavy process not bringing much for now
# cat("analysing redundant entries \n")
# redundancydf  <- inhouseDb %>%
#   filter(!is.na(organismLowestTaxon) &
#            !is.na(inchikeySanitized)) %>%
#   distinct(organismLowestTaxon,
#            inchikeySanitized,
#            database,
#            .keep_all = TRUE) %>%
#   group_by(organismLowestTaxon, inchikeySanitized) %>%
#   add_count() %>%
#   ungroup() %>%
#   filter(n >= 5) %>%
#   select(
#     database,
#     structureOriginal_inchi,
#     structureOriginal_smiles,
#     structureOriginal_nominal,
#     organismOriginal,
#     organismLowestTaxon,
#     inchikeySanitized,
#     n
#   )

# exporting
# # stats
# ## structures by kingdom

# write.table(
#   x = inhouseStructuresByKingdom,
#   file = pathDataInterimTablesAnalysedStructuresByKingdom,
#   row.names = FALSE,
#   quote = FALSE,
#   sep = "\t",
#   fileEncoding = "UTF-8"
# )

### unique structures per species
# write.table(
#   x = inhouseUniqueStructuresPerSpecies,
#   file = pathDataInterimTablesAnalysedUniqueStructuresBySpecies,
#   row.names = FALSE,
#   quote = FALSE,
#   sep = "\t",
#   fileEncoding = "UTF-8"
# )

### widespread metabolites
# write.table(
#   x = openDbWidespread,
#   file = pathDataInterimTablesAnalysedWidespreadStructures,
#   row.names = FALSE,
#   quote = FALSE,
#   sep = "\t",
#   fileEncoding = "UTF-8"
# )

### mismatched genera
# write.table(
#   x = mismatchedGenera,
#   file = pathDataInterimTablesAnalysedMismatchedGenera,
#   row.names = TRUE,
#   quote = FALSE,
#   sep = "\t",
#   fileEncoding = "UTF-8"
# )

### redundancy table
# write.table(
#   x = redundancydf,
#   file = pathDataInterimTablesAnalysedRedundancyTable,
#   row.names = FALSE,
#   quote = FALSE,
#   sep = "\t",
#   fileEncoding = "UTF-8"
# )

if (mode == "FULL" | mode == "full") {
  write.table(
    x = dataset,
    file = "../docs/dataset.tsv",
    row.names = FALSE,
    quote = FALSE,
    sep = "\t",
    fileEncoding = "UTF-8"
  )

  cat(
    paste(
      nrow(pairsOpenDb),
      "unique (2D) referenced structure-organism pairs. \n",
      "They consist of \n",
      nrow(openDbStructure),
      "unique (2D) curated structures and \n",
      nrow(openDbOrganism),
      "unique organisms,\n",
      "originating from \n",
      nrow(pairsOpenDb %>% distinct(database)),
      "initial open databases. \n",
      "\n",
      "Among structures, \n",
      tableOrganisms[1, 1],
      "are present in only 1 organism, \n",
      tableOrganisms[2, 1],
      "are present in between 1 and 10 organisms, \n",
      tableOrganisms[3, 1],
      "are present in between 10 and 100 organisms, \n",
      tableOrganisms[4, 1],
      "are present in more than 100 organisms. \n",
      "\n",
      "Among organisms, \n",
      tableStructures[1, 1],
      "contain only 1 structure, \n",
      tableStructures[2, 1],
      "contain between 1 and 10 structures, \n",
      tableStructures[3, 1],
      "contain between 10 and 100 structures, \n",
      tableStructures[4, 1],
      "contain more than 100 structures. \n",
      sep = " "
    ),
    file = "../docs/metrics.adoc"
  )
}

end <- Sys.time()

cat("Script finished in", format(end - start), "\n")
