cat("This script outputs some metrics related to the DB \n")

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")

cat("... libraries \n")
library(tidyverse)
source("r/vroom_safe.R")

cat("loading ... \n")
cat("... validated db, if running fullmode, this may take a while \n")
openDb <-
  vroom_read_safe(path = pathDataInterimTablesAnalysedPlatinum) %>%
  data.frame()

cat("... dnp db \n")
dnpDb <-
  vroom_read_safe(path = file.path(pathDataInterimTablesAnalysed, "dnp.tsv.gz")) %>%
  data.frame()

inhouseDb <- bind_rows(dnpDb, openDb)

pairsOpenDb <- openDb %>%
  filter(!is.na(organismCleaned) &
    !is.na(structureCleanedInchikey2D)) %>%
  distinct(structureCleanedInchikey2D,
    organismCleaned,
    .keep_all = TRUE
  )

pairsOutsideDnp <- inhouseDb %>%
  filter(!is.na(organismCleaned) &
    !is.na(structureCleanedInchikey2D)) %>%
  distinct(structureCleanedInchikey2D,
    organismCleaned,
    .keep_all = TRUE
  ) %>%
  filter(database != "dnp_1")

pairsFull <- bind_rows(openDb, dnpDb) %>%
  filter(!is.na(organismCleaned) &
    !is.na(structureCleanedInchikey2D)) %>%
  distinct(structureCleanedInchikey2D,
    organismCleaned,
    .keep_all = TRUE
  )

pairsDNP <- dnpDb %>%
  distinct(structureCleanedInchikey2D,
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
  filter(!is.na(structureCleanedInchikey2D)) %>%
  distinct(structureCleanedInchikey2D, .keep_all = TRUE)

### inhouseDB
inhouseDbStructure <- inhouseDb %>%
  filter(!is.na(structureCleanedInchikey2D)) %>%
  distinct(structureCleanedInchikey2D, .keep_all = TRUE)

cat(paste(
  "inhouse:",
  nrow(inhouseDbStructure),
  "distinct structures \n",
  sep = " "
))

### DNP
dnpDbStructure <- dnpDb %>%
  filter(!is.na(structureCleanedInchikey2D)) %>%
  distinct(structureCleanedInchikey2D, .keep_all = TRUE)

cat(paste("dnp:", nrow(dnpDbStructure), "distinct structures \n", sep = " "))

structuresPerOrganism <- pairsOpenDb %>%
  filter(grepl(pattern = "species", x = organismCleaned_dbTaxoTaxonRanks)) %>%
  distinct(organismCleaned, structureCleanedInchikey2D) %>%
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
  distinct(organismCleaned, structureCleanedInchikey2D) %>%
  group_by(structureCleanedInchikey2D) %>%
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
#   distinct(structureCleanedInchikey3D, .keep_all = TRUE) %>%
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
#   distinct(structureCleanedInchikey3D, organismCleaned_dbTaxo_1kingdom, .keep_all = TRUE) %>%
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
#   group_by(structureCleanedInchikey3D) %>%
#   add_count(structureCleanedInchikey3D) %>%
#   ungroup() %>%
#   filter(n == 1) %>%
#   select(-n) %>%
#   group_by(organismCleaned_dbTaxo_1kingdom) %>%
#   distinct(structureCleanedInchikey3D, organismCleaned_dbTaxo_1kingdom, .keep_all = TRUE) %>%
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
#   group_by(structureCleanedInchikey3D) %>%
#   add_count(structureCleanedInchikey3D) %>%
#   ungroup() %>%
#   filter(n == 1) %>%
#   select(-n) %>%
#   group_by(organismCleaned_dbTaxo_7species) %>%
#   distinct(structureCleanedInchikey3D, organismCleaned_dbTaxo_7species, .keep_all = TRUE) %>%
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
#   group_by(structureCleanedInchikey3D) %>%
#   filter(!is.na(organismCleaned_dbTaxo_1kingdom)) %>%
#   distinct(organismCleaned_dbTaxo_1kingdom, .keep_all = TRUE) %>%
#   add_count() %>%
#   ungroup() %>%
#   filter(n >= 6) %>%
#   arrange(structureCleanedInchikey3D)

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

if (mode == "FULL") {
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
