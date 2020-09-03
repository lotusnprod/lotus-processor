# title: "Unique pairs compileR"

# loading functions
source("functions.R")
source("paths.R")

# loading files
print("loading db, if running fullmode, this may take a while \n")

## inhouseDb
inhouseDbMinimal <- read_delim(
  file = gzfile(pathDataInterimTablesCuratedTable),
  col_types = cols(.default = "c"),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

openDbMinimalFiltered <- inhouseDbMinimal %>%
  filter(!is.na(organismCleaned) &
           !is.na(structureCleanedInchikey3D)) %>%
  filter(
    !is.na(referenceCleanedDoi) |
      !is.na(referenceCleanedPmcid) |
      !is.na(referenceCleanedPmid)
  ) %>%
  distinct(
    database,
    organismCleaned,
    organismCleaned_dbTaxo,
    organismCleaned_dbTaxoTaxonIds,
    organismCleaned_dbTaxoTaxonRanks,
    organismCleaned_dbTaxoTaxonomy,
    structureCleanedSmiles,
    structureCleanedInchi,
    structureCleanedInchikey3D,
    referenceCleanedDoi,
    referenceCleanedPmcid,
    referenceCleanedPmid,
    .keep_all = TRUE
  )

dnpDb <- inhouseDbMinimal %>%
  filter(database == "dnp_1") %>%
  filter(!is.na(organismCleaned) &
           !is.na(structureCleanedInchikey3D)) %>%
  distinct(
    database,
    organismCleaned,
    organismCleaned_dbTaxo,
    organismCleaned_dbTaxoTaxonIds,
    organismCleaned_dbTaxoTaxonRanks,
    organismCleaned_dbTaxoTaxonomy,
    structureCleanedSmiles,
    structureCleanedInchi,
    structureCleanedInchikey3D,
    .keep_all = TRUE
  ) %>%
  mutate(
    referenceCleanedDoi = NA,
    referenceCleanedPmcid = NA,
    referenceCleanedPmid = NA
  ) %>% mutate_all(as.character)

rm(inhouseDbMinimal)

openDbMinimalFilteredRef <- openDbMinimalFiltered %>%
  distinct(organismCleaned,
           referenceCleanedDoi,
           referenceCleanedPmcid,
           referenceCleanedPmid)


## reference metadata
referenceTableFull <- read_delim(
  file = gzfile(pathDataInterimDictionariesReferenceOrganismDictionary),
  col_types = cols(.default = "c"),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

# selecting reference metadata
referenceMetadataFiltered <- referenceTableFull %>%
  filter(
    !is.na(referenceCleanedDoi) |
      !is.na(referenceCleanedPmcid) |
      !is.na(referenceCleanedPmid)
  ) %>%
  select(
    -referenceCleaned_title,
    -referenceCleaned_journal,
    -referenceCleaned_date,
    -referenceCleaned_author
  ) %>%
  mutate(
    referenceCleaned_score_crossref = as.integer(referenceCleaned_score_crossref),
    referenceCleaned_score_distance = as.integer(referenceCleaned_score_distance),
    referenceCleaned_score_titleOrganism = as.integer(referenceCleaned_score_titleOrganism)
  )

# joining inhousedb minimal and reference metadata
openDbRef <- left_join(openDbMinimalFilteredRef,
                       referenceMetadataFiltered) %>%
  filter(
    !is.na(referenceCleanedPmid) |
      !is.na(referenceCleanedPmcid) |
      !is.na(referenceCleanedDoi)
  ) %>%
  arrange(desc(referenceCleanedPmid)) %>%
  arrange(desc(referenceCleanedPmcid)) %>%
  arrange(desc(referenceCleanedDoi)) %>%
  arrange(desc(referenceCleaned_score_crossref)) %>%
  arrange(referenceCleaned_score_distance) %>%
  arrange(desc(referenceCleaned_score_titleOrganism)) #very important to keep references

openDb <- right_join(openDbRef, openDbMinimalFiltered) %>%
  distinct(
    database,
    organismCleaned,
    organismCleaned_dbTaxo,
    organismCleaned_dbTaxoTaxonIds,
    organismCleaned_dbTaxoTaxonRanks,
    organismCleaned_dbTaxoTaxonomy,
    structureCleanedInchi,
    structureCleanedInchikey3D,
    structureCleanedSmiles,
    referenceCleanedDoi,
    referenceCleanedPmcid,
    referenceCleanedPmid,
    .keep_all = TRUE
  ) %>%
  select(
    database,
    organismOriginal,
    structureType,
    structureValue,
    referenceType,
    referenceValue,
    organismCleaned,
    organismCleaned_dbTaxo,
    organismCleaned_dbTaxoTaxonIds,
    organismCleaned_dbTaxoTaxonRanks,
    organismCleaned_dbTaxoTaxonomy,
    structureCleanedInchi,
    structureCleanedInchikey3D,
    structureCleanedSmiles,
    referenceCleanedDoi,
    referenceCleanedPmcid,
    referenceCleanedPmid
  )

inhouseDb <- bind_rows(dnpDb, openDb)

pairsOpenDb <- openDb %>%
  filter(!is.na(organismCleaned) &
           !is.na(structureCleanedInchikey3D)) %>%
  distinct(structureCleanedInchikey3D,
           organismCleaned,
           .keep_all = TRUE)

pairsOutsideDnp <- inhouseDb %>%
  filter(!is.na(organismCleaned) &
           !is.na(structureCleanedInchikey3D)) %>%
  distinct(structureCleanedInchikey3D,
           organismCleaned,
           .keep_all = TRUE) %>%
  filter(database != "dnp_1")

pairsFull <- bind_rows(openDb, dnpDb) %>%
  filter(!is.na(organismCleaned) &
           !is.na(structureCleanedInchikey3D)) %>%
  distinct(structureCleanedInchikey3D,
           organismCleaned,
           .keep_all = TRUE)

pairsDNP <- dnpDb %>%
  distinct(structureCleanedInchikey3D,
           organismCleaned,
           .keep_all = TRUE)

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

cat(paste("open:", nrow(openDbOrganism), "distinct organisms \n", sep = " "))

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
  filter(!is.na(structureCleanedInchikey3D)) %>%
  distinct(structureCleanedInchikey3D, .keep_all = TRUE)

cat(paste("open:", nrow(openDbStructure), "distinct structures \n", sep = " "))

### inhouseDB
inhouseDbStructure <- inhouseDb %>%
  filter(!is.na(structureCleanedInchikey3D)) %>%
  distinct(structureCleanedInchikey3D, .keep_all = TRUE)

cat(paste(
  "inhouse:",
  nrow(inhouseDbStructure),
  "distinct structures \n",
  sep = " "
))

### DNP
dnpDbStructure <- dnpDb %>%
  filter(!is.na(structureCleanedInchikey3D)) %>%
  distinct(structureCleanedInchikey3D, .keep_all = TRUE)

cat(paste("dnp:", nrow(dnpDbStructure), "distinct structures \n", sep = " "))

# writing tabular stats
## species by kingdom
# cat("analysing species by kingdom \n")
# inhouseSpeciesByKingdom <- inhouseDbPairs %>%
#   group_by(organism_1_kingdom) %>%
#   distinct(organism_7_species, .keep_all = TRUE) %>%
#   count(organism_1_kingdom) %>%
#   ungroup() %>%
#   mutate(speciesPercent = 100 * n / sum(n)) %>%
#   select(kingdom = organism_1_kingdom,
#          species = n,
#          speciesPercent) %>%
#   arrange(desc(speciesPercent)) %>%
#   head(10)

## structures by class
# inhouseStructuresByClass <- inhouseDbPairs %>%
#   group_by(structure_03_class) %>%
#   distinct(inchikeySanitized, .keep_all = TRUE) %>%
#   count(structure_03_class) %>%
#   ungroup() %>%
#   mutate(structuresPercent = 100 * n / sum(n)) %>%
#   select(class = structure_03_class,
#          structures = n,
#          structuresPercent) %>%
#   arrange(desc(structuresPercent)) %>%
#   head(10)

## structures by kingdom
# cat("analysing structures by kingdom \n")
# inhouseStructuresByOrganismKingdom <- inhouseDbPairs %>%
#   group_by(organism_1_kingdom) %>%
#   distinct(inchikeySanitized, organism_1_kingdom, .keep_all = TRUE) %>%
#   count(organism_1_kingdom) %>%
#   ungroup() %>%
#   mutate(structuresPercent = 100 * n / sum(n)) %>%
#   select(kingdom = organism_1_kingdom,
#          structures = n,
#          structuresPercent) %>%
#   arrange(desc(structuresPercent)) %>%
#   head(10)

## unique structures per kingdom
# cat("analysing unique structures by kingdom \n")
# inhouseUniqueStructuresPerKingdom <- inhouseDbPairs %>%
#   group_by(inchikeySanitized) %>%
#   add_count(inchikeySanitized) %>%
#   ungroup() %>%
#   filter(n == 1) %>%
#   select(-n) %>%
#   group_by(organism_1_kingdom) %>%
#   distinct(inchikeySanitized, organism_1_kingdom, .keep_all = TRUE) %>%
#   count(organism_1_kingdom) %>%
#   ungroup() %>%
#   select(kingdom = organism_1_kingdom,
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

##unique structures per species
# cat("analysing unique structures by species \n")
# inhouseUniqueStructuresPerSpecies <- inhouseDbPairs %>%
#   group_by(inchikeySanitized) %>%
#   add_count(inchikeySanitized) %>%
#   ungroup() %>%
#   filter(n == 1) %>%
#   select(-n) %>%
#   group_by(organism_7_species) %>%
#   distinct(inchikeySanitized, organism_7_species, .keep_all = TRUE) %>%
#   count(organism_7_species) %>%
#   ungroup() %>%
#   select(species = organism_7_species,
#          specificStructures = n) %>%
#   arrange(desc(specificStructures)) %>%
#   filter(!is.na(species)) %>%
#   head(10)

##widespread metabolites
# cat("analysing widespread metabolites \n")
# openDbWidespread <- openDbPairs %>%
#   group_by(inchikeySanitized) %>%
#   filter(!is.na(organism_1_kingdom)) %>%
#   distinct(organism_1_kingdom, .keep_all = TRUE) %>%
#   add_count() %>%
#   ungroup() %>%
#   filter(n >= 6) %>%
#   arrange(inchikeySanitized)

##word(species,1) != genus
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

#exporting
cat("exporting, may take a while if running full mode \n")
## creating directories if they do not exist
ifelse(
  !dir.exists(pathDataInterimTablesAnalysed),
  dir.create(pathDataInterimTablesAnalysed),
  FALSE
)

##open
write.table(
  x = openDb,
  file = gzfile(
    description = pathDataInterimTablesAnalysedOpenDbTriplets,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

##inhouse
write.table(
  x = inhouseDb,
  file = gzfile(
    description = pathDataInterimTablesAnalysedInhouseDbTriplets,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

##dnp
write.table(
  x = dnpDb,
  file = gzfile(
    description = pathDataInterimTablesAnalysedDnpDbTriplets,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

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

# ## unique structures per species
# write.table(
#   x = inhouseUniqueStructuresPerSpecies,
#   file = pathDataInterimTablesAnalysedUniqueStructuresBySpecies,
#   row.names = FALSE,
#   quote = FALSE,
#   sep = "\t",
#   fileEncoding = "UTF-8"
# )

###widespread metabolites
# write.table(
#   x = openDbWidespread,
#   file = pathDataInterimTablesAnalysedWidespreadStructures,
#   row.names = FALSE,
#   quote = FALSE,
#   sep = "\t",
#   fileEncoding = "UTF-8"
# )

###mismatched genera
# write.table(
#   x = mismatchedGenera,
#   file = pathDataInterimTablesAnalysedMismatchedGenera,
#   row.names = TRUE,
#   quote = FALSE,
#   sep = "\t",
#   fileEncoding = "UTF-8"
# )

###redundancy table
# write.table(
#   x = redundancydf,
#   file = pathDataInterimTablesAnalysedRedundancyTable,
#   row.names = FALSE,
#   quote = FALSE,
#   sep = "\t",
#   fileEncoding = "UTF-8"
# )
