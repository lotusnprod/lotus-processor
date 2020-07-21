# title: "Unique pairs compileR"

# loading functions
source("functions.R")
source("paths.R")

# loading files
print(x = "loading db, if running fullmode, this may take a while")

## inhouseDb
inhouseDb <- read_delim(
  file = gzfile(pathDataInterimTablesCuratedTable),
  col_types = cols(.default = "c"),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate(referenceCleanedTranslationScore = as.integer(referenceCleanedTranslationScore)) %>%
  arrange(desc(referenceOriginalExternal)) %>%
  arrange(desc(referenceCleanedTranslationScore)) %>%
  arrange(desc(referenceCleanedDoi)) %>% #very important to keep references
  data.frame()

## openDB
openDb <- inhouseDb %>%
  filter(database != "dnp_1")

## DNP
dnpDb <- inhouseDb %>%
  filter(database == "dnp_1")

pairsOutsideDnp <- rbind(dnpDb, openDb) %>%
  filter(!is.na(organismLowestTaxon) &
           !is.na(inchikeySanitized)) %>%
  distinct(inchikeySanitized,
           organismLowestTaxon,
           .keep_all = TRUE) %>%
  filter(database != "dnp_1")

pairsFull <- rbind(openDb, dnpDb) %>%
  filter(!is.na(organismLowestTaxon) &
           !is.na(inchikeySanitized)) %>%
  distinct(inchikeySanitized,
           organismLowestTaxon,
           .keep_all = TRUE)

# warning: keeping only one ref per triplet

tripletsOutsideDnpStrict <- pairsOutsideDnp %>%
  distinct(
    inchikeySanitized,
    organismLowestTaxon,
    referenceCleanedDoi,
    referenceOriginalExternal,
    referenceCleanedTitle,
    .keep_all = TRUE
  ) %>%
  filter(!is.na(referenceCleanedDoi) &
           referenceCleanedTranslationScore == 100)

tripletsOverlapDnpStrict <- pairsFull %>%
  distinct(
    inchikeySanitized,
    organismLowestTaxon,
    referenceCleanedDoi,
    referenceOriginalExternal,
    referenceCleanedTitle,
    .keep_all = TRUE
  ) %>%
  filter(!is.na(referenceCleanedDoi) &
           referenceCleanedTranslationScore == 100)

tripletsWithDnpStrict <- pairsFull %>%
  distinct(
    inchikeySanitized,
    organismLowestTaxon,
    referenceCleanedDoi,
    referenceOriginalExternal,
    referenceCleanedTitle,
    .keep_all = TRUE
  ) %>%
  filter(
    !is.na(referenceCleanedDoi) &
      referenceCleanedTranslationScore == 100 |
      referenceOriginalExternal == "DNP"
  )

tripletsDNP <- dnpDb %>%
  distinct(inchikeySanitized,
           organismLowestTaxon,
           .keep_all = TRUE)

tripletsOverlapDnpMedium <- pairsFull %>%
  distinct(
    inchikeySanitized,
    organismLowestTaxon,
    referenceCleanedDoi,
    referenceOriginalExternal,
    referenceCleanedTitle,
    .keep_all = TRUE
  ) %>%
  filter(
    !is.na(referenceCleanedDoi) |
      !is.na(referenceCleanedTitle) |
      !is.na(referenceOriginalExternal)
  ) %>%
  filter(referenceOriginalExternal != "DNP" |
           referenceCleanedTranslationScore == 100)

stats <- pairsOutsideDnp %>%
  group_by(database) %>%
  count()

# unique
print(x = "analysing unique organisms per db")
## biological taxa
### open NP DB
print(x = "open")
openDbOrganism <- openDb %>%
  filter(!is.na(organismLowestTaxon)) %>%
  distinct(organismLowestTaxon)

### inhouseDB
print(x = "inhouse")
inhouseDbOrganism <- inhouseDb %>%
  filter(!is.na(organismLowestTaxon)) %>%
  
  distinct(organismLowestTaxon)

### DNP
print(x = "dnp")
dnpDbOrganism <- dnpDb %>%
  filter(!is.na(organismLowestTaxon)) %>%
  
  distinct(organismLowestTaxon)

## structures
print(x = "analysing unique structures per db")
### open NP DB
print(x = "open")
openDbStructure <- openDb %>%
  filter(!is.na(inchikeySanitized)) %>%
  distinct(inchikeySanitized, .keep_all = TRUE)

### inhouseDB
print(x = "inhouse")
inhouseDbStructure <- inhouseDb %>%
  filter(!is.na(inchikeySanitized)) %>%
  distinct(inchikeySanitized, .keep_all = TRUE)

### DNP
print(x = "dnp")
dnpDbStructure <- dnpDb %>%
  filter(!is.na(inchikeySanitized)) %>%
  distinct(inchikeySanitized, .keep_all = TRUE)

## references
### open NP DB
openDbReference <- openDb %>%
  filter(!is.na(referenceCleanedDoi)) %>%
  distinct(referenceCleanedDoi, .keep_all = TRUE)

### inhouseDB
inhouseDbReference <- inhouseDb %>%
  filter(!is.na(referenceOriginalExternal) |
           !is.na(referenceCleanedDoi)) %>%
  distinct(referenceOriginalExternal, referenceCleanedDoi, .keep_all = TRUE)

### DNP
dnpDbReference <- dnpDb %>%
  filter(!is.na(referenceOriginalExternal) |
           !is.na(referenceCleanedDoi)) %>%
  distinct(referenceOriginalExternal, referenceCleanedDoi, .keep_all = TRUE)

## triplets
print(x = "analysing triplets, this may take a while")
print(x = "open")
###open NP DB
openDbTriplets <- openDb %>%
  filter(!is.na(referenceCleanedDoi) |
           !is.na(referenceOriginalExternal)) %>%
  filter(!is.na(inchikeySanitized) &
           !is.na(organismLowestTaxon)) %>%
  distinct(
    inchikeySanitized,
    referenceCleanedDoi,
    referenceOriginalExternal,
    organismLowestTaxon,
    organismTaxonId,
    # here could be modif
    .keep_all = TRUE
  )

print(x = "inhouse")
### inhouseDB
inhouseDbTriplets <- inhouseDb %>%
  filter(!is.na(referenceCleanedDoi) |
           !is.na(referenceOriginalExternal)) %>%
  filter(!is.na(inchikeySanitized) &
           !is.na(organismLowestTaxon)) %>%
  distinct(
    inchikeySanitized,
    referenceCleanedDoi,
    referenceOriginalExternal,
    organismLowestTaxon,
    organismTaxonId,
    # here could be modif
    .keep_all = TRUE
  )

print(x = "dnp")
### DNP
dnpDbTriplets <- dnpDb %>%
  filter(!is.na(referenceCleanedDoi) |
           !is.na(referenceOriginalExternal)) %>%
  filter(!is.na(inchikeySanitized) &
           !is.na(organismLowestTaxon)) %>%
  distinct(
    inchikeySanitized,
    referenceCleanedDoi,
    referenceOriginalExternal,
    organismLowestTaxon,
    organismTaxonId,
    # here could be modif
    .keep_all = TRUE
  )

## pairs
print(x = "analysing pairs, this should be faster")
### open NP DB
print(x = "open")
openDbPairs <- openDbTriplets %>%
  filter(referenceCleanedTranslationScore >= 30 |
           !is.na(referenceOriginalExternal)) %>%
  distinct(inchikeySanitized, organismLowestTaxon, .keep_all = TRUE)

### inhouseDB
print(x = "inhouse")
inhouseDbPairs <- inhouseDbTriplets %>%
  filter(referenceCleanedTranslationScore >= 30 |
           !is.na(referenceOriginalExternal)) %>%
  distinct(inchikeySanitized, organismLowestTaxon, .keep_all = TRUE)

### DNP
print(x = "dnp")
dnpDbPairs <- dnpDbTriplets %>%
  filter(referenceCleanedTranslationScore >= 30 |
           !is.na(referenceOriginalExternal)) %>%
  distinct(inchikeySanitized, organismLowestTaxon, .keep_all = TRUE)

# writing tabular stats
## species by kingdom
# print(x = "analysing species by kingdom")
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
# print(x = "analysing structures by kingdom")
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
# print(x = "analysing unique structures by kingdom")
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
# print(x = "joining")
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
# print(x = "analysing unique structures by species")
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
# print(x = "analysing widespread metabolites")
# openDbWidespread <- openDbPairs %>%
#   group_by(inchikeySanitized) %>%
#   filter(!is.na(organism_1_kingdom)) %>%
#   distinct(organism_1_kingdom, .keep_all = TRUE) %>%
#   add_count() %>%
#   ungroup() %>%
#   filter(n >= 6) %>%
#   arrange(inchikeySanitized)

##word(species,1) != genus
# print(x = "analysing mismatched genera")
# mismatchedGenera <- inhouseDbOrganism %>%
#   filter(word(organism_7_species, 1) != organism_6_genus)

##redundancy table
print(x = "analysing redundant entries")
redundancydf  <- inhouseDb %>%
  filter(!is.na(organismLowestTaxon) &
           !is.na(inchikeySanitized)) %>%
  distinct(organismLowestTaxon,
           inchikeySanitized,
           database,
           .keep_all = TRUE) %>%
  group_by(organismLowestTaxon, inchikeySanitized) %>%
  add_count() %>%
  ungroup() %>%
  filter(n >= 5) %>%
  select(
    database,
    structureOriginalInchi,
    structureOriginalSmiles,
    structureOriginalNominal,
    organismOriginal,
    organismLowestTaxon,
    inchikeySanitized,
    n
  )

#exporting
print(x = "exporting, may take a while if running full mode")
## creating directories if they do not exist
ifelse(
  !dir.exists(pathDataInterimTablesAnalysed),
  dir.create(pathDataInterimTablesAnalysed),
  FALSE
)

##open
write.table(
  x = openDbTriplets,
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
  x = inhouseDbTriplets,
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
  x = dnpDbTriplets,
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

##stats
###structures by kingdom
# write.table(
#   x = inhouseStructuresByKingdom,
#   file = pathDataInterimTablesAnalysedStructuresByKingdom,
#   row.names = FALSE,
#   quote = FALSE,
#   sep = "\t",
#   fileEncoding = "UTF-8"
# )

###unique structures per species
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
write.table(
  x = redundancydf,
  file = pathDataInterimTablesAnalysedRedundancyTable,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
