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
  file = "../docs/dataset.csv",
  delim = ","
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

cleaned_stats_3D <- inhouseDbMinimal %>%
  group_by(database) %>%
  distinct(
    organismCleaned,
    structureCleanedInchikey,
    referenceCleanedTitle
  ) %>%
  count(name = "cleaned documented structure-organism pairs")

final_stats_3D <- openDb %>%
  group_by(database) %>%
  distinct(
    organismCleaned,
    structureCleanedInchikey,
    referenceCleanedTitle
  ) %>%
  count(name = "pairs validated for wikidata export")

stats_table <- left_join(initial_stats, cleaned_stats_3D) %>%
  left_join(., final_stats_3D)

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

pairsOpenDb_3D <- openDb %>%
  filter(!is.na(organismCleaned) &
    !is.na(structureCleanedInchikey)) %>%
  distinct(structureCleanedInchikey,
    organismCleaned,
    .keep_all = TRUE
  )

pairsOpenDb_2D <- openDb %>%
  filter(!is.na(organismCleaned) &
    !is.na(structureCleaned_inchikey2D)) %>%
  distinct(structureCleaned_inchikey2D,
    organismCleaned,
    .keep_all = TRUE
  )

pairsOutsideDnp_3D <- inhouseDb %>%
  filter(!is.na(organismCleaned) &
    !is.na(structureCleanedInchikey)) %>%
  distinct(structureCleanedInchikey,
    organismCleaned,
    .keep_all = TRUE
  ) %>%
  filter(database != "dnp")

pairsOutsideDnp_2D <- inhouseDb %>%
  filter(!is.na(organismCleaned) &
    !is.na(structureCleaned_inchikey2D)) %>%
  distinct(structureCleaned_inchikey2D,
    organismCleaned,
    .keep_all = TRUE
  ) %>%
  filter(database != "dnp")

pairsFull_3D <- bind_rows(openDb, dnpDb) %>%
  filter(!is.na(organismCleaned) &
    !is.na(structureCleanedInchikey)) %>%
  distinct(structureCleanedInchikey,
    organismCleaned,
    .keep_all = TRUE
  )

pairsFull_2D <- bind_rows(openDb, dnpDb) %>%
  filter(!is.na(organismCleaned) &
    !is.na(structureCleaned_inchikey2D)) %>%
  distinct(structureCleaned_inchikey2D,
    organismCleaned,
    .keep_all = TRUE
  )

pairsDNP_3D <- dnpDb %>%
  distinct(structureCleanedInchikey,
    organismCleaned,
    .keep_all = TRUE
  )

pairsDNP_2D <- dnpDb %>%
  distinct(structureCleaned_inchikey2D,
    organismCleaned,
    .keep_all = TRUE
  )

stats_3D <- pairsOutsideDnp_3D %>%
  group_by(database) %>%
  count() %>%
  arrange(desc(n))

stats_2D <- pairsOutsideDnp_2D %>%
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
cat("analysing unique structures (2D) per db \n")
### open NP DB
openDbStructure_3D <- openDb %>%
  filter(!is.na(structureCleanedInchikey)) %>%
  distinct(structureCleanedInchikey, .keep_all = TRUE)

### inhouseDB
inhouseDbStructure_3D <- inhouseDb %>%
  filter(!is.na(structureCleanedInchikey)) %>%
  distinct(structureCleanedInchikey, .keep_all = TRUE)

cat(paste(
  "inhouse:",
  nrow(inhouseDbStructure_3D),
  "distinct structures (3D) \n",
  sep = " "
))

cat("analysing unique structures (2D) per db \n")
### open NP DB
openDbStructure_2D <- openDb %>%
  filter(!is.na(structureCleaned_inchikey2D)) %>%
  distinct(structureCleaned_inchikey2D, .keep_all = TRUE)

### inhouseDB
inhouseDbStructure_2D <- inhouseDb %>%
  filter(!is.na(structureCleaned_inchikey2D)) %>%
  distinct(structureCleaned_inchikey2D, .keep_all = TRUE)

cat(paste(
  "inhouse:",
  nrow(inhouseDbStructure_2D),
  "distinct structures (2D) \n",
  sep = " "
))

### DNP
dnpDbStructure_3D <- dnpDb %>%
  filter(!is.na(structureCleanedInchikey)) %>%
  distinct(structureCleanedInchikey, .keep_all = TRUE)

cat(paste(
  "dnp:",
  nrow(dnpDbStructure_3D),
  "distinct structures (3D) \n",
  sep = " "
))

dnpDbStructure_2D <- dnpDb %>%
  filter(!is.na(structureCleaned_inchikey2D)) %>%
  distinct(structureCleaned_inchikey2D, .keep_all = TRUE)

cat(paste(
  "dnp:",
  nrow(dnpDbStructure_2D),
  "distinct structures (2D) \n",
  sep = " "
))

structuresPerOrganism_3D <- pairsOpenDb_3D %>%
  filter(grepl(pattern = "species", x = organismCleaned_dbTaxoTaxonRanks)) %>%
  distinct(organismCleaned, structureCleanedInchikey) %>%
  group_by(organismCleaned) %>%
  count()

structuresPerOrganism_2D <- pairsOpenDb_2D %>%
  filter(grepl(pattern = "species", x = organismCleaned_dbTaxoTaxonRanks)) %>%
  distinct(organismCleaned, structureCleaned_inchikey2D) %>%
  group_by(organismCleaned) %>%
  count()

tableStructures_3D <-
  c(
    "only001_structure" = sum(structuresPerOrganism_3D$n == 1),
    "between001and010_structures" = sum(
      structuresPerOrganism_3D$n >= 1 &
        structuresPerOrganism_3D$n <= 9
    ),
    "between010and100_structures" = sum(
      structuresPerOrganism_3D$n >= 10 &
        structuresPerOrganism_3D$n <= 99
    ),
    "above100_structures" = sum(structuresPerOrganism_3D$n >= 100)
  ) %>%
  data.frame()
colnames(tableStructures_3D)[1] <- "organisms"

tableStructures_2D <-
  c(
    "only001_structure" = sum(structuresPerOrganism_2D$n == 1),
    "between001and010_structures" = sum(
      structuresPerOrganism_2D$n >= 1 &
        structuresPerOrganism_2D$n <= 9
    ),
    "between010and100_structures" = sum(
      structuresPerOrganism_2D$n >= 10 &
        structuresPerOrganism_2D$n <= 99
    ),
    "above100_structures" = sum(structuresPerOrganism_2D$n >= 100)
  ) %>%
  data.frame()
colnames(tableStructures_2D)[1] <- "organisms"

organismsPerStructure_3D <- pairsOpenDb_3D %>%
  filter(grepl(pattern = "species", x = organismCleaned_dbTaxoTaxonRanks)) %>%
  distinct(organismCleaned, structureCleanedInchikey) %>%
  group_by(structureCleanedInchikey) %>%
  count()

organismsPerStructure_2D <- pairsOpenDb_2D %>%
  filter(grepl(pattern = "species", x = organismCleaned_dbTaxoTaxonRanks)) %>%
  distinct(organismCleaned, structureCleaned_inchikey2D) %>%
  group_by(structureCleaned_inchikey2D) %>%
  count()

tableOrganisms_3D <-
  c(
    "only001_organism" = sum(organismsPerStructure_3D$n == 1),
    "between001and010_organisms" = sum(
      organismsPerStructure_3D$n >= 1 &
        organismsPerStructure_3D$n <= 9
    ),
    "between010and100_organisms" = sum(
      organismsPerStructure_3D$n >= 10 &
        organismsPerStructure_3D$n <= 99
    ),
    "above100_organisms" = sum(organismsPerStructure_3D$n >= 100)
  ) %>%
  data.frame()
colnames(tableOrganisms_3D)[1] <- "structures"

tableOrganisms_2D <-
  c(
    "only001_organism" = sum(organismsPerStructure_2D$n == 1),
    "between001and010_organisms" = sum(
      organismsPerStructure_2D$n >= 1 &
        organismsPerStructure_2D$n <= 9
    ),
    "between010and100_organisms" = sum(
      organismsPerStructure_2D$n >= 10 &
        organismsPerStructure_2D$n <= 99
    ),
    "above100_organisms" = sum(organismsPerStructure_2D$n >= 100)
  ) %>%
  data.frame()
colnames(tableOrganisms_2D)[1] <- "structures"

if (mode == "FULL" | mode == "full") {
  write.table(
    x = dataset,
    file = "../docs/dataset.csv",
    row.names = FALSE,
    quote = FALSE,
    sep = ",",
    fileEncoding = "UTF-8"
  )

  cat(
    paste(
      nrow(pairsOpenDb_3D),
      "|",
      nrow(pairsOpenDb_2D),
      "(3D|2D)",
      "unique referenced structure-organism pairs. \n",
      "They consist of \n",
      nrow(openDbStructure_3D),
      "|",
      nrow(openDbStructure_2D),
      "(3D|2D)",
      "unique curated structures and \n",
      nrow(openDbOrganism),
      "unique organisms,\n",
      "originating from \n",
      nrow(pairsOpenDb_2D %>% distinct(database)),
      "initial open databases. \n",
      "\n",
      "Among 2D structures, \n",
      tableOrganisms_2D[1, 1],
      "are present in only 1 organism, \n",
      tableOrganisms_2D[2, 1],
      "are present in between 1 and 10 organisms, \n",
      tableOrganisms_2D[3, 1],
      "are present in between 10 and 100 organisms, \n",
      tableOrganisms_2D[4, 1],
      "are present in more than 100 organisms. \n",
      "\n",
      "Among organisms, \n",
      tableStructures_2D[1, 1],
      "contain only 1 2D structure, \n",
      tableStructures_2D[2, 1],
      "contain between 1 and 10 2D structures, \n",
      tableStructures_2D[3, 1],
      "contain between 10 and 100 2D structures, \n",
      tableStructures_2D[4, 1],
      "contain more than 100 2D structures. \n",
      sep = " "
    ),
    file = "../docs/metrics.adoc"
  )
}

fwrite(
  x = structuresPerOrganism_3D,
  file = file.path(
    pathDataProcessed,
    "structuresPerOrganism.tsv.gz"
  )
)

fwrite(
  x = organismsPerStructure_3D,
  file = file.path(
    pathDataProcessed,
    "organismsPerStructure.tsv.gz"
  )
)

end <- Sys.time()

cat("Script finished in", format(end - start), "\n")
