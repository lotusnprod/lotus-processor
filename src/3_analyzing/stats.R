source("r/log_debug.R")
start <- Sys.time()

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")

log_debug("... libraries")
library(dplyr)
library(readr)

"%ni%" <- Negate("%in%")

log_debug("loading db, if running fullmode, this may take a while")
inhouseDbMinimal <-
  read_delim(
    file = pathDataInterimTablesCuratedTable,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  ) %>%
  filter(database %ni% forbidden_export) %>%
  data.frame()

frozen_metadata <-
  read_csv(file = file.path(
    pathDataProcessed,
    pathLastFrozen
  ))

structureMetadata <-
  read_delim(
    file = pathDataInterimDictionariesStructureMetadata,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  ) %>%
  filter(structureCleaned_stereocenters_unspecified == 0)

referenceOrganismDictionary <-
  read_delim(
    file = pathDataInterimDictionariesReferenceOrganismDictionary,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  ) %>%
  filter(referenceCleaned_score_titleOrganism == 1)

## structures
inchi <- inhouseDbMinimal %>%
  filter(structureType == "inchi") %>%
  distinct(structureValue, .keep_all = TRUE)

smiles <- inhouseDbMinimal %>%
  filter(structureType == "smiles") %>%
  distinct(structureValue, .keep_all = TRUE)

nominal <- inhouseDbMinimal %>%
  filter(structureType == "nominal") %>%
  distinct(structureValue, .keep_all = TRUE)

inchi_cleaned <- inchi %>%
  distinct(structureCleanedInchikey, .keep_all = TRUE)

smiles_cleaned <- smiles %>%
  distinct(structureCleanedInchikey, .keep_all = TRUE)

nominal_cleaned <- nominal %>%
  distinct(structureCleanedInchikey, .keep_all = TRUE)

inchi_cleaned_3D <- left_join(inchi_cleaned, structureMetadata) %>%
  filter(structureCleaned_stereocenters_unspecified == 0)

smiles_cleaned_3D <-
  left_join(smiles_cleaned, structureMetadata) %>%
  filter(structureCleaned_stereocenters_unspecified == 0)

nominal_cleaned_3D <-
  left_join(nominal_cleaned, structureMetadata) %>%
  filter(structureCleaned_stereocenters_unspecified == 0)

inchi_frozen <-
  semi_join(
    inchi_cleaned,
    frozen_metadata %>% distinct(structure_inchikey),
    by = c("structureCleanedInchikey" = "structure_inchikey")
  ) %>%
  distinct(structureCleanedInchikey)

smiles_frozen <-
  semi_join(
    smiles_cleaned,
    frozen_metadata %>% distinct(structure_inchikey),
    by = c("structureCleanedInchikey" = "structure_inchikey")
  ) %>%
  distinct(structureCleanedInchikey)

nominal_frozen <-
  semi_join(
    nominal_cleaned,
    frozen_metadata %>% distinct(structure_inchikey),
    by = c("structureCleanedInchikey" = "structure_inchikey")
  ) %>%
  distinct(structureCleanedInchikey)

structures_3D <- frozen_metadata %>%
  filter(structure_stereocenters_unspecified == 0) %>%
  distinct(structure_inchikey)

structures_2D <- frozen_metadata %>%
  filter(structure_stereocenters_unspecified != 0) %>%
  distinct(structure_inchikey)

## organisms
organism_clean <- inhouseDbMinimal %>%
  filter(organismType == "clean") %>%
  distinct(organismValue, .keep_all = TRUE)

organism_dirty <- inhouseDbMinimal %>%
  filter(organismType == "dirty") %>%
  distinct(organismValue, .keep_all = TRUE)

organism_clean_cleaned <- organism_clean %>%
  distinct(organismCleaned, .keep_all = TRUE)

organism_dirty_cleaned <- organism_dirty %>%
  distinct(organismCleaned, .keep_all = TRUE)

organism_clean_frozen <-
  semi_join(
    organism_clean_cleaned,
    frozen_metadata %>% distinct(organism_name),
    by = c("organismCleaned" = "organism_name")
  ) %>%
  distinct(organismCleaned)

organism_dirty_frozen <-
  semi_join(
    organism_dirty_cleaned,
    frozen_metadata %>% distinct(organism_name),
    by = c("organismCleaned" = "organism_name")
  ) %>%
  distinct(organismCleaned)

## references
reference_original <- inhouseDbMinimal %>%
  filter(referenceType == "original") %>%
  distinct(organismValue, .keep_all = TRUE)

reference_pubmed <- inhouseDbMinimal %>%
  filter(referenceType == "pubmed") %>%
  distinct(organismValue, .keep_all = TRUE)

reference_doi <- inhouseDbMinimal %>%
  filter(referenceType == "doi") %>%
  distinct(organismValue, .keep_all = TRUE)

reference_title <- inhouseDbMinimal %>%
  filter(referenceType == "title") %>%
  distinct(organismValue, .keep_all = TRUE)

reference_split <- inhouseDbMinimal %>%
  filter(referenceType == "split") %>%
  distinct(organismValue, .keep_all = TRUE)

reference_publishingDetails <- inhouseDbMinimal %>%
  filter(referenceType == "publishingDetails") %>%
  distinct(organismValue, .keep_all = TRUE)

reference_original_cleaned <- left_join(
  reference_original,
  referenceOrganismDictionary
) %>%
  distinct(referenceCleanedTitle, .keep_all = TRUE)

reference_pubmed_cleaned <- left_join(
  reference_pubmed,
  referenceOrganismDictionary
) %>%
  distinct(referenceCleanedTitle, .keep_all = TRUE)

reference_doi_cleaned <- left_join(
  reference_doi,
  referenceOrganismDictionary
) %>%
  distinct(referenceCleanedTitle, .keep_all = TRUE)

reference_title_cleaned <- left_join(
  reference_title,
  referenceOrganismDictionary
) %>%
  distinct(referenceCleanedTitle, .keep_all = TRUE)

reference_split_cleaned <- left_join(
  reference_split,
  referenceOrganismDictionary
) %>%
  distinct(referenceCleanedTitle, .keep_all = TRUE)

reference_publishingDetails_cleaned <-
  left_join(
    reference_publishingDetails,
    referenceOrganismDictionary
  ) %>%
  distinct(referenceCleanedTitle, .keep_all = TRUE)

reference_original_cleaned_plus <- reference_original_cleaned %>%
  filter(referenceCleaned_score_titleOrganism == 1)

reference_pubmed_cleaned_plus <- reference_pubmed_cleaned %>%
  filter(referenceCleaned_score_titleOrganism == 1)

reference_doi_cleaned_plus <- reference_doi_cleaned %>%
  filter(referenceCleaned_score_titleOrganism == 1)

reference_title_cleaned_plus <- reference_title_cleaned %>%
  filter(referenceCleaned_score_titleOrganism == 1)

reference_split_cleaned_plus <- reference_split_cleaned %>%
  filter(referenceCleaned_score_titleOrganism == 1)

reference_publishingDetails_cleaned_plus <-
  reference_publishingDetails_cleaned %>%
  filter(referenceCleaned_score_titleOrganism == 1)

reference_original_frozen <-
  semi_join(
    reference_original_cleaned_plus,
    frozen_metadata %>% distinct(reference_doi),
    by = c("referenceCleanedDoi" = "reference_doi")
  ) %>%
  distinct(referenceCleanedDoi)

reference_pubmed_frozen <-
  semi_join(
    reference_pubmed_cleaned_plus,
    frozen_metadata %>% distinct(reference_doi),
    by = c("referenceCleanedDoi" = "reference_doi")
  ) %>%
  distinct(referenceCleanedDoi)

reference_doi_frozen <-
  semi_join(
    reference_doi_cleaned_plus,
    frozen_metadata %>% distinct(reference_doi),
    by = c("referenceCleanedDoi" = "reference_doi")
  ) %>%
  distinct(referenceCleanedDoi)

reference_title_frozen <-
  semi_join(
    reference_title_cleaned_plus,
    frozen_metadata %>% distinct(reference_doi),
    by = c("referenceCleanedDoi" = "reference_doi")
  ) %>%
  distinct(referenceCleanedDoi)

reference_split_frozen <-
  semi_join(
    reference_split_cleaned_plus,
    frozen_metadata %>% distinct(reference_doi),
    by = c("referenceCleanedDoi" = "reference_doi")
  ) %>%
  distinct(referenceCleanedDoi)

reference_publishingDetails_frozen <-
  semi_join(
    reference_publishingDetails_cleaned_plus,
    frozen_metadata %>% distinct(reference_doi),
    by = c("referenceCleanedDoi" = "reference_doi")
  ) %>%
  distinct(referenceCleanedDoi)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
