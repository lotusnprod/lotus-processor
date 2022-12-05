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
  readr::read_delim(
    file = pathDataInterimTablesCuratedTable,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  ) |>
  dplyr::filter(database %ni% forbidden_export) |>
  data.frame()

frozen_metadata <-
  readr::read_csv(file = file.path(
    pathDataProcessed,
    pathLastFrozen
  ))

structureMetadata <-
  readr::read_delim(
    file = pathDataInterimDictionariesStructureMetadata,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  ) |>
  dplyr::filter(structureCleaned_stereocenters_unspecified == 0)

referenceOrganismDictionary <-
  readr::read_delim(
    file = pathDataInterimDictionariesReferenceOrganismDictionary,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  ) |>
  dplyr::filter(referenceCleaned_score_titleOrganism == 1)

## structures
inchi <- inhouseDbMinimal |>
  dplyr::filter(structureType == "inchi") |>
  dplyr::distinct(structureValue, .keep_all = TRUE)

smiles <- inhouseDbMinimal |>
  dplyr::filter(structureType == "smiles") |>
  dplyr::distinct(structureValue, .keep_all = TRUE)

nominal <- inhouseDbMinimal |>
  dplyr::filter(structureType == "nominal") |>
  dplyr::distinct(structureValue, .keep_all = TRUE)

inchi_cleaned <- inchi |>
  dplyr::distinct(structureCleanedInchikey, .keep_all = TRUE)

smiles_cleaned <- smiles |>
  dplyr::distinct(structureCleanedInchikey, .keep_all = TRUE)

nominal_cleaned <- nominal |>
  dplyr::distinct(structureCleanedInchikey, .keep_all = TRUE)

inchi_cleaned_3D <- left_join(inchi_cleaned, structureMetadata) |>
  dplyr::filter(structureCleaned_stereocenters_unspecified == 0)

smiles_cleaned_3D <-
  dplyr::left_join(smiles_cleaned, structureMetadata) |>
  dplyr::filter(structureCleaned_stereocenters_unspecified == 0)

nominal_cleaned_3D <-
  dplyr::left_join(nominal_cleaned, structureMetadata) |>
  dplyr::filter(structureCleaned_stereocenters_unspecified == 0)

inchi_frozen <-
  dplyr::semi_join(
    inchi_cleaned,
    frozen_metadata |>
      dplyr::distinct(structure_inchikey),
    by = c("structureCleanedInchikey" = "structure_inchikey")
  ) |>
  dplyr::distinct(structureCleanedInchikey)

smiles_frozen <-
  dplyr::semi_join(
    smiles_cleaned,
    frozen_metadata |>
      dplyr::distinct(structure_inchikey),
    by = c("structureCleanedInchikey" = "structure_inchikey")
  ) |>
  dplyr::distinct(structureCleanedInchikey)

nominal_frozen <-
  dplyr::semi_join(
    nominal_cleaned,
    frozen_metadata |>
      dplyr::distinct(structure_inchikey),
    by = c("structureCleanedInchikey" = "structure_inchikey")
  ) |>
  dplyr::distinct(structureCleanedInchikey)

structures_3D <- frozen_metadata |>
  dplyr::filter(structure_stereocenters_unspecified == 0) |>
  dplyr::distinct(structure_inchikey)

structures_2D <- frozen_metadata |>
  dplyr::filter(structure_stereocenters_unspecified != 0) |>
  dplyr::distinct(structure_inchikey)

## organisms
organism_clean <- inhouseDbMinimal |>
  dplyr::filter(organismType == "clean") |>
  dplyr::distinct(organismValue, .keep_all = TRUE)

organism_dirty <- inhouseDbMinimal |>
  dplyr::filter(organismType == "dirty") |>
  dplyr::distinct(organismValue, .keep_all = TRUE)

organism_clean_cleaned <- organism_clean |>
  dplyr::distinct(organismCleaned, .keep_all = TRUE)

organism_dirty_cleaned <- organism_dirty |>
  dplyr::distinct(organismCleaned, .keep_all = TRUE)

organism_clean_frozen <-
  dplyr::semi_join(
    organism_clean_cleaned,
    frozen_metadata |>
      dplyr::distinct(organism_name),
    by = c("organismCleaned" = "organism_name")
  ) |>
  dplyr::distinct(organismCleaned)

organism_dirty_frozen <-
  dplyr::semi_join(
    organism_dirty_cleaned,
    frozen_metadata |>
      dplyr::distinct(organism_name),
    by = c("organismCleaned" = "organism_name")
  ) |>
  dplyr::distinct(organismCleaned)

## references
reference_original <- inhouseDbMinimal |>
  dplyr::filter(referenceType == "original") |>
  dplyr::distinct(organismValue, .keep_all = TRUE)

reference_pubmed <- inhouseDbMinimal |>
  dplyr::filter(referenceType == "pubmed") |>
  dplyr::distinct(organismValue, .keep_all = TRUE)

reference_doi <- inhouseDbMinimal |>
  dplyr::filter(referenceType == "doi") |>
  dplyr::distinct(organismValue, .keep_all = TRUE)

reference_title <- inhouseDbMinimal |>
  dplyr::filter(referenceType == "title") |>
  dplyr::distinct(organismValue, .keep_all = TRUE)

reference_split <- inhouseDbMinimal |>
  dplyr::filter(referenceType == "split") |>
  dplyr::distinct(organismValue, .keep_all = TRUE)

reference_publishingDetails <- inhouseDbMinimal |>
  dplyr::filter(referenceType == "publishingDetails") |>
  dplyr::distinct(organismValue, .keep_all = TRUE)

reference_original_cleaned <- dplyr::left_join(
  reference_original,
  referenceOrganismDictionary
) |>
  dplyr::distinct(referenceCleanedTitle, .keep_all = TRUE)

reference_pubmed_cleaned <- dplyr::left_join(
  reference_pubmed,
  referenceOrganismDictionary
) |>
  dplyr::distinct(referenceCleanedTitle, .keep_all = TRUE)

reference_doi_cleaned <- dplyr::left_join(
  reference_doi,
  referenceOrganismDictionary
) |>
  dplyr::distinct(referenceCleanedTitle, .keep_all = TRUE)

reference_title_cleaned <- dplyr::left_join(
  reference_title,
  referenceOrganismDictionary
) |>
  dplyr::distinct(referenceCleanedTitle, .keep_all = TRUE)

reference_split_cleaned <- dplyr::left_join(
  reference_split,
  referenceOrganismDictionary
) |>
  dplyr::distinct(referenceCleanedTitle, .keep_all = TRUE)

reference_publishingDetails_cleaned <-
  left_join(
    reference_publishingDetails,
    referenceOrganismDictionary
  ) |>
  dplyr::distinct(referenceCleanedTitle, .keep_all = TRUE)

reference_original_cleaned_plus <- reference_original_cleaned |>
  dplyr::filter(referenceCleaned_score_titleOrganism == 1)

reference_pubmed_cleaned_plus <- reference_pubmed_cleaned |>
  dplyr::filter(referenceCleaned_score_titleOrganism == 1)

reference_doi_cleaned_plus <- reference_doi_cleaned |>
  dplyr::filter(referenceCleaned_score_titleOrganism == 1)

reference_title_cleaned_plus <- reference_title_cleaned |>
  dplyr::filter(referenceCleaned_score_titleOrganism == 1)

reference_split_cleaned_plus <- reference_split_cleaned |>
  dplyr::filter(referenceCleaned_score_titleOrganism == 1)

reference_publishingDetails_cleaned_plus <-
  reference_publishingDetails_cleaned |>
  dplyr::filter(referenceCleaned_score_titleOrganism == 1)

reference_original_frozen <-
  dplyr::semi_join(
    reference_original_cleaned_plus,
    frozen_metadata |>
      dplyr::distinct(reference_doi),
    by = c("referenceCleanedDoi" = "reference_doi")
  ) |>
  dplyr::distinct(referenceCleanedDoi)

reference_pubmed_frozen <-
  dplyr::semi_join(
    reference_pubmed_cleaned_plus,
    frozen_metadata |>
      dplyr::distinct(reference_doi),
    by = c("referenceCleanedDoi" = "reference_doi")
  ) |>
  dplyr::distinct(referenceCleanedDoi)

reference_doi_frozen <-
  dplyr::semi_join(
    reference_doi_cleaned_plus,
    frozen_metadata |>
      dplyr::distinct(reference_doi),
    by = c("referenceCleanedDoi" = "reference_doi")
  ) |>
  dplyr::distinct(referenceCleanedDoi)

reference_title_frozen <-
  dplyr::semi_join(
    reference_title_cleaned_plus,
    frozen_metadata |>
      dplyr::distinct(reference_doi),
    by = c("referenceCleanedDoi" = "reference_doi")
  ) |>
  dplyr::distinct(referenceCleanedDoi)

reference_split_frozen <-
  dplyr::semi_join(
    reference_split_cleaned_plus,
    frozen_metadata |>
      dplyr::distinct(reference_doi),
    by = c("referenceCleanedDoi" = "reference_doi")
  ) |>
  dplyr::distinct(referenceCleanedDoi)

reference_publishingDetails_frozen <-
  dplyr::semi_join(
    reference_publishingDetails_cleaned_plus,
    frozen_metadata |>
      dplyr::distinct(reference_doi),
    by = c("referenceCleanedDoi" = "reference_doi")
  ) |>
  dplyr::distinct(referenceCleanedDoi)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
