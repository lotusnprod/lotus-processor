source("r/log_debug.R")
start <- Sys.time()

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")

log_debug("... libraries")
library(dplyr)
library(splitstackshape)
library(readr)

log_debug("loading the LOTUS, this may take a while")
frozen_lotus <- readr::read_csv(
  file = file.path(
    pathDataProcessed,
    pathLastFrozen
  )
) |>
  dplyr::mutate(
    structure_inchikey_2D = substring(
      text = structure_inchikey,
      first = 1,
      last = 14
    )
  )

TAX_LEVEL_1 <- "organism_taxonomy_01domain"
TAX_LEVEL_2 <- "organism_taxonomy_02kingdom"
# TAX_LEVEL <- "organism_taxonomy_03phylum"
# TAX_LEVEL <- "organism_taxonomy_04class"
# TAX_LEVEL <- "organism_taxonomy_05order"
# TAX_LEVEL <- "organism_taxonomy_06family"
# TAX_LEVEL <- "organism_taxonomy_07tribe"
# TAX_LEVEL <- "organism_taxonomy_08genus"
# TAX_LEVEL <- "organism_taxonomy_09species"
# TAX_LEVEL <- "organism_taxonomy_10varietas"

## BEWARE some taxa are not in OTL!
domain <- frozen_lotus |>
  dplyr::filter(
    !is.na(!!as.name(TAX_LEVEL_1)) |
      !is.na(!!as.name(TAX_LEVEL_2))
  ) |>
  splitstackshape::cSplit(
    "structure_taxonomy_npclassifier_01pathway",
    sep = "|",
    direction = "long"
  ) |>
  splitstackshape::cSplit(
    "structure_taxonomy_npclassifier_02superclass",
    sep = "|",
    direction = "long"
  ) |>
  splitstackshape::cSplit(
    "structure_taxonomy_npclassifier_03class",
    sep = "|",
    direction = "long"
  ) |>
  dplyr::filter(
    !is.na(structure_inchikey_2D) &
      !is.na(organism_name)
  ) |>
  dplyr::mutate(
    Group = paste(
      organism_taxonomy_01domain,
      organism_taxonomy_02kingdom,
      sep = "_"
    )
  ) |>
  dplyr::filter(
    Group == "Eukaryota_Archaeplastida" |
      Group == "Eukaryota_Fungi" |
      Group == "Eukaryota_Metazoa" |
      Group == "Bacteria_NA"
  ) |>
  dplyr::mutate(
    Group = dplyr::if_else(
      condition = Group == "Eukaryota_Archaeplastida",
      true = "Plantae",
      false = dplyr::if_else(
        condition = Group == "Eukaryota_Fungi",
        true = "Fungi",
        false = dplyr::if_else(
          condition = Group == "Eukaryota_Metazoa",
          true = "Animalia",
          false = "Bacteria"
        )
      )
    )
  )

domain_unique <- domain |>
  dplyr::distinct(structure_inchikey_2D, organism_name, .keep_all = TRUE)

domain_unique_organisms <- domain_unique |>
  dplyr::distinct(organism_name, Group) |>
  dplyr::group_by(Group) |>
  dplyr::count(
    name = "Organisms",
    sort = TRUE
  )

domain_unique_pairs <- domain_unique |>
  dplyr::group_by(Group) |>
  dplyr::count(
    name = "2D Structure-Organism Pairs",
    sort = TRUE
  )

domain_unique_structures_2D <- domain_unique |>
  dplyr::distinct(Group, structure_inchikey_2D) |>
  dplyr::group_by(Group) |>
  dplyr::count(
    name = "2D Chemical Structures",
    sort = TRUE
  )

domain_unique_classes <- domain_unique |>
  dplyr::distinct(Group, structure_taxonomy_npclassifier_03class) |>
  dplyr::group_by(Group) |>
  dplyr::count(
    name = "Chemical Classes",
    sort = TRUE
  )

domain_unique_structures_2D_specific <- domain_unique |>
  dplyr::distinct(Group, structure_inchikey_2D) |>
  dplyr::group_by(structure_inchikey_2D) |>
  dplyr::add_count() |>
  dplyr::filter(n == 1) |>
  dplyr::ungroup() |>
  dplyr::group_by(Group) |>
  dplyr::count(
    name = "Specific 2D Chemical Structures",
    sort = TRUE
  )

domain_unique_classes_specific <- domain_unique |>
  dplyr::distinct(Group, structure_taxonomy_npclassifier_03class) |>
  dplyr::group_by(structure_taxonomy_npclassifier_03class) |>
  dplyr::add_count() |>
  dplyr::filter(n == 1) |>
  dplyr::ungroup() |>
  dplyr::group_by(Group) |>
  dplyr::count(
    name = "Specific Chemical Classes",
    sort = TRUE
  )

domain <-
  dplyr::left_join(domain_unique_organisms, domain_unique_pairs) |>
  dplyr::left_join(domain_unique_structures_2D) |>
  dplyr::left_join(domain_unique_structures_2D_specific) |>
  dplyr::left_join(domain_unique_classes) |>
  dplyr::left_join(domain_unique_classes_specific) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    `Specific 2D Chemical Structures` = paste0(
      `Specific 2D Chemical Structures`,
      " (",
      round(
        x = 100 * `Specific 2D Chemical Structures` / `2D Chemical Structures`,
        digits = 0
      ),
      "%)"
    ),
    `Specific Chemical Classes` = paste0(
      `Specific Chemical Classes`,
      " (",
      round(
        x = 100 * `Specific Chemical Classes` / `Chemical Classes`,
        digits = 0
      ),
      "%)"
    )
  )

if (safety == TRUE) {
  log_debug(
    "Exporting to",
    "../docs/repartition.csv"
  )

  readr::write_csv(
    x = domain,
    file = "../docs/repartition.csv",
    na = ""
  )
}

end <- Sys.time()

log_debug("Script finished in", format(end - start))
