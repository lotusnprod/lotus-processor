source("r/log_debug.R")
start <- Sys.time()

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")

log_debug("... libraries")
library(dplyr)
library(splitstackshape)
library(readr)

safety <- FALSE

log_debug("loading the LOTUS, this may take a while")
frozen_lotus <- read_csv(file = file.path(
  pathDataProcessed,
  pathLastFrozen
)) %>%
  mutate(structure_inchikey_2D = substring(
    text = structure_inchikey,
    first = 1,
    last = 14
  ))

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
domain <- frozen_lotus %>%
  filter(!is.na(!!as.name(TAX_LEVEL_1)) |
    !is.na(!!as.name(TAX_LEVEL_2))) %>%
  cSplit(
    "structure_taxonomy_npclassifier_01pathway",
    sep = "|",
    direction = "long"
  ) %>%
  cSplit(
    "structure_taxonomy_npclassifier_02superclass",
    sep = "|",
    direction = "long"
  ) %>%
  cSplit("structure_taxonomy_npclassifier_03class",
    sep = "|",
    direction = "long"
  ) %>%
  filter(!is.na(structure_inchikey_2D) & !is.na(organism_name)) %>%
  mutate(Group = paste(
    organism_taxonomy_01domain,
    organism_taxonomy_02kingdom,
    sep = "_"
  )) %>%
  filter(
    Group == "Eukaryota_Archaeplastida" |
      Group == "Eukaryota_Fungi" |
      Group == "Eukaryota_Metazoa" |
      Group == "Bacteria_NA"
  ) %>%
  mutate(
    Group = if_else(
      condition = Group == "Eukaryota_Archaeplastida",
      true = "Plantae",
      false = if_else(
        condition = Group == "Eukaryota_Fungi",
        true = "Fungi",
        false = if_else(
          condition = Group == "Eukaryota_Metazoa",
          true = "Animalia",
          false = "Bacteria"
        )
      )
    )
  )

domain_unique <- domain %>%
  distinct(structure_inchikey, organism_name, .keep_all = TRUE)

domain_unique_organisms <- domain_unique %>%
  distinct(organism_name, Group) %>%
  group_by(Group) %>%
  count(
    name = "Organisms",
    sort = TRUE
  )

domain_unique_pairs <- domain_unique %>%
  group_by(Group) %>%
  count(
    name = "Structure-Organism Pairs",
    sort = TRUE
  )

domain_unique_structures_3D <- domain_unique %>%
  distinct(Group, structure_inchikey) %>%
  group_by(Group) %>%
  count(
    name = "Chemical Structures",
    sort = TRUE
  )

domain_unique_classes <- domain_unique %>%
  distinct(Group, structure_taxonomy_npclassifier_03class) %>%
  group_by(Group) %>%
  count(
    name = "Chemical Classes",
    sort = TRUE
  )

domain_unique_structures_3D_specific <- domain_unique %>%
  distinct(Group, structure_inchikey) %>%
  group_by(structure_inchikey) %>%
  add_count() %>%
  filter(n == 1) %>%
  ungroup() %>%
  group_by(Group) %>%
  count(
    name = "Specific Chemical Structures",
    sort = TRUE
  )

domain_unique_classes_specific <- domain_unique %>%
  distinct(Group, structure_taxonomy_npclassifier_03class) %>%
  group_by(structure_taxonomy_npclassifier_03class) %>%
  add_count() %>%
  filter(n == 1) %>%
  ungroup() %>%
  group_by(Group) %>%
  count(
    name = "Specific Chemical Classes",
    sort = TRUE
  )

domain <-
  left_join(domain_unique_organisms, domain_unique_pairs) %>%
  left_join(., domain_unique_structures_3D) %>%
  left_join(., domain_unique_structures_3D_specific) %>%
  left_join(., domain_unique_classes) %>%
  left_join(., domain_unique_classes_specific) %>%
  ungroup() %>%
  mutate(
    `Specific Chemical Structures` = paste0(
      `Specific Chemical Structures`,
      " (",
      round(
        x = 100 * `Specific Chemical Structures` / `Chemical Structures`,
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

  write_csv(x = domain, file = "../docs/repartition.csv")
}

end <- Sys.time()

log_debug("Script finished in", format(end - start))
