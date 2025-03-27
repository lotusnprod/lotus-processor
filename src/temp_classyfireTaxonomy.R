source("r/log_debug.R")
log_debug("This script outputs tabular classyfire taxonomy from json")

start <- Sys.time()

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")

log_debug("... libraries")
library(dplyr)
library(jsonlite)
library(readr)
library(tidyr)

# classyfire_old <-
#   read_delim(file = pathDataInterimDictionariesStructureDictionaryClassyfireFile) %>%
#   distinct(structureCleanedInchikey3D, .keep_all = TRUE) %>%
#   select(
#     structureCleanedInchikey = structureCleanedInchikey3D,
#     structureCleaned_classyfire_1kingdom = structureCleaned_1kingdom,
#     structureCleaned_classyfire_2superclass = structureCleaned_2superclass,
#     structureCleaned_classyfire_3class = structureCleaned_3class,
#     structureCleaned_classyfire_4subclass = structureCleaned_4subclass,
#     structureCleaned_classyfire_5directParent = structureCleaned_5directParent
#   )

classyfire_json <-
  jsonlite::fromJSON(
    txt = "../data/external/taxonomySource/structure/classyfire/tax_nodes.json"
  )

classyfire_direct_parent <-
  readr::read_delim(
    file = "../data/interim/dictionaries_full/structure/classyfire/direct_parent.tsv.gz",
    delim = "\t",
    locale = locales
  ) |>
  dplyr::distinct(inchikey, chemontId = directParent)

# classyfire_alternative_parent <-
#   readr::read_delim(file = "../data/interim/dictionaries/structure/classyfire/alternative_parents.tsv.gz",
#              delim = "\t") |>
#   dplyr::distinct(inchikey,
#            chemontId)
#
# classyfire <- dplyr::bind_rows(classyfire_direct_parent,
#                         classyfire_alternative_parent) |>
#   dplyr::distinct()

parent <- classyfire_json |>
  dplyr::filter(is.na(parent_chemont_id))

children <- classyfire_json |>
  dplyr::filter(!is.na(parent_chemont_id))

parents <- classyfire_json |>
  dplyr::distinct(parent_chemont_id)

classyfire_taxonomy_long <-
  dplyr::left_join(
    parent,
    children,
    by = c("chemont_id" = "parent_chemont_id"),
    suffix = c("", "_01kingdom")
  ) |>
  dplyr::left_join(
    children,
    by = c("chemont_id_01kingdom" = "parent_chemont_id"),
    suffix = c("", "_02superclass")
  ) |>
  dplyr::left_join(
    children,
    by = c("chemont_id_02superclass" = "parent_chemont_id"),
    suffix = c("", "_03class")
  ) |>
  dplyr::left_join(
    children,
    by = c("chemont_id_03class" = "parent_chemont_id"),
    suffix = c("", "_04subclass")
  ) |>
  dplyr::left_join(
    children,
    by = c("chemont_id_04subclass" = "parent_chemont_id"),
    suffix = c("", "_05parent")
  ) |>
  dplyr::left_join(
    children,
    by = c("chemont_id_05parent" = "parent_chemont_id"),
    suffix = c("", "_06parent")
  ) |>
  dplyr::left_join(
    children,
    by = c("chemont_id_06parent" = "parent_chemont_id"),
    suffix = c("", "_07parent")
  ) |>
  dplyr::left_join(
    children,
    by = c("chemont_id_07parent" = "parent_chemont_id"),
    suffix = c("", "_08parent")
  ) |>
  dplyr::left_join(
    children,
    by = c("chemont_id_08parent" = "parent_chemont_id"),
    suffix = c("", "_09parent")
  ) |>
  dplyr::left_join(
    children,
    by = c("chemont_id_09parent" = "parent_chemont_id"),
    suffix = c("", "_10parent")
  ) |>
  dplyr::left_join(
    children,
    by = c("chemont_id_10parent" = "parent_chemont_id"),
    suffix = c("", "_11parent")
  ) |>
  dplyr::select(-c(1:3)) |>
  dplyr::select(
    dplyr::starts_with("name"),
    dplyr::starts_with("chemont_id")
  )

kingdom <- dplyr::left_join(
  parents,
  classyfire_taxonomy_long,
  by = c("parent_chemont_id" = "chemont_id_01kingdom")
) |>
  dplyr::filter(!is.na(name_01kingdom)) |>
  dplyr::distinct(chemont_id = parent_chemont_id, name_01kingdom)

superclass <- dplyr::left_join(
  parents,
  classyfire_taxonomy_long,
  by = c("parent_chemont_id" = "chemont_id_02superclass")
) |>
  dplyr::filter(!is.na(name_02superclass)) |>
  dplyr::distinct(
    chemont_id = parent_chemont_id,
    name_01kingdom,
    name_02superclass
  )

class <- dplyr::left_join(
  parents,
  classyfire_taxonomy_long,
  by = c("parent_chemont_id" = "chemont_id_03class")
) |>
  dplyr::filter(!is.na(name_03class)) |>
  dplyr::distinct(
    chemont_id = parent_chemont_id,
    name_01kingdom,
    name_02superclass,
    name_03class
  )

subclass <- dplyr::left_join(
  parents,
  classyfire_taxonomy_long,
  by = c("parent_chemont_id" = "chemont_id_04subclass")
) |>
  dplyr::filter(!is.na(name_04subclass)) |>
  dplyr::distinct(
    chemont_id = parent_chemont_id,
    name_01kingdom,
    name_02superclass,
    name_03class,
    name_04subclass
  )

parent_05 <- dplyr::left_join(
  parents,
  classyfire_taxonomy_long,
  by = c("parent_chemont_id" = "chemont_id_05parent")
) |>
  dplyr::filter(!is.na(name_05parent)) |>
  dplyr::distinct(
    chemont_id = parent_chemont_id,
    name_01kingdom,
    name_02superclass,
    name_03class,
    name_04subclass,
    name_05parent
  )

parent_06 <- dplyr::left_join(
  parents,
  classyfire_taxonomy_long,
  by = c("parent_chemont_id" = "chemont_id_06parent")
) |>
  dplyr::filter(!is.na(name_06parent)) |>
  dplyr::distinct(
    chemont_id = parent_chemont_id,
    name_01kingdom,
    name_02superclass,
    name_03class,
    name_04subclass,
    name_05parent,
    name_06parent
  )

parent_07 <- dplyr::left_join(
  parents,
  classyfire_taxonomy_long,
  by = c("parent_chemont_id" = "chemont_id_07parent")
) |>
  dplyr::filter(!is.na(name_07parent)) |>
  dplyr::distinct(
    chemont_id = parent_chemont_id,
    name_01kingdom,
    name_02superclass,
    name_03class,
    name_04subclass,
    name_05parent,
    name_06parent,
    name_07parent
  )

parent_08 <- dplyr::left_join(
  parents,
  classyfire_taxonomy_long,
  by = c("parent_chemont_id" = "chemont_id_08parent")
) |>
  dplyr::filter(!is.na(name_08parent)) |>
  dplyr::distinct(
    chemont_id = parent_chemont_id,
    name_01kingdom,
    name_02superclass,
    name_03class,
    name_04subclass,
    name_05parent,
    name_06parent,
    name_07parent,
    name_08parent
  )

parent_09 <- dplyr::left_join(
  parents,
  classyfire_taxonomy_long,
  by = c("parent_chemont_id" = "chemont_id_09parent")
) |>
  dplyr::filter(!is.na(name_09parent)) |>
  dplyr::distinct(
    chemont_id = parent_chemont_id,
    name_01kingdom,
    name_02superclass,
    name_03class,
    name_04subclass,
    name_05parent,
    name_06parent,
    name_07parent,
    name_08parent,
    name_09parent
  )

parent_10 <- dplyr::left_join(
  parents,
  classyfire_taxonomy_long,
  by = c("parent_chemont_id" = "chemont_id_10parent")
) |>
  dplyr::filter(!is.na(name_10parent)) |>
  dplyr::distinct(
    chemont_id = parent_chemont_id,
    name_01kingdom,
    name_02superclass,
    name_03class,
    name_04subclass,
    name_05parent,
    name_06parent,
    name_07parent,
    name_08parent,
    name_09parent,
    name_10parent
  )

classyfire_full <- dplyr::bind_rows(
  classyfire_taxonomy_long,
  kingdom,
  superclass,
  class,
  subclass,
  parent_05,
  parent_06,
  parent_07,
  parent_08,
  parent_09,
  parent_10
) %>%
  mutate(
    chemont_id = apply(
      .[, grepl(
        pattern = "chemont_id",
        x = colnames(.),
        fixed = TRUE
      )],
      1,
      function(x) {
        tail(na.omit(x), 1)
      }
    )
  ) %>%
  dplyr::select(
    chemont_id,
    starts_with("name")
  ) %>%
  tidyr::pivot_longer(2:ncol(.), names_prefix = "name_") %>%
  dplyr::filter(!is.na(value))

classyfire_taxonomy_wide <- classyfire_full |>
  tidyr::pivot_wider()

classy_temp <-
  dplyr::left_join(
    classyfire_direct_parent,
    classyfire_taxonomy_wide,
    by = c("chemontId" = "chemont_id")
  ) |>
  dplyr::mutate(
    chemont_id = gsub(
      pattern = "CHEMONTID:",
      replacement = "",
      x = chemontId,
      fixed = TRUE
    )
  ) |>
  dplyr::select(
    inchikey,
    chemont_id,
    dplyr::everything(),
    -chemontId
  )

chemical_taxonomy_2 <<- classy_temp %>%
  dplyr::mutate(
    direct_parent = apply(.[, 3:13], 1, function(x) {
      tail(na.omit(x), 1)
    })
  ) %>%
  dplyr::distinct(
    structure_inchikey = inchikey,
    structure_taxonomy_classyfire_chemontid = chemont_id,
    structure_taxonomy_classyfire_01kingdom = `01kingdom`,
    structure_taxonomy_classyfire_02superclass = `02superclass`,
    structure_taxonomy_classyfire_03class = `03class`,
    structure_taxonomy_classyfire_04directparent = direct_parent,
  )

end <- Sys.time()

log_debug("Script finished in", format(end - start))
