cat("This script outputs tabular classyfire taxonomy from json \n")

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")

cat("... libraries \n")
library(data.table)
library(jsonlite)
library(tidyverse)
source("r/vroom_safe.R")

classyfire_old <-
  vroom_read_safe(path = pathDataInterimDictionariesStructureDictionaryClassyfireFile) %>%
  distinct(structureCleanedInchikey3D, .keep_all = TRUE) %>%
  select(
    structureCleanedInchikey = structureCleanedInchikey3D,
    structureCleaned_classyfire_1kingdom = structureCleaned_1kingdom,
    structureCleaned_classyfire_2superclass = structureCleaned_2superclass,
    structureCleaned_classyfire_3class = structureCleaned_3class,
    structureCleaned_classyfire_4subclass = structureCleaned_4subclass,
    structureCleaned_classyfire_5directParent = structureCleaned_5directParent
  ) %>%
  tibble()

classyfire_json <-
  fromJSON(txt = "../data/external/taxonomySource/structure/classyfire/tax_nodes.json")

classyfire_new <-
  vroom_read_safe(path = "../data/interim/dictionaries/structure/classyfire/direct_parent.tsv.gz")

parent <- classyfire_json %>% filter(is.na(parent_chemont_id))

children <- classyfire_json %>% filter(!is.na(parent_chemont_id))

parents <- classyfire_json %>% distinct(parent_chemont_id)

classyfire_taxonomy_long <-
  left_join(
    parent,
    children,
    by = c("chemont_id" = "parent_chemont_id"),
    suffix = c("", "_01kingdom")
  ) %>%
  left_join(
    .,
    children,
    by = c("chemont_id_01kingdom" = "parent_chemont_id"),
    suffix = c("", "_02superclass")
  ) %>%
  left_join(
    .,
    children,
    by = c("chemont_id_02superclass" = "parent_chemont_id"),
    suffix = c("", "_03class")
  ) %>%
  left_join(
    .,
    children,
    by = c("chemont_id_03class" = "parent_chemont_id"),
    suffix = c("", "_04subclass")
  ) %>%
  left_join(
    .,
    children,
    by = c("chemont_id_04subclass" = "parent_chemont_id"),
    suffix = c("", "_05parent")
  ) %>%
  left_join(
    .,
    children,
    by = c("chemont_id_05parent" = "parent_chemont_id"),
    suffix = c("", "_06parent")
  ) %>%
  left_join(
    .,
    children,
    by = c("chemont_id_06parent" = "parent_chemont_id"),
    suffix = c("", "_07parent")
  ) %>%
  left_join(
    .,
    children,
    by = c("chemont_id_07parent" = "parent_chemont_id"),
    suffix = c("", "_08parent")
  ) %>%
  left_join(
    .,
    children,
    by = c("chemont_id_08parent" = "parent_chemont_id"),
    suffix = c("", "_09parent")
  ) %>%
  left_join(
    .,
    children,
    by = c("chemont_id_09parent" = "parent_chemont_id"),
    suffix = c("", "_10parent")
  ) %>%
  left_join(
    .,
    children,
    by = c("chemont_id_10parent" = "parent_chemont_id"),
    suffix = c("", "_11parent")
  ) %>%
  select(-c(1:3)) %>%
  select(
    starts_with("name"),
    starts_with("chemont_id")
  )

kingdom <- left_join(
  parents,
  classyfire_taxonomy_long,
  by = c("parent_chemont_id" = "chemont_id_01kingdom")
) %>%
  filter(!is.na(name_01kingdom)) %>%
  distinct(chemont_id = parent_chemont_id, name_01kingdom)

superclass <- left_join(
  parents,
  classyfire_taxonomy_long,
  by = c("parent_chemont_id" = "chemont_id_02superclass")
) %>%
  filter(!is.na(name_02superclass)) %>%
  distinct(
    chemont_id = parent_chemont_id,
    name_01kingdom,
    name_02superclass
  )

class <- left_join(parents,
  classyfire_taxonomy_long,
  by = c("parent_chemont_id" = "chemont_id_03class")
) %>%
  filter(!is.na(name_03class)) %>%
  distinct(
    chemont_id = parent_chemont_id,
    name_01kingdom,
    name_02superclass,
    name_03class
  )

subclass <- left_join(
  parents,
  classyfire_taxonomy_long,
  by = c("parent_chemont_id" = "chemont_id_04subclass")
) %>%
  filter(!is.na(name_04subclass)) %>%
  distinct(
    chemont_id = parent_chemont_id,
    name_01kingdom,
    name_02superclass,
    name_03class,
    name_04subclass
  )

parent_05 <- left_join(
  parents,
  classyfire_taxonomy_long,
  by = c("parent_chemont_id" = "chemont_id_05parent")
) %>%
  filter(!is.na(name_05parent)) %>%
  distinct(
    chemont_id = parent_chemont_id,
    name_01kingdom,
    name_02superclass,
    name_03class,
    name_04subclass,
    name_05parent
  )

parent_06 <- left_join(
  parents,
  classyfire_taxonomy_long,
  by = c("parent_chemont_id" = "chemont_id_06parent")
) %>%
  filter(!is.na(name_06parent)) %>%
  distinct(
    chemont_id = parent_chemont_id,
    name_01kingdom,
    name_02superclass,
    name_03class,
    name_04subclass,
    name_05parent,
    name_06parent
  )

parent_07 <- left_join(
  parents,
  classyfire_taxonomy_long,
  by = c("parent_chemont_id" = "chemont_id_07parent")
) %>%
  filter(!is.na(name_07parent)) %>%
  distinct(
    chemont_id = parent_chemont_id,
    name_01kingdom,
    name_02superclass,
    name_03class,
    name_04subclass,
    name_05parent,
    name_06parent,
    name_07parent
  )

parent_08 <- left_join(
  parents,
  classyfire_taxonomy_long,
  by = c("parent_chemont_id" = "chemont_id_08parent")
) %>%
  filter(!is.na(name_08parent)) %>%
  distinct(
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

parent_09 <- left_join(
  parents,
  classyfire_taxonomy_long,
  by = c("parent_chemont_id" = "chemont_id_09parent")
) %>%
  filter(!is.na(name_09parent)) %>%
  distinct(
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

parent_10 <- left_join(
  parents,
  classyfire_taxonomy_long,
  by = c("parent_chemont_id" = "chemont_id_10parent")
) %>%
  filter(!is.na(name_10parent)) %>%
  distinct(
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

classyfire_full <- bind_rows(
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
  mutate(chemont_id = apply(.[, grepl(
    pattern = "chemont_id",
    x = colnames(.),
    fixed = TRUE
  )], 1, function(x) {
    tail(na.omit(x), 1)
  })) %>%
  select(
    chemont_id,
    starts_with("name")
  ) %>%
  pivot_longer(2:ncol(.), names_prefix = "name_") %>%
  filter(!is.na(value))

classyfire_taxonomy_wide <- classyfire_full %>%
  pivot_wider()

classy_temp <-
  left_join(classyfire_new,
    classyfire_taxonomy_wide,
    by = c("directParent" = "chemont_id")
  ) %>%
  mutate(chemont_id = gsub(
    pattern = "CHEMONTID:",
    replacement = "",
    x = directParent,
    fixed = TRUE
  )) %>%
  select(
    inchikey,
    chemont_id, everything(), -directParent
  )

end <- Sys.time()

cat("Script finished in", format(end - start), "\n")
