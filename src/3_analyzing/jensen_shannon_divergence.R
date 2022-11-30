source("r/log_debug.R")
log_debug(
  "This script calculcates Jensen-Shannon Divergence at given chemical and biological levels"
)

start <- Sys.time()

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")
source("r/get_jsd.R")

log_debug("... libraries")
library(data.table)
library(dplyr)
library(philentropy)
library(readr)
library(tidyr)

log_debug("loading the LOTUS, this may take a while")
table <- read_csv(file = file.path(
  pathDataProcessed,
  pathLastFrozen
)) %>%
  select(
    structure_smiles_2D,
    structure_taxonomy_npclassifier_01pathway,
    structure_taxonomy_npclassifier_02superclass,
    structure_taxonomy_npclassifier_03class,
    organism_name,
    organism_taxonomy_ottid,
    organism_taxonomy_01domain,
    organism_taxonomy_02kingdom,
    organism_taxonomy_03phylum,
    organism_taxonomy_04class,
    organism_taxonomy_05order,
    organism_taxonomy_06family,
    organism_taxonomy_07tribe,
    organism_taxonomy_08genus,
    organism_taxonomy_09species,
    organism_taxonomy_10varietas
  ) %>%
  distinct()

log_debug("counting occurrences")
table_counted <- table %>%
  group_by(
    structure_taxonomy_npclassifier_01pathway,
    organism_taxonomy_02kingdom
  ) %>%
  add_count(name = "structure_taxonomy_npclassifier_01pathway_organism_taxonomy_02kingdom") %>%
  group_by(
    structure_taxonomy_npclassifier_02superclass,
    organism_taxonomy_02kingdom
  ) %>%
  add_count(name = "structure_taxonomy_npclassifier_02superclass_organism_taxonomy_02kingdom") %>%
  group_by(
    structure_taxonomy_npclassifier_03class,
    organism_taxonomy_02kingdom
  ) %>%
  add_count(name = "structure_taxonomy_npclassifier_03class_organism_taxonomy_02kingdom") %>%
  group_by(
    structure_taxonomy_npclassifier_01pathway,
    organism_taxonomy_06family
  ) %>%
  add_count(name = "structure_taxonomy_npclassifier_01pathway_organism_taxonomy_06family") %>%
  group_by(
    structure_taxonomy_npclassifier_02superclass,
    organism_taxonomy_06family
  ) %>%
  add_count(name = "structure_taxonomy_npclassifier_02superclass_organism_taxonomy_06family") %>%
  group_by(
    structure_taxonomy_npclassifier_03class,
    organism_taxonomy_06family
  ) %>%
  add_count(name = "structure_taxonomy_npclassifier_03class_organism_taxonomy_06family") %>%
  group_by(
    structure_taxonomy_npclassifier_01pathway,
    organism_taxonomy_08genus
  ) %>%
  add_count(name = "structure_taxonomy_npclassifier_01pathway_organism_taxonomy_08genus") %>%
  group_by(
    structure_taxonomy_npclassifier_02superclass,
    organism_taxonomy_08genus
  ) %>%
  add_count(name = "structure_taxonomy_npclassifier_02superclass_organism_taxonomy_08genus") %>%
  group_by(
    structure_taxonomy_npclassifier_03class,
    organism_taxonomy_08genus
  ) %>%
  add_count(name = "structure_taxonomy_npclassifier_03class_organism_taxonomy_08genus") %>%
  ungroup() %>%
  distinct(
    structure_taxonomy_npclassifier_01pathway,
    structure_taxonomy_npclassifier_02superclass,
    structure_taxonomy_npclassifier_03class,
    organism_taxonomy_02kingdom,
    organism_taxonomy_06family,
    organism_taxonomy_08genus,
    structure_taxonomy_npclassifier_01pathway_organism_taxonomy_02kingdom,
    structure_taxonomy_npclassifier_02superclass_organism_taxonomy_02kingdom,
    structure_taxonomy_npclassifier_03class_organism_taxonomy_02kingdom,
    structure_taxonomy_npclassifier_01pathway_organism_taxonomy_06family,
    structure_taxonomy_npclassifier_02superclass_organism_taxonomy_06family,
    structure_taxonomy_npclassifier_03class_organism_taxonomy_06family,
    structure_taxonomy_npclassifier_01pathway_organism_taxonomy_08genus,
    structure_taxonomy_npclassifier_02superclass_organism_taxonomy_08genus,
    structure_taxonomy_npclassifier_03class_organism_taxonomy_08genus
  ) %>%
  arrange(
    structure_taxonomy_npclassifier_01pathway,
    structure_taxonomy_npclassifier_02superclass,
    structure_taxonomy_npclassifier_03class
  )

log_debug("computing JSD at the chemical class level ...")
chem_level <- "structure_taxonomy_npclassifier_03class"

log_debug("... at the biological genus level")
bio_level <- "organism_taxonomy_08genus"
Y <-
  unique(table_counted[, chem_level][!is.na(table_counted[, chem_level]) &
    !is.na(table_counted[, bio_level])])
X <- seq_along(Y)

class_1 <- invisible(lapply(
  FUN = get_jsd,
  X = X
)) %>%
  rbindlist()

log_debug("... at the biological family level")
bio_level <- "organism_taxonomy_06family"
Y <-
  unique(table_counted[, chem_level][!is.na(table_counted[, chem_level]) &
    !is.na(table_counted[, bio_level])])
X <- seq_along(Y)

class_2 <- invisible(lapply(
  FUN = get_jsd,
  X = X
)) %>%
  rbindlist()

log_debug("... at the biological kingdom level")
bio_level <- "organism_taxonomy_02kingdom"
Y <-
  unique(table_counted[, chem_level][!is.na(table_counted[, chem_level]) &
    !is.na(table_counted[, bio_level])])
X <- seq_along(Y)

class_3 <- invisible(lapply(
  FUN = get_jsd,
  X = X
)) %>%
  rbindlist()

class <- class_1 %>%
  full_join(class_2) %>%
  full_join(class_3) %>%
  data.frame()

colnames(class)[grepl(pattern = "organism", x = colnames(class))] <-
  paste(colnames(class)[grepl(pattern = "organism", x = colnames(class))], "JSD", sep = "_")

log_debug("... at the chemical superclass level ...")
chem_level <- "structure_taxonomy_npclassifier_02superclass"

log_debug("... at the biological genus level")
bio_level <- "organism_taxonomy_08genus"
Y <-
  unique(table_counted[, chem_level][!is.na(table_counted[, chem_level]) &
    !is.na(table_counted[, bio_level])])
X <- seq_along(Y)

superclass_1 <- invisible(lapply(
  FUN = get_jsd,
  X = X
)) %>%
  rbindlist()

log_debug("... at the biological family level")
bio_level <- "organism_taxonomy_06family"
Y <-
  unique(table_counted[, chem_level][!is.na(table_counted[, chem_level]) &
    !is.na(table_counted[, bio_level])])
X <- seq_along(Y)

superclass_2 <- invisible(lapply(
  FUN = get_jsd,
  X = X
)) %>%
  rbindlist()

log_debug("... at the biological kingdom level")
bio_level <- "organism_taxonomy_02kingdom"
Y <-
  unique(table_counted[, chem_level][!is.na(table_counted[, chem_level]) &
    !is.na(table_counted[, bio_level])])
X <- seq_along(Y)

superclass_3 <- invisible(lapply(
  FUN = get_jsd,
  X = X
)) %>%
  rbindlist()

superclass <- superclass_1 %>%
  full_join(superclass_2) %>%
  full_join(superclass_3) %>%
  data.frame()

colnames(superclass)[grepl(pattern = "organism", x = colnames(superclass))] <-
  paste(colnames(superclass)[grepl(pattern = "organism", x = colnames(superclass))], "JSD", sep = "_")

log_debug("... at the chemical pathway level ...")
chem_level <- "structure_taxonomy_npclassifier_01pathway"

log_debug("... at the biological genus level")
bio_level <- "organism_taxonomy_08genus"
Y <-
  unique(table_counted[, chem_level][!is.na(table_counted[, chem_level]) &
    !is.na(table_counted[, bio_level])])
X <- seq_along(Y)

pathway_1 <- invisible(lapply(
  FUN = get_jsd,
  X = X
)) %>%
  rbindlist()

log_debug("... at the biological family level")
bio_level <- "organism_taxonomy_06family"
Y <-
  unique(table_counted[, chem_level][!is.na(table_counted[, chem_level]) &
    !is.na(table_counted[, bio_level])])
X <- seq_along(Y)

pathway_2 <- invisible(lapply(
  FUN = get_jsd,
  X = X
)) %>%
  rbindlist()

log_debug("... at the biological kingdom level")
bio_level <- "organism_taxonomy_02kingdom"
Y <-
  unique(table_counted[, chem_level][!is.na(table_counted[, chem_level]) &
    !is.na(table_counted[, bio_level])])
X <- seq_along(Y)

pathway_3 <- invisible(lapply(
  FUN = get_jsd,
  X = X
)) %>%
  rbindlist()

pathway <- pathway_1 %>%
  full_join(pathway_2) %>%
  full_join(pathway_3) %>%
  data.frame()

colnames(pathway)[grepl(pattern = "organism", x = colnames(pathway))] <-
  paste(colnames(pathway)[grepl(pattern = "organism", x = colnames(pathway))], "JSD", sep = "_")

log_debug("pivoting results")
class_pivoted <- class %>%
  pivot_longer(cols = 2:4, values_to = c("value_score")) %>%
  mutate(
    name_structure = "structure_taxonomy_npclassifier_03class",
    value_structure = structure_taxonomy_npclassifier_03class,
    name_score = name
  ) %>%
  select(
    name_structure,
    value_structure,
    name_score,
    value_score
  ) %>%
  filter(!is.na(value_score))

superclass_pivoted <- superclass %>%
  pivot_longer(cols = 2:4, values_to = c("value_score")) %>%
  mutate(
    name_structure = "structure_taxonomy_npclassifier_02superclass",
    value_structure = structure_taxonomy_npclassifier_02superclass,
    name_score = name
  ) %>%
  select(
    name_structure,
    value_structure,
    name_score,
    value_score
  ) %>%
  filter(!is.na(value_score))

pathway_pivoted <- pathway %>%
  pivot_longer(cols = 2:4, values_to = c("value_score")) %>%
  mutate(
    name_structure = "structure_taxonomy_npclassifier_01pathway",
    value_structure = structure_taxonomy_npclassifier_01pathway,
    name_score = name
  ) %>%
  select(
    name_structure,
    value_structure,
    name_score,
    value_score
  ) %>%
  filter(!is.na(value_score))

final <- rbind(pathway_pivoted, superclass_pivoted, class_pivoted)

log_debug("exporting")
write_delim(
  x = final,
  file = file.path(pathDataProcessed, "jsd_full.tsv"),
  delim = "\t",
  na = ""
)

write_delim(
  x = class,
  file = file.path(pathDataProcessed, "jsd_class.tsv"),
  delim = "\t",
  na = ""
)

write_delim(
  x = superclass,
  file = file.path(pathDataProcessed, "jsd_superclass.tsv"),
  delim = "\t",
  na = ""
)

write_delim(
  x = pathway,
  file = file.path(pathDataProcessed, "jsd_pathway.tsv"),
  delim = "\t",
  na = ""
)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
