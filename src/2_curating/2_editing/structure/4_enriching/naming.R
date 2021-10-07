source("r/log_debug.R")
log_debug("This script adds chemical names to structures dictionary")

start <- Sys.time()

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")

log_debug("... libraries")
library(dplyr)
library(readr)

log_debug("loading files ...")
log_debug("...  counted structures")
structureCounted <-
  read_delim(file = pathDataInterimTablesProcessedStructureStereoCounted)

log_debug("keeping smiles only ...")
smilesDictionary <- structureCounted %>%
  distinct(smilesSanitized) %>%
  select(smiles = smilesSanitized)

log_debug("writing the smiles table")
write_delim(
  x = smilesDictionary,
  delim = "\t",
  file = pathDataInterimTablesProcessedStructureSmiles
)

if (works_locally_only == FALSE) {
  log_debug("submitting to molconvert (traditional names) (no worries...running long)")
  system(
    command = paste(
      "bash",
      molconvertPath,
      "name:t",
      pathDataInterimTablesProcessedStructureSmiles,
      "-o",
      pathDataInterimTablesProcessedStructureSmiles_1,
      "-g"
    )
  )

  log_debug("submitting to molconvert (iupac) (no worries...running long)")
  system(
    command = paste(
      "bash",
      molconvertPath,
      "name:i",
      pathDataInterimTablesProcessedStructureSmiles,
      "-o",
      pathDataInterimTablesProcessedStructureSmiles_2,
      "-g"
    )
  )

  # log_debug("submitting to molconvert (common name)")
  # system(
  #   command = paste(
  #     "bash",
  #     molconvertPath,
  #     "name:common",
  #     pathDataInterimTablesProcessedStructureSmiles,
  #     "-o",
  #     pathDataInterimTablesProcessedStructureSmiles_3,
  #     "-g"
  #   )
  # )

  # log_debug("submitting to molconvert (all common name)")
  # system(
  #   command = paste(
  #     "bash",
  #     molconvertPath,
  #     "name:common,all",
  #     pathDataInterimTablesProcessedStructureSmiles,
  #     "-o",
  #     pathDataInterimTablesProcessedStructureSmiles_4,
  #     "-g"
  #   )
  # )

  log_debug("loading files ...")
  structureNamesTraditional <- read_delim(
    file = pathDataInterimTablesProcessedStructureSmiles_1,
    delim = "\t",
    col_types = cols(.default = "c"),
    escape_double = FALSE,
    col_names = FALSE,
    skip_empty_rows = FALSE,
    trim_ws = TRUE
  ) %>%
    select(structureCleaned_nameTraditional = X1)

  structureNamesIupac <- read_delim(
    file = pathDataInterimTablesProcessedStructureSmiles_2,
    delim = "\t",
    col_types = cols(.default = "c"),
    escape_double = FALSE,
    col_names = FALSE,
    skip_empty_rows = FALSE,
    trim_ws = TRUE
  ) %>%
    select(structureCleaned_nameIupac = X1)

  # structureNamesCommon <- read_delim(
  #   file = pathDataInterimTablesProcessedStructureSmiles_3,
  #   delim = "\t",
  #   col_types = cols(.default = "c"),
  #   escape_double = FALSE,
  #   col_names = FALSE,
  #   skip_empty_rows = FALSE,
  #   trim_ws = TRUE
  # ) %>%
  #   select(structureCleanedNameCommon = X1)

  # structureNamesCommonAll <- read_delim(
  #   file = pathDataInterimTablesProcessedStructureSmiles_4,
  #   delim = "\t",
  #   col_types = cols(.default = "c"),
  #   escape_double = FALSE,
  #   col_names = FALSE,
  #   skip_empty_rows = FALSE,
  #   trim_ws = TRUE
  # ) %>%
  #   select(structureCleanedNameCommonAll = X1)
}

if (works_locally_only == TRUE) {
  structureNamesTraditional <-
    data.frame(structureCleaned_nameTraditional = NA)

  structureNamesIupac <-
    data.frame(structureCleaned_nameIupac = NA)
}

smilesFilled <- bind_cols(
  smilesDictionary,
  structureNamesTraditional,
  structureNamesIupac
)

structureNamed <-
  left_join(structureCounted,
    smilesFilled,
    by = c("smilesSanitized" = "smiles")
  )

structureNamed_defined <- structureNamed %>%
  filter(count_unspecified_atomic_stereocenters == 0) %>%
  mutate_all(as.character)

structureNamed_undefined <- structureNamed %>%
  filter(count_unspecified_atomic_stereocenters != 0) %>%
  mutate(
    structureCleaned_nameTraditional = gsub(
      pattern = "(-)-",
      replacement = "",
      x = structureCleaned_nameTraditional,
      fixed = TRUE
    )
  ) %>%
  mutate(
    structureCleaned_nameTraditional = gsub(
      pattern = "(+)-",
      replacement = "",
      x = structureCleaned_nameTraditional,
      fixed = TRUE
    )
  ) %>%
  mutate(
    structureCleaned_nameTraditional = gsub(
      pattern = "(-)",
      replacement = "",
      x = structureCleaned_nameTraditional,
      fixed = TRUE
    )
  ) %>%
  mutate(
    structureCleaned_nameTraditional = gsub(
      pattern = "(+)",
      replacement = "",
      x = structureCleaned_nameTraditional,
      fixed = TRUE
    )
  ) %>%
  mutate_all(as.character)

structureNamed_cleaned <-
  bind_rows(structureNamed_defined, structureNamed_undefined)

log_debug("ensuring directories exist")
log_debug("exporting ...")
log_debug(pathDataInterimTablesProcessedStructureNamed)

write_delim(
  x = structureNamed_cleaned,
  delim = "\t",
  file = pathDataInterimTablesProcessedStructureNamed
)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
