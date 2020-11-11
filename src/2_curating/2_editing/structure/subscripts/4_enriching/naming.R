cat("This script adds chemical names to structures dictionary \n")

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")

cat("... functions \n")
source("functions/analysis.R")
source("functions/helpers.R")

cat("loading files ... \n")
cat("...  counted structures \n")
structureCounted <- read_delim(
  file = gzfile(description = pathDataInterimTablesCleanedStructureStereoCounted),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
)

cat("keeping smiles only ... \n")
smilesDictionary <- structureCounted %>%
  distinct(smilesSanitized) %>%
  select(smiles = smilesSanitized)

cat("writing the smiles table \n")
write.table(
  x = smilesDictionary,
  file = gzfile(
    description = pathDataInterimTablesCleanedStructureSmiles,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

if (works_locally_only == FALSE) {
  cat("submitting to molconvert (traditional names) (no worries...running long) \n")
  system(
    command = paste(
      "bash",
      molconvertPath,
      "name:t",
      pathDataInterimTablesCleanedStructureSmiles,
      "-o",
      pathDataInterimTablesCleanedStructureSmiles_1,
      "-g"
    )
  )

  cat("submitting to molconvert (iupac) (no worries...running long) \n")
  system(
    command = paste(
      "bash",
      molconvertPath,
      "name:i",
      pathDataInterimTablesCleanedStructureSmiles,
      "-o",
      pathDataInterimTablesCleanedStructureSmiles_2,
      "-g"
    )
  )

  # cat("submitting to molconvert (common name) \n")
  # system(
  #   command = paste(
  #     "bash",
  #     molconvertPath,
  #     "name:common",
  #     pathDataInterimTablesCleanedStructureSmiles,
  #     "-o",
  #     pathDataInterimTablesCleanedStructureSmiles_3,
  #     "-g"
  #   )
  # )

  # cat("submitting to molconvert (all common name) \n")
  # system(
  #   command = paste(
  #     "bash",
  #     molconvertPath,
  #     "name:common,all",
  #     pathDataInterimTablesCleanedStructureSmiles,
  #     "-o",
  #     pathDataInterimTablesCleanedStructureSmiles_4,
  #     "-g"
  #   )
  # )
}

cat("loading files ... \n")
structureNamesTraditional <- read_delim(
  file = pathDataInterimTablesCleanedStructureSmiles_1,
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  col_names = FALSE,
  skip_empty_rows = FALSE,
  trim_ws = TRUE
) %>%
  select(structureCleaned_nameTraditional = X1)

structureNamesIupac <- read_delim(
  file = pathDataInterimTablesCleanedStructureSmiles_2,
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  col_names = FALSE,
  skip_empty_rows = FALSE,
  trim_ws = TRUE
) %>%
  select(structureCleaned_nameIupac = X1)

# structureNamesCommon <- read_delim(
#   file = pathDataInterimTablesCleanedStructureSmiles_3,
#   delim = "\t",
#   col_types = cols(.default = "c"),
#   escape_double = FALSE,
#   col_names = FALSE,
#   skip_empty_rows = FALSE,
#   trim_ws = TRUE
# ) %>%
#   select(structureCleanedNameCommon = X1)

# structureNamesCommonAll <- read_delim(
#   file = pathDataInterimTablesCleanedStructureSmiles_4,
#   delim = "\t",
#   col_types = cols(.default = "c"),
#   escape_double = FALSE,
#   col_names = FALSE,
#   skip_empty_rows = FALSE,
#   trim_ws = TRUE
# ) %>%
#   select(structureCleanedNameCommonAll = X1)

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


cat("ensuring directories exist \n")
cat("exporting ... \n")
cat(pathDataInterimTablesCleanedStructureNamed, "\n")

write.table(
  x = structureNamed,
  file = gzfile(
    description = pathDataInterimTablesCleanedStructureNamed,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

end <- Sys.time()

cat("Script finished in", format(end - start), "\n")