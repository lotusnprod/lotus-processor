#title: "Structure (sanitized) curatoR"

#writing paths
##input
inpath <- # "outputs/tables/2_sanitized/sanitizedStructure.tsv.zip"
  "outputs/tables/2_sanitized/sanitizedStructure.tsv" 

##output
outpath <- "outputs/tables/3_curated/curatedStructure.tsv.zip"
outpath2 <- "outputs/tables/problematicStructure.tsv.zip"

#loading file
dataSanitizedStructure <- read_delim(
  file = inpath, 
  # file = gzfile(inpath),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  select(structureTranslated,
         validatorLog,
         smiles = smiles_sanitized,
         inchi = inchi_sanitized,
         inchikey = inchikeySanitized,
         inchikey2D = shortikSanitized,
         molecularFormula = formulaSanitized,
         exactMass = exactmassSanitized,
         xlogP = xlogpSanitized
  ) %>%
  data.frame()

# selecting structures without problematic fragments
# problematic fragments are the ones where molecules 
# do have a fragment but do not have a corresponding log. 
# This means the fragments were not recognized by molVS. 
# Therefore we consider it problematic and drop the entry.

dataCuratedStructure <- dataSanitizedStructure %>%
  mutate(
    structureCurated = ifelse(
      test = !grepl(pattern = "\\.",
                    x = structureTranslated) |
        grepl(pattern = ".C6H3N3O7",
              x = structureTranslated) |
        #keeping picrate counter ion
        grepl(pattern = ".Fe",
              x = structureTranslated) |
        #keeping Fe counter ion
        grepl(pattern = ".[0-9]Fe",
              x = structureTranslated) |
        #keeping Fe counter ion
        grepl(pattern = ".Cu",
              x = structureTranslated) |
        #keeping Cu counter ion
        grepl(pattern = ".[0-9]Cu",
              x = structureTranslated) |
        #keeping Cu counter ion
        grepl(pattern = ".Ni",
              x = structureTranslated) |
        #keeping Ni counter ion
        grepl(pattern = ".[0-9]Ni",
              x = structureTranslated) |
        #keeping Ni counter ion
        grepl(pattern = ".Ga",
              x = structureTranslated) |
        #keeping Ga counter ion
        grepl(pattern = ".[0-9]Ga",
              x = structureTranslated) |
        #keeping Ga counter ion
        grepl(pattern = ".Mn",
              x = structureTranslated) |
        #keeping Mn counter ion
        grepl(pattern = ".[0-9]Mn",
              x = structureTranslated) |
        #keeping Mn counter ion
        grepl(pattern = ".Co",
              x = structureTranslated) |
        #keeping cyanocobalamine counter ion
        validatorLog != "[]",
      yes = inchi,
      no = NA
    )
  )

problematicStructure <- dataCuratedStructure %>%
  filter(is.na(structureCurated))

#exporting
#curated
write.table(
  x = dataCuratedStructure,
  file = gzfile(
    description = outpath,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

##problematics
write.table(
  x = problematicStructure,
  file = gzfile(
    description = outpath2,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
