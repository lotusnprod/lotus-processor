#title: "Open NP DB (sanitized) curatoR"

#loading functions
source("functions.R")

#writing paths
##inputs
###sanitized table
inpathSanitized <- "outputs/tables/2_sanitized/sanitizedTable.tsv.zip"

###curated bio
inpathOrganism <- "outputs/tables/3_curated/curatedOrganism.tsv.zip"

###(temporarily) sanitized chemo
inpathStructure <- "outputs/tables/4_additional_metadata/classyfire/first_run_classy_4.tsv"

###(TODO) sanitized ref
inpathReference <- "outputs/tables/3_curated/curatedReference.tsv.zip"

##output
outpath <- "outputs/tables/3_curated/curatedTable.tsv.zip"

#loading files
##sanitized table
dataSanitized <- read_delim(
  file = gzfile(inpathSanitized),
  col_types = cols(.default = "c"),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  data.frame()

##curated organisms
dataCuratedOrganism <- read_delim(
  file = gzfile(inpathOrganism),
  col_types = cols(.default = "c"),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

##structure
dataCuratedStructure <- read_delim(
  file = inpathStructure,
  col_types = cols(.default = "c"),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character) %>%
  data.frame()

##reference
dataCuratedReference <- read_delim(
  file = gzfile(inpathReference),
  col_types = cols(.default = "c"),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character) %>%
  data.frame()

#joining
##structure
dataCurated <- left_join(dataSanitized,
                         dataCuratedStructure)

#TO UPDATE WHEN CLASSY PROPERLY DONE

##organism
dataCurated <- left_join(dataCurated,
                         dataCuratedOrganism)

##ref
dataCurated <- left_join(dataCurated,
                         dataCuratedReference)

#list of InChIs (normally classyfied)

# inchi <- dataCurated %>%
#   filter(!is.na(structureInchi)) %>% 
#   distinct(structureInchi, .keep_all = TRUE) %>%
#   select(
#     inchi = structureInchi,
#     structure_chemontid,
#     structure_lowertaxon,
#     structure_01_kingdom,
#     structure_02_superclass,
#     structure_03_class,
#     structure_04_subclass,
#     structure_05_parent,
#     structure_06_subparent_1,
#     structure_07_subparent_2,
#     structure_08_subparent_3,
#     structure_09_subparent_4,
#     structure_10_subparent_5,
#     structure_11_subparent_6
#   )

#exporting
##curated table
write.table(
  x = dataCurated,
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

##InChI
# write.table(
#   x = inchi,
#   file = "outputs/tables/inchi.tsv",
#   row.names = FALSE,
#   quote = FALSE,
#   sep = "\t",
#   fileEncoding = "UTF-8"
# )
