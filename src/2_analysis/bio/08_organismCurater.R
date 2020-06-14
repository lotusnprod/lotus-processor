#title: "Organisms (sanitized) curatoR"

#loading functions
source("functions.R")

#writing paths
##input
inpath <- "outputs/tables/2_sanitized/sanitizedOrganism.tsv.zip"

##output
outpath <- "outputs/tables/3_curated/curatedOrganism.tsv.zip"

outpath2 <-
  "outputs/tables/curatedOrganismsDifferentSpecies.tsv.zip"

#loading file
dataSanitizedOrganism <- read_delim(
  file = gzfile(inpath),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  select(
    organismOriginal,
    organismTranslated,
    organismSanitized,
    organism_database,
    organism_1_kingdom = kingdom,
    organism_2_phylum = phylum,
    organism_3_class = class,
    organism_4_order = order,
    organism_5_family = family,
    organism_6_genus = genus,
    organism_7_species = species
  ) %>%
  data.frame()

#curating taxonomy
##auto
dataCuratedOrganismAuto <-
  taxo_cleaning_auto(dfsel = dataSanitizedOrganism)

##manual
dataCuratedOrganism <-
  taxo_cleaning_manual(dfsel = dataCuratedOrganismAuto)

#outputing lower taxon
dataCuratedOrganism$organism_lowertaxon <-
  dataCuratedOrganism$organismCurated

#outputing differences in species names
diff <-
  dataCuratedOrganism %>% filter(organism_lowertaxon != organismCurated)

realDiff <-
  dataCuratedOrganism %>% filter(organism_lowertaxon != organismSanitized &
                                   !is.na(organism_6_genus)) %>%
  distinct(organismSanitized,
           organismCurated,
           .keep_all = TRUE) %>%
  group_by(organism_6_genus) %>%
  add_count() %>%
  ungroup() %>%
  select(organismSanitized,
         organismCurated,
         n)

#selecting
dataCuratedOrganism <- dataCuratedOrganism %>%
  select(
    organismOriginal,
    organismTranslated,
    organismSanitized,
    organismCurated,
    organism_database,
    organism_modified_taxonomy_auto,
    organism_modified_taxonomy_manual,
    organism_lowertaxon,
    organism_1_kingdom,
    organism_2_phylum,
    organism_3_class,
    organism_4_order,
    organism_5_family,
    organism_6_genus,
    organism_7_species
  )

#exporting
#curated
write.table(
  x = dataCuratedOrganism,
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

##differences
write.table(
  x = realDiff,
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

