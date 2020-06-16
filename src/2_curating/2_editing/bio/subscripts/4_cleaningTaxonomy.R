# title: "treating bio"

# setting working directory
setwd("~/GitLab/opennaturalproductsdb/src/")

# loading paths
source("paths.R")

# loading functions
source("functions.R")

#curating taxonomy
##auto
dataCuratedOrganismAuto <-
  taxo_cleaning_auto(dfsel = dataCleanedOrganismManipulated)

##manual
dataCuratedOrganism <-
  taxo_cleaning_manual(dfsel = dataCuratedOrganismAuto)

#outputing lowest taxon
dataCuratedOrganism$organismLowestTaxon <-
  dataCuratedOrganism$organismCurated

#outputing differences in species names
diff <-
  dataCuratedOrganism %>% filter(organismLowestTaxon != organismCurated)

realDiff <-
  dataCuratedOrganism %>% filter(organismLowestTaxon != organismCleaned &
                                   !is.na(organism_6_genus)) %>%
  distinct(organismCleaned,
           organismCurated,
           .keep_all = TRUE) %>%
  group_by(organism_6_genus) %>%
  add_count() %>%
  ungroup() %>%
  select(organismCleaned,
         organismCurated,
         n)

#selecting
dataCuratedOrganism <- dataCuratedOrganism %>%
  select(
    organismOriginal,
    organismCleaned,
    organismCurated,
    organismLowestTaxon,
    #duplicate of organism curated, choose
    organismDbTaxo,
    organismDbTaxoQuality,
    organismModifiedTaxonomyAuto = organism_modified_taxonomy_auto,
    organismModifiedTaxonomyManual = organism_modified_taxonomy_manual,
    organismTaxonId,
    organism_1_kingdom,
    organism_2_phylum,
    organism_3_class,
    organism_4_order,
    organism_5_family,
    organism_6_genus,
    organism_7_species,
    organism_8_variety
  )

#exporting
#curated
write.table(
  x = dataCuratedOrganism,
  file = gzfile(
    description = pathCuratedOrganism,
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
    description = pathCuratedOrganismRealDiff,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
