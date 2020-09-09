cat("This script performs taxonomy alignment. \n")

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")

cat("... functions \n")
source("functions/log.R")
source("functions/helpers.R")
source("functions/bio.R")
source("2_curating/2_editing/organism/functions/manipulating_taxo.R")

log_debug(" Step 4")
cat("defining function \n")
taxo_cleaning_auto <- function(dfsel) {
  test <- dfsel %>%
    filter(!is.na(organismCleaned)) %>%
    distinct(organismCleaned,
             organismDbTaxo,
             .keep_all = TRUE)
  
  test$organism_1_kingdom <-
    y_as_na(x = test$organism_1_kingdom,
            y = "Not assigned")
  
  test$organism_2_phylum <-
    y_as_na(x = test$organism_2_phylum,
            y = "Not assigned")
  
  test$organism_3_class <-
    y_as_na(x = test$organism_3_class,
            y = "Not assigned")
  
  test$organism_4_order <-
    y_as_na(x = test$organism_4_order,
            y = "Not assigned")
  
  test$organism_5_family <-
    y_as_na(x = test$organism_5_family,
            y = "Not assigned")
  
  test$organism_6_genus <-
    y_as_na(x = test$organism_6_genus,
            y = "Not assigned")
  
  test$organism_7_species <-
    y_as_na(x = test$organism_7_species,
            y = "Not assigned")
  
  test$organism_8_variety <-
    y_as_na(x = test$organism_8_variety,
            y = "Not assigned")
  
  variety <- test %>%
    filter(!is.na(organism_8_variety)) %>%
    arrange(match(x = organismDbTaxo,
                  table =  c("NCBI"))) %>%
    group_by(organism_1_kingdom,
             organism_8_variety) %>%
    mutate(
      organism_2_phylum = na.locf(
        object = organism_2_phylum,
        na.rm = FALSE,
        fromLast = TRUE
      ),
      organism_3_class = na.locf(
        object = organism_3_class,
        na.rm = FALSE,
        fromLast = TRUE
      ),
      organism_4_order = na.locf(
        object = organism_4_order,
        na.rm = FALSE,
        fromLast = TRUE
      ),
      organism_5_family = na.locf(
        object = organism_5_family,
        na.rm = FALSE,
        fromLast = TRUE
      ),
      organism_6_genus = na.locf(
        object = organism_6_genus,
        na.rm = FALSE,
        fromLast = TRUE
      ),
      organism_7_species = na.locf(
        object = organism_7_species,
        na.rm = FALSE,
        fromLast = TRUE
      )
    ) %>%
    ungroup() %>%
    distinct(organism_8_variety,
             .keep_all = TRUE) %>%
    select(-organismOriginal, -organismCleaned, -organismTaxonId)
  
  variety_fill <- test %>%
    filter(!is.na(organism_8_variety)) %>%
    select(organismOriginal,
           organismCleaned,
           organismTaxonId,
           organism_8_variety)
  
  variety_full <- left_join(variety_fill, variety)
  
  unvariety <- test %>%
    filter(is.na(organism_8_variety))
  
  species_1 <- rbind(variety_full, unvariety)
  
  species <- test %>%
    filter(!is.na(organism_7_species)) %>%
    arrange(match(x = organismDbTaxo,
                  table =  c("NCBI"))) %>%
    group_by(organism_1_kingdom,
             organism_7_species) %>%
    mutate(
      organism_2_phylum = na.locf(
        object = organism_2_phylum,
        na.rm = FALSE,
        fromLast = TRUE
      ),
      organism_3_class = na.locf(
        object = organism_3_class,
        na.rm = FALSE,
        fromLast = TRUE
      ),
      organism_4_order = na.locf(
        object = organism_4_order,
        na.rm = FALSE,
        fromLast = TRUE
      ),
      organism_5_family = na.locf(
        object = organism_5_family,
        na.rm = FALSE,
        fromLast = TRUE
      ),
      organism_6_genus = na.locf(
        object = organism_6_genus,
        na.rm = FALSE,
        fromLast = TRUE
      ),
    ) %>%
    ungroup() %>%
    distinct(organism_7_species,
             .keep_all = TRUE) %>%
    select(-organismOriginal,
           -organismCleaned,
           -organismTaxonId,
           -organism_8_variety)
  
  species_fill <- test %>%
    filter(!is.na(organism_7_species)) %>%
    select(
      organismOriginal,
      organismCleaned,
      organismTaxonId,
      organism_8_variety,
      organism_7_species
    )
  
  species_full <- left_join(species_fill, species)
  
  unspecies <- test %>%
    filter(is.na(organism_7_species))
  
  genus_1 <- rbind(species_full, unspecies)
  
  genus <- genus_1 %>%
    filter(!is.na(organism_6_genus)) %>%
    arrange(match(x = organismDbTaxo,
                  table = c("NCBI"))) %>%
    group_by(organism_1_kingdom,
             organism_6_genus) %>%
    mutate(
      organism_2_phylum = na.locf(
        object = organism_2_phylum,
        na.rm = FALSE,
        fromLast = TRUE
      ),
      organism_3_class = na.locf(
        object = organism_3_class,
        na.rm = FALSE,
        fromLast = TRUE
      ),
      organism_4_order = na.locf(
        object = organism_4_order,
        na.rm = FALSE,
        fromLast = TRUE
      ),
      organism_5_family = na.locf(
        object = organism_5_family,
        na.rm = FALSE,
        fromLast = TRUE
      )
    ) %>%
    ungroup() %>%
    distinct(organism_6_genus,
             .keep_all = TRUE) %>%
    select(
      -organismOriginal,
      -organismCleaned,
      -organismTaxonId,
      -organism_8_variety,
      -organism_7_species
    )
  
  genus_fill <- genus_1 %>%
    filter(!is.na(organism_6_genus)) %>%
    select(
      organismOriginal,
      organismCleaned,
      organismTaxonId,
      organism_8_variety,
      organism_7_species,
      organism_6_genus
    )
  
  genus_full <- left_join(genus_fill, genus)
  
  ungenus <- genus_1 %>%
    filter(is.na(organism_6_genus))
  
  family_1 <- rbind(genus_full, ungenus)
  
  family <- family_1 %>%
    filter(!is.na(organism_5_family)) %>%
    arrange(match(x = organismDbTaxo,
                  table = c("NCBI"))) %>%
    group_by(organism_1_kingdom,
             organism_5_family) %>%
    mutate(
      organism_2_phylum = na.locf(
        object = organism_2_phylum,
        na.rm = FALSE,
        fromLast = TRUE
      ),
      organism_3_class = na.locf(
        object = organism_3_class,
        na.rm = FALSE,
        fromLast = TRUE
      ),
      organism_4_order = na.locf(
        object = organism_4_order,
        na.rm = FALSE,
        fromLast = TRUE
      )
    ) %>%
    ungroup() %>%
    distinct(organism_5_family,
             .keep_all = TRUE) %>%
    select(
      -organismOriginal,
      -organismCleaned,
      -organismTaxonId,
      -organism_8_variety,
      -organism_7_species,
      -organism_6_genus
    )
  
  family_fill <- family_1 %>%
    filter(!is.na(organism_5_family)) %>%
    select(
      organismOriginal,
      organismCleaned,
      organismTaxonId,
      organism_8_variety,
      organism_7_species,
      organism_6_genus,
      organism_5_family
    )
  
  family_full <- left_join(family_fill, family)
  
  unfamily <- family_1 %>%
    filter(is.na(organism_5_family))
  
  order_1 <- rbind(family_full, unfamily)
  
  order <- order_1 %>%
    filter(!is.na(organism_4_order)) %>%
    arrange(match(x = organismDbTaxo,
                  table = c("NCBI"))) %>%
    group_by(organism_1_kingdom,
             organism_4_order) %>%
    mutate(
      organism_2_phylum = na.locf(
        object = organism_2_phylum,
        na.rm = FALSE,
        fromLast = TRUE
      ),
      organism_3_class = na.locf(
        object = organism_3_class,
        na.rm = FALSE,
        fromLast = TRUE
      )
    ) %>%
    ungroup() %>%
    ungroup() %>%
    distinct(organism_4_order,
             .keep_all = TRUE) %>%
    select(
      -organismOriginal,
      -organismCleaned,
      -organismTaxonId,
      -organism_8_variety,
      -organism_7_species,
      -organism_6_genus,
      -organism_5_family
    )
  
  order_fill <- order_1 %>%
    filter(!is.na(organism_4_order)) %>%
    select(
      organismOriginal,
      organismCleaned,
      organismTaxonId,
      organism_8_variety,
      organism_7_species,
      organism_6_genus,
      organism_5_family,
      organism_4_order
    )
  
  order_full <- left_join(order_fill, order)
  
  unorder <- order_1 %>%
    filter(is.na(organism_4_order))
  
  class_1 <- rbind(order_full, unorder)
  
  class <- class_1 %>%
    filter(!is.na(organism_3_class)) %>%
    arrange(match(x = organismDbTaxo,
                  table = c("NCBI"))) %>%
    group_by(organism_1_kingdom,
             organism_3_class) %>%
    mutate(organism_2_phylum = na.locf(
      object = organism_2_phylum,
      na.rm = FALSE,
      fromLast = TRUE
    )) %>%
    ungroup() %>%
    distinct(organism_3_class,
             .keep_all = TRUE) %>%
    select(
      -organismOriginal,
      -organismCleaned,
      -organismTaxonId,
      -organism_8_variety,
      -organism_7_species,
      -organism_6_genus,
      -organism_5_family,
      -organism_4_order,
      
    )
  
  class_fill <- class_1 %>%
    filter(!is.na(organism_3_class)) %>%
    select(
      organismOriginal,
      organismCleaned,
      organismTaxonId,
      organism_8_variety,
      organism_7_species,
      organism_6_genus,
      organism_5_family,
      organism_4_order,
      organism_3_class
    )
  
  class_full <- left_join(class_fill, class)
  
  unclass <- class_1 %>%
    filter(is.na(organism_3_class))
  
  phylum_1 <- rbind(class_full, unclass)
  
  phylum <- phylum_1 %>%
    filter(!is.na(organism_2_phylum)) %>%
    arrange(match(x = organismDbTaxo,
                  table = c("NCBI"))) %>%
    group_by(organism_2_phylum) %>%
    mutate(organism_1_kingdom = na.locf(
      object = organism_1_kingdom,
      na.rm = FALSE,
      fromLast = TRUE
    )) %>%
    ungroup() %>%
    distinct(organism_2_phylum,
             .keep_all = TRUE) %>%
    select(
      -organismOriginal,
      -organismCleaned,
      -organismTaxonId,
      -organism_8_variety,
      -organism_7_species,
      -organism_6_genus,
      -organism_5_family,
      -organism_4_order,
      -organism_3_class
    )
  
  phylum_fill <- phylum_1 %>%
    filter(!is.na(organism_2_phylum)) %>%
    select(
      organismOriginal,
      organismCleaned,
      organismTaxonId,
      organism_8_variety,
      organism_7_species,
      organism_6_genus,
      organism_5_family,
      organism_4_order,
      organism_3_class,
      organism_2_phylum
    )
  
  phylum_full <- left_join(phylum_fill, phylum)
  
  unphylum <- phylum_1 %>%
    filter(is.na(organism_2_phylum))
  
  variety_tojoin <- variety_full %>%
    select(everything())
  
  kingdom_1 <- rbind(phylum_full, unphylum)
  
  kingdom_tojoin <- kingdom_1 %>%
    select(-organismOriginal)
  
  tojoin <- dfsel
  
  newdf = left_join(tojoin,
                    kingdom_tojoin,
                    by = c("organismCleaned" = "organismCleaned")) %>%
    mutate(organism_modified_taxonomy_auto = if_else(
      condition = paste(
        organism_1_kingdom.x,
        organism_2_phylum.x,
        organism_3_class.x,
        organism_4_order.x,
        organism_5_family.x,
        organism_6_genus.x,
        organism_7_species.x,
        organism_8_variety.x
      ) ==
        paste(
          organism_1_kingdom.y,
          organism_2_phylum.y,
          organism_3_class.y,
          organism_4_order.y,
          organism_5_family.y,
          organism_6_genus.y,
          organism_7_species.y,
          organism_8_variety.y
        ),
      true = "y",
      false = ""
    )) %>%
    select(
      organismOriginal,
      organismCleaned,
      organismDbTaxo = organismDbTaxo.y,
      organismDbTaxoQuality = organismDbTaxoQuality.y,
      organismTaxonId = organismTaxonId.y,
      organism_1_kingdom = organism_1_kingdom.y,
      organism_2_phylum = organism_2_phylum.y,
      organism_3_class = organism_3_class.y,
      organism_4_order = organism_4_order.y,
      organism_5_family = organism_5_family.y,
      organism_6_genus = organism_6_genus.y,
      organism_7_species = organism_7_species.y,
      organism_8_variety = organism_8_variety.y,
      organism_modified_taxonomy_auto
    ) %>%
    distinct(organismOriginal,
             organismCleaned,
             .keep_all = TRUE)
  
  newdf$organism_modified_taxonomy_auto <-
    y_as_na(x = newdf$organism_modified_taxonomy_auto,
            y = "")
  
  return(newdf)
}

cat("Step 4 \n")
cat("loading cleaned organisms \n")
dataCleanedOrganismManipulated <- read_delim(
  file = gzfile(pathDataInterimTablesCleanedOrganismTranslatedTable),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = FALSE
) %>%
  relocate(organismCleaned, .after = organismTaxonomy)

# curating taxonomy
## auto
# dataCuratedOrganismAuto <-
#   taxo_cleaning_auto(dfsel = dataCleanedOrganismManipulated)

## manual
# dataCuratedOrganism <-
#   taxo_cleaning_manual(dfsel = dataCleanedOrganismManipulated)

cat("keeping lowest taxon \n")
dataCuratedOrganism <- dataCleanedOrganismManipulated %>%
  mutate(organismCleaned =  as.character(apply(dataCleanedOrganismManipulated[7:15], 1, function(x)
    tail(na.omit(x), 1))))

dataCuratedOrganism$organismCleaned <-
  y_as_na(x = dataCuratedOrganism$organismCleaned,
          y = "character(0)")

dataCuratedOrganism$organismCleaned <-
  y_as_na(x = dataCuratedOrganism$organismCleaned,
          y = "NA")

cat("selecting \n")
dataCuratedOrganism[setdiff(
  x = c(
    "organismOriginal",
    "organismCleaned",
    "organismDbTaxo",
    "organismDbTaxoQuality",
    "organismTaxonIds",
    "organismTaxonRanks",
    "organismTaxonomy",
    "organism_1_kingdom",
    "organism_2_phylum",
    "organism_3_class",
    "organism_4_order",
    "organism_5_family",
    "organism_6_genus",
    "organism_7_species",
    "organism_8_quality"
  ),
  y = names(dataCuratedOrganism)
)] <- NA

dataCuratedOrganism <- dataCuratedOrganism %>%
  select(
    organismOriginal,
    organismCleaned,
    organismDbTaxo,
    organismDbTaxoQuality,
    # organismModifiedTaxonomyAuto = organism_modified_taxonomy_auto,
    # organismModifiedTaxonomyManual = organism_modified_taxonomy_manual,
    organismTaxonIds,
    organismTaxonRanks,
    organismTaxonomy,
    organism_1_kingdom,
    organism_2_phylum,
    organism_3_class,
    organism_4_order,
    organism_5_family,
    organism_6_genus,
    organism_7_species,
    organism_8_variety
  )

cat("exporting ... \n")
cat("pathDataInterimTablesCleanedOrganismFinal \n")
write.table(
  x = dataCuratedOrganism,
  file = gzfile(
    description = pathDataInterimTablesCleanedOrganismFinal,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

end <- Sys.time()

cat("Script finished in", end - start , "seconds \n")
