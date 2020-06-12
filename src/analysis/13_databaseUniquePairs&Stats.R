#title: "Unique pairs compileR"

#loading functions
source("functions.R")

inpath <- "outputs/tables/3_curated/curatedTable.tsv.zip"

#outpaths
outpathOpen <- "outputs/tables/openDbPairs.tsv.zip"

outpathInhouse <- "outputs/tables/inhouseDbPairs.tsv.zip"

outpathDNP <- "outputs/tables/dnpDbPairs.tsv.zip"

outpathKingdom <- "outputs/tables/structuresByKingdom.tsv"

outpathSpecies <- "outputs/tables/uniqueStructuresBySpecies.tsv"

outpathWidespread <- "outputs/tables/widespreadStructures.tsv"

outpathMismatched <- "outputs/tables/mismatchedGenera.tsv"

outpathRedundancy <- "outputs/tables/redundancyTable.tsv"

#loading files
##fullDB
inhouseDb <- read_delim(
  file = gzfile(inpath),
  col_types = cols(.default = "c"),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  data.frame() %>%
  select(
    database,
    nameSanitized,
    structureOriginalInchi,
    structureOriginalSmiles,
    structureOriginalNominal,
    organismOriginal,
    referenceOriginal,
    organismCurated,
    structureCurated,
    sanitizedTitle,
    sanitizedDoi,
    inchi,
    smiles,
    inchikey,
    inchikey2D,
    molecularFormula,
    exactMass,
    xlogP,
    chemo_lowertaxon,
    chemo_01_kingdom,
    chemo_02_superclass,
    chemo_03_class,
    chemo_04_subclass,
    chemo_05_parent,
    organism_lowertaxon,
    organism_1_kingdom,
    organism_2_phylum,
    organism_3_class,
    organism_4_order,
    organism_5_family,
    organism_6_genus,
    organism_7_species
  )
#filter(!is.na(organism_lowertaxon))

##openDB
openDb <- inhouseDb %>%
  filter(database != "dnp_1")

##DNP
dnpDb <- inhouseDb %>%
  filter(database == "dnp_1")

#unique
##biological taxa
###open NP DB
openDbOrganism <- distinct_biosources(x = openDb)

###inhouseDB
inhouseDbOrganism <- distinct_biosources(x = inhouseDb)

###DNP
dnpDbOrganism <- distinct_biosources(x = dnpDb)

##structures
###open NP DB
openDbStructure <- openDb %>%
  filter(!is.na(structureCurated)) %>%
  distinct(structureCurated, .keep_all = TRUE)

###inhouseDB
inhouseDbStructure <- inhouseDb %>%
  filter(!is.na(structureCurated)) %>%
  distinct(structureCurated, .keep_all = TRUE)

###DNP
dnpDbStructure <- dnpDb %>%
  filter(!is.na(structureCurated)) %>%
  distinct(structureCurated, .keep_all = TRUE)

##pairs
###open NP DB
openDbPairs <- distinct_pairs(x = openDb)

###inhouseDB
inhouseDbPairs <- distinct_pairs(x = inhouseDb)

###DNP
dnpDbPairs <- distinct_pairs(dnpDb)

#writing tabular stats
##species by kingdom
inhouseSpeciesByKingdom <- inhouseDbPairs %>%
  group_by(organism_1_kingdom) %>%
  distinct(organism_7_species, .keep_all = TRUE) %>%
  count(organism_1_kingdom) %>%
  ungroup() %>%
  mutate(speciesPercent = 100 * n / sum(n)) %>%
  select(kingdom = organism_1_kingdom,
         species = n,
         speciesPercent) %>%
  arrange(desc(speciesPercent)) %>%
  head(10)

##structures by class
inhouseStructuresByClass <- inhouseDbPairs %>%
  group_by(chemo_03_class) %>%
  distinct(structureCurated, .keep_all = TRUE) %>%
  count(chemo_03_class) %>%
  ungroup() %>%
  mutate(structuresPercent = 100 * n / sum(n)) %>%
  select(class = chemo_03_class,
         structures = n,
         structuresPercent) %>%
  arrange(desc(structuresPercent)) %>%
  head(10)

##structures by kingdom
inhouseStructuresByOrganismKingdom <- inhouseDbPairs %>%
  group_by(organism_1_kingdom) %>%
  distinct(structureCurated, organism_1_kingdom, .keep_all = TRUE) %>%
  count(organism_1_kingdom) %>%
  ungroup() %>%
  mutate(structuresPercent = 100 * n / sum(n)) %>%
  select(kingdom = organism_1_kingdom,
         structures = n,
         structuresPercent) %>%
  arrange(desc(structuresPercent)) %>%
  head(10)

##unique structures per kingdom
inhouseUniqueStructuresPerKingdom <- inhouseDbPairs %>%
  group_by(structureCurated) %>%
  add_count(structureCurated) %>%
  ungroup() %>%
  filter(n == 1) %>%
  select(-n) %>%
  group_by(organism_1_kingdom) %>%
  distinct(structureCurated, organism_1_kingdom, .keep_all = TRUE) %>%
  count(organism_1_kingdom) %>%
  ungroup() %>%
  select(kingdom = organism_1_kingdom,
         specificStructures = n) %>%
  arrange(desc(specificStructures)) %>%
  head(10)

##structures by kingdom
inhouseStructuresByKingdom <-
  full_join(inhouseSpeciesByKingdom,
            inhouseStructuresByOrganismKingdom) %>%
  mutate(strucuresPerSpecies = structures / species)

inhouseStructuresByKingdom <-
  full_join(inhouseStructuresByKingdom,
            inhouseUniqueStructuresPerKingdom) %>%
  mutate(structuresSpecificity = 100 * specificStructures / structures) %>%
  select(1, 2, 3, 4, 7, 8, 5, 6)

##unique structures per species
inhouseUniqueStructuresPerSpecies <- inhouseDbPairs %>%
  group_by(structureCurated) %>%
  add_count(structureCurated) %>%
  ungroup() %>%
  filter(n == 1) %>%
  select(-n) %>%
  group_by(organism_7_species) %>%
  distinct(structureCurated, organism_7_species, .keep_all = TRUE) %>%
  count(organism_7_species) %>%
  ungroup() %>%
  select(species = organism_7_species,
         specificStructures = n) %>%
  arrange(desc(specificStructures)) %>%
  head(10)

##widespread metabolites
openDbWidespread <- openDbPairs %>%
  group_by(structureCurated) %>%
  filter(!is.na(organism_1_kingdom)) %>%
  distinct(organism_1_kingdom, .keep_all = TRUE) %>%
  add_count() %>%
  ungroup() %>%
  filter(n == 7) %>%
  arrange(structureCurated)

##word(species,1) != genus
mismatchedGenera <- inhouseDbOrganism %>%
  filter(word(organism_7_species, 1) != organism_6_genus)

##redundancy table
redundancydf  <- inhouseDb %>%
  filter(!is.na(organismCurated) &
           !is.na(structureCurated) &
           !is.na(organism_7_species)) %>%
  distinct(organismCurated,
           structureCurated,
           database,
           .keep_all = TRUE) %>%
  group_by(organismCurated, structureCurated) %>%
  add_count() %>%
  ungroup() %>%
  filter(n >= 5) %>%
  select(
    database,
    structureOriginalInchi,
    structureOriginalSmiles,
    structureOriginalNominal,
    organismOriginal,
    organismCurated,
    structureCurated,
    n
  )

##stereo trial
# most_studied_families <- inhouseDbPairs %>%
#   group_by(organism_5_family) %>%
#   count() %>%
#   ungroup() %>%
#   arrange(desc(n)) %>%
#   head(10)
#
# #Asteraceae
# most_studied_Asteraceae <- inhouseDbPairs %>%
#   filter(organism_5_family == "Asteraceae") %>%
#   group_by(organism_6_genus) %>%
#   count() %>%
#   ungroup() %>%
#   arrange(desc(n)) %>%
#   head(10)
#
# #Artemisia
# most_studied_Artemisia <- inhouseDbPairs %>%
#   filter(organism_6_genus == "Artemisia") %>%
#   filter(!is.na(organism_7_species)) %>%
#   group_by(organism_7_species) %>%
#   count() %>%
#   ungroup() %>%
#   arrange(desc(n)) %>%
#   head(10)
#
# a_annua_3D <- inhouseDbPairs %>%
#   filter(organism_7_species == "Artemisia annua",!is.na(structureCurated)) %>%
#   select(structure_inchikey) %>%
#   data.table()
#
# a_capillaris_3D <- inhouseDbPairs %>%
#   filter(organism_7_species == "Artemisia capillaris",!is.na(structureCurated)) %>%
#   select(structure_inchikey) %>%
#   data.table()
#
# a_annua_2D <- inhouseDbPairs %>%
#   filter(organism_7_species == "Artemisia annua",!is.na(structureCurated)) %>%
#   select(structure_shortinchikey) %>%
#   data.table()
#
# a_capillaris_2D <- inhouseDbPairs %>%
#   filter(organism_7_species == "Artemisia capillaris",!is.na(structureCurated)) %>%
#   select(structure_shortinchikey) %>%
#   data.table()
#
# overlap_artemisia_3D <- inner_join(a_annua_3D, a_capillaris_3D)
# overlap_artemisia_2D <- inner_join(a_annua_2D, a_capillaris_2D)
#
# overlap_percent_artemisia_3D <-
#   100 * nrow(overlap_artemisia_3D) / nrow(full_join(a_annua_3D, a_capillaris_3D))
#
# overlap_percent_artemisia_2D <-
#   100 * nrow(overlap_artemisia_2D) / nrow(full_join(a_annua_2D, a_capillaris_2D))
#
# ratio_artemisia_3D_2D <-
#   100 * overlap_percent_artemisia_3D / overlap_percent_artemisia_2D
#
# #Senecio
# most_studied_Senecio <- inhouseDbPairs %>%
#   filter(organism_6_genus == "Senecio") %>%
#   filter(!is.na(organism_7_species)) %>%
#   group_by(organism_7_species) %>%
#   count() %>%
#   ungroup() %>%
#   arrange(desc(n)) %>%
#   head(10)
#
# s_scandens_3D <- inhouseDbPairs %>%
#   filter(organism_7_species == "Senecio scandens",!is.na(structureCurated)) %>%
#   select(structure_inchikey) %>%
#   data.table()
#
# s_squalidus_3D <- inhouseDbPairs %>%
#   filter(organism_7_species == "Senecio squalidus",!is.na(structureCurated)) %>%
#   select(structure_inchikey) %>%
#   data.table()
#
# s_scandens_2D <- inhouseDbPairs %>%
#   filter(organism_7_species == "Senecio scandens",!is.na(structureCurated)) %>%
#   select(structure_shortinchikey) %>%
#   data.table()
#
# s_squalidus_2D <- inhouseDbPairs %>%
#   filter(organism_7_species == "Senecio squalidus",!is.na(structureCurated)) %>%
#   select(structure_shortinchikey) %>%
#   data.table()
#
# overlap_senecio_3D <- inner_join(s_scandens_3D, s_squalidus_3D)
# overlap_senecio_2D <- inner_join(s_scandens_2D, s_squalidus_2D)
#
# overlap_percent_senecio_3D <-
#   100 * nrow(overlap_senecio_3D) / nrow(full_join(s_scandens_3D, s_squalidus_3D))
#
# overlap_percent_senecio_2D <-
#   100 * nrow(overlap_senecio_2D) / nrow(full_join(s_scandens_2D, s_squalidus_2D))
#
# ratio_senecio_3D_2D <-
#   100 * overlap_percent_senecio_3D / overlap_percent_senecio_2D
#
# #Fabaceae
# most_studied_Fabaceae <- inhouseDbPairs %>%
#   filter(organism_5_family == "Fabaceae") %>%
#   group_by(organism_6_genus) %>%
#   count() %>%
#   ungroup() %>%
#   arrange(desc(n)) %>%
#   head(10)
#
# #Glycyrrhiza
# most_studied_Glycyrrhiza <- inhouseDbPairs %>%
#   filter(organism_6_genus == "Glycyrrhiza") %>% #Trifolium first but less diversity
#   filter(!is.na(organism_7_species)) %>%
#   group_by(organism_7_species) %>%
#   count() %>%
#   ungroup() %>%
#   arrange(desc(n)) %>%
#   head(10)
#
# g_uralensis_3D <- inhouseDbPairs %>%
#   filter(organism_7_species == "Glycyrrhiza uralensis",!is.na(structureCurated)) %>%
#   select(structure_inchikey) %>%
#   data.table()
#
# g_glabra_3D <- inhouseDbPairs %>%
#   filter(organism_7_species == "Glycyrrhiza glabra",!is.na(structureCurated)) %>%
#   select(structure_inchikey) %>%
#   data.table()
#
# g_uralensis_2D <- inhouseDbPairs %>%
#   filter(organism_7_species == "Glycyrrhiza uralensis",!is.na(structureCurated)) %>%
#   select(structure_shortinchikey) %>%
#   data.table()
#
# g_glabra_2D <- inhouseDbPairs %>%
#   filter(organism_7_species == "Glycyrrhiza glabra",!is.na(structureCurated)) %>%
#   select(structure_shortinchikey) %>%
#   data.table()
#
# overlap_glycyrrhiza_3D <- inner_join(g_uralensis_3D, g_glabra_3D)
# overlap_glycyrrhiza_2D <- inner_join(g_uralensis_2D, g_glabra_2D)
#
# overlap_percent_glycyrrhiza_3D <-
#   100 * nrow(overlap_glycyrrhiza_3D) / nrow(full_join(g_uralensis_3D, g_glabra_3D))
#
# overlap_percent_glycyrrhiza_2D <-
#   100 * nrow(overlap_glycyrrhiza_2D) / nrow(full_join(g_uralensis_2D, g_glabra_2D))
#
# ratio_glycyrrhiza_3D_2D <-
#   100 * overlap_percent_glycyrrhiza_3D / overlap_percent_glycyrrhiza_2D
#
# #Sophora
# most_studied_Sophora <- inhouseDbPairs %>%
#   filter(organism_6_genus == "Sophora") %>%
#   filter(!is.na(organism_7_species)) %>%
#   group_by(organism_7_species) %>%
#   count() %>%
#   ungroup() %>%
#   arrange(desc(n)) %>%
#   head(10)
#
# s_flavescens_3D <- inhouseDbPairs %>%
#   filter(organism_7_species == "Sophora flavescens",!is.na(structureCurated)) %>%
#   select(structure_inchikey) %>%
#   data.table()
#
# s_japonica_3D <- inhouseDbPairs %>%
#   filter(organism_7_species == "Sophora japonica",!is.na(structureCurated)) %>%
#   select(structure_inchikey) %>%
#   data.table()
#
# s_flavescens_2D <- inhouseDbPairs %>%
#   filter(organism_7_species == "Sophora flavescens",!is.na(structureCurated)) %>%
#   select(structure_shortinchikey) %>%
#   data.table()
#
# s_japonica_2D <- inhouseDbPairs %>%
#   filter(organism_7_species == "Sophora japonica",!is.na(structureCurated)) %>%
#   select(structure_shortinchikey) %>%
#   data.table()
#
# overlap_sophora_3D <- inner_join(s_flavescens_3D, s_japonica_3D)
# overlap_sophora_2D <- inner_join(s_flavescens_2D, s_japonica_2D)
#
# overlap_percent_sophora_3D <-
#   100 * nrow(overlap_sophora_3D) / nrow(full_join(s_flavescens_3D, s_japonica_3D))
#
# overlap_percent_sophora_2D <-
#   100 * nrow(overlap_sophora_2D) / nrow(full_join(s_flavescens_2D, s_japonica_2D))
#
# ratio_sophora_3D_2D <-
#   100 * overlap_percent_sophora_3D / overlap_percent_sophora_2D
#
# artemisia_3D <- inhouseDbPairs %>%
#   filter(organism_6_genus == "Artemisia",!is.na(structureCurated)) %>%
#   select(structure_inchikey) %>%
#   data.table()
#
# artemisia_2D <- inhouseDbPairs %>%
#   filter(organism_6_genus == "Artemisia",!is.na(structureCurated)) %>%
#   select(structure_shortinchikey) %>%
#   data.table()
#
# senecio_3D <- inhouseDbPairs %>%
#   filter(organism_6_genus == "Senecio",!is.na(structureCurated)) %>%
#   select(structure_inchikey) %>%
#   data.table()
#
# senecio_2D <- inhouseDbPairs %>%
#   filter(organism_6_genus == "Senecio",!is.na(structureCurated)) %>%
#   select(structure_shortinchikey) %>%
#   data.table()
#
# overlap_asteraceae_3D <- inner_join(artemisia_3D, senecio_3D)
# overlap_asteraceae_2D <- inner_join(artemisia_2D, senecio_2D)
#
# overlap_percent_asteraceae_3D <-
#   100 * nrow(overlap_asteraceae_3D) / nrow(full_join(artemisia_3D, senecio_3D))
#
# overlap_percent_asteraceae_2D <-
#   100 * nrow(overlap_asteraceae_2D) / nrow(full_join(artemisia_2D, senecio_2D))
#
# ratio_asteraceae_3D_2D <-
#   100 * overlap_percent_asteraceae_3D / overlap_percent_asteraceae_2D
#
# glycyrrhiza_3D <- inhouseDbPairs %>%
#   filter(organism_6_genus == "Glycyrrhiza",!is.na(structureCurated)) %>%
#   select(structure_inchikey) %>%
#   data.table()
#
# glycyrrhiza_2D <- inhouseDbPairs %>%
#   filter(organism_6_genus == "Glycyrrhiza",!is.na(structureCurated)) %>%
#   select(structure_shortinchikey) %>%
#   data.table()
#
# sophora_3D <- inhouseDbPairs %>%
#   filter(organism_6_genus == "Sophora",!is.na(structureCurated)) %>%
#   select(structure_inchikey) %>%
#   data.table()
#
# sophora_2D <- inhouseDbPairs %>%
#   filter(organism_6_genus == "Sophora",!is.na(structureCurated)) %>%
#   select(structure_shortinchikey) %>%
#   data.table()
#
# overlap_fabaceae_3D <- inner_join(glycyrrhiza_3D, sophora_3D)
# overlap_fabaceae_2D <- inner_join(glycyrrhiza_2D, sophora_2D)
#
# overlap_percent_fabaceae_3D <-
#   100 * nrow(overlap_fabaceae_3D) / nrow(full_join(glycyrrhiza_3D, sophora_3D))
#
# overlap_percent_fabaceae_2D <-
#   100 * nrow(overlap_fabaceae_2D) / nrow(full_join(glycyrrhiza_2D, sophora_2D))
#
# ratio_fabaceae_3D_2D <-
#   100 * overlap_percent_fabaceae_3D / overlap_percent_fabaceae_2D

# get_stereo_ratio <-
#   function(data,
#            biological_level,
#            biological_filter_value = NULL,
#            biological_filter_level = NULL,
#            chemical_filter_value = NULL,
#            chemical_filter_level = NULL) {
#     if (biological_filter_level == biological_level)
#       unique_ik_data <-
#         unique(data[c("structure_inchikey",
#                       biological_filter_level,
#                       chemical_filter_level)])
#
#     if (!is.null(biological_filter_level) &
#         biological_filter_level != biological_level)
#       unique_ik_data <-
#         unique(data[c(
#           "structure_inchikey",
#           biological_level,
#           biological_filter_level,
#           chemical_filter_level
#         )])
#
#     if (is.null(biological_filter_level))
#       unique_ik_data <-
#         unique(data[c("structure_inchikey",
#                       biological_level,
#                       chemical_filter_level)])
#
#     sik_data <- data %>%
#       distinct(structure_inchikey, structure_shortinchikey)
#
#     unique_ik_sik_data <- left_join(unique_ik_data, sik_data)
#
#     table_3D <- data.frame(unique_ik_sik_data)
#
#     if (!is.null(biological_filter_value))
#       table_3D <-
#       table_3D[table_3D[, biological_filter_level] %in% biological_filter_value, ]
#
#     if (!is.null(chemical_filter_value))
#       table_3D <-
#       table_3D[table_3D[, chemical_filter_level] %in% chemical_filter_value, ]
#
#     table_3D <-
#       table_3D[!is.na(table_3D[, biological_level]), ]
#
#     table_3D <- table_3D[c("structure_inchikey", biological_level)]
#
#     int_3D <-
#       mat.or.vec(nr = length(unique(table_3D[, biological_level])),
#                  nc = length(unique(table_3D[, biological_level])))
#
#     colnames(int_3D) <- unique(table_3D[, biological_level])
#     rownames(int_3D) <- unique(table_3D[, biological_level])
#
#     un_3D <- int_3D
#
#     for (col in colnames(int_3D)) {
#       for (row in colnames(int_3D)) {
#         int_3D[row, col] <-
#           length(intersect(table_3D[table_3D[, biological_level] == col, "structure_inchikey"],
#                            table_3D[table_3D[, biological_level] == row, "structure_inchikey"]))
#       }
#     }
#
#     for (col in colnames(un_3D)) {
#       for (row in colnames(un_3D)) {
#         un_3D[row, col] <-
#           sum(length(union(table_3D[table_3D[, biological_level] == col, "structure_inchikey"],
#                            table_3D[table_3D[, biological_level] == row, "structure_inchikey"])),
#               ifelse(
#                 test = length(union(table_3D[table_3D[, biological_level] == col, "structure_inchikey"],
#                                     table_3D[table_3D[, biological_level] == row, "structure_inchikey"])) == 2,
#                 yes = length(union(table_3D[table_3D[, biological_level] == col, "structure_inchikey"],
#                                    table_3D[table_3D[, biological_level] == row, "structure_inchikey"])[[2]]),
#                 no = 0
#               ))
#       }
#     }
#
#     int_3D <- data.frame(int_3D)
#     un_3D <- data.frame(un_3D)
#     ratio_3D <- data.frame(int_3D / un_3D * 100)
#
#     table_2D <- data.frame(unique_ik_sik_data)
#
#     if (!is.null(biological_filter_value))
#       table_2D <-
#       table_2D[table_2D[, biological_filter_level] %in% biological_filter_value, ]
#
#     if (!is.null(chemical_filter_value))
#       table_2D <-
#       table_2D[table_2D[, chemical_filter_level] %in% chemical_filter_value, ]
#
#     table_2D <-
#       table_2D[!is.na(table_2D[, biological_level]), ]
#
#     table_2D <-
#       table_2D[c("structure_shortinchikey", biological_level)]
#
#     int_2D <-
#       mat.or.vec(nr = length(unique(table_2D[, biological_level])),
#                  nc = length(unique(table_2D[, biological_level])))
#
#     colnames(int_2D) <- unique(table_2D[, biological_level])
#     rownames(int_2D) <- unique(table_2D[, biological_level])
#
#     un_2D <- int_2D
#
#     for (col in colnames(int_2D)) {
#       for (row in colnames(int_2D)) {
#         int_2D[row, col] <-
#           length(intersect(table_2D[table_2D[, biological_level] == col, "structure_shortinchikey"],
#                            table_2D[table_2D[, biological_level] ==  row, "structure_shortinchikey"]))
#       }
#     }
#
#     for (col in colnames(un_2D)) {
#       for (row in colnames(un_2D)) {
#         un_2D[row, col] <-
#           sum(length(union(table_2D[table_2D[, biological_level] == col, "structure_shortinchikey"],
#                            table_2D[table_2D[, biological_level] == row, "structure_shortinchikey"])),
#               ifelse(
#                 test = length(union(table_2D[table_2D[, biological_level] == col, "structure_shortinchikey"],
#                                     table_2D[table_2D[, biological_level] == row, "structure_shortinchikey"])) == 2,
#                 yes =  length(union(table_2D[table_2D[, biological_level] == col, "structure_shortinchikey"],
#                                     table_2D[table_2D[, biological_level] == row, "structure_shortinchikey"])[[2]]),
#                 no = 0
#               ))
#       }
#     }
#
#     int_2D <- data.frame(int_2D)
#     un_2D <- data.frame(un_2D)
#     ratio_2D <- data.frame(int_2D / un_2D * 100)
#
#     ratio_stereo <- ratio_3D / ratio_2D * 100
#
#     result <-
#       list(int_2D, int_3D, un_2D, un_3D, ratio_2D, ratio_3D, ratio_stereo)
#
#     names(result) <-
#       c(
#         "intersect2D",
#         "intersect3D",
#         "union2D",
#         "union3D",
#         "ratio2D",
#         "ratio3D",
#         "ratiostereo"
#       )
#
#     return(result)
#   }

#list of good examples
##at organism_6_genus level: structure_05_parent == "Carotenoids"
##at organism_6_genus level: structure_05_parent == "Monosaccharides"
###(Armoracia, Leipidium, Erysimum (all Brassicaceae) highly different but not Brassica (intermediary))
##at organism_7_species level: structure_04_subclass == "Terpene lactones"
###(almost no common structures except quassinoids -> led to ref issue discovery)

organismFilterDf <- openDbPairs %>%
  filter(!is.na(organism_7_species)) %>%
  filter(chemo_04_subclass == "Terpene lactones") %>%
  group_by(organism_7_species) %>%
  add_count() %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  distinct(organism_7_species) %>%
  head(24)

organismFilter <- organismFilterDf$organism_7_species

structureFilterDf <- openDbPairs %>%
  filter(!is.na(chemo_04_subclass)) %>%
  group_by(chemo_04_subclass) %>%
  add_count() %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  distinct(chemo_04_subclass) %>%
  filter(chemo_04_subclass == "Terpene lactones")

structureFilter <- structureFilterDf$chemo_04_subclass

# stereolist <- get_stereo_ratio(
#   data = open_db_pairs,
#   biological_level = "organism_7_species",
#   biological_filter_value = organism_filter,
#   biological_filter_level = "organism_7_species",
#   chemical_filter_value = structure_filter,
#   chemical_filter_level = "chemo_04_subclass"
# )

#exporting
##open
write.table(
  x = openDbPairs,
  file = gzfile(
    description = outpathOpen,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

##inhouse
write.table(
  x = inhouseDbPairs,
  file = gzfile(
    description = outpathInhouse,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

##dnp
write.table(
  x = dnpDbPairs,
  file = gzfile(
    description = outpathDNP,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

##stats
###structures by kingdom
write.table(
  x = inhouseStructuresByKingdom,
  file = outpathKingdom,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

###unique structures per species
write.table(
  x = inhouseUniqueStructuresPerSpecies,
  file = outpathSpecies,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

###widespread metabolites
write.table(
  x = openDbWidespread,
  file = outpathWidespread,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

###mismatched genera
write.table(
  x = mismatchedGenera,
  file = outpathMismatched,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

###redundancy table
write.table(
  x = redundancydf,
  file = outpathRedundancy,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
