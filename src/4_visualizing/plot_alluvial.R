cat("This script plots the alluvial plot ... \n")

start <- Sys.time()
library(data.table)
library(ggalluvial)
library(ggfittext)
library(splitstackshape)
library(tidyverse)
source(file = "paths.R")
source(file = "r/vroom_safe.R")

cat("... open DB \n")
openDbMetaValidated <-
  vroom_read_safe(path = pathDataInterimTablesAnalysedPlatinum) %>%
  filter(
    !is.na(structureCleanedInchikey) &
      !is.na(organismCleaned) &
      !is.na(referenceCleanedDoi)
  ) %>%
  mutate(validation = "validated") %>%
  tibble()

openDbMaximal <-
  vroom_read_safe(path = pathDataInterimTablesCuratedTableMaximal) %>%
  tibble()

full <- left_join(openDbMaximal, openDbMetaValidated) %>%
  mutate(
    validation = ifelse(
      test = !is.na(validation),
      yes = validation,
      no = "not_validated"
    ),
    cleaned_reference = ifelse(
      test = !is.na(referenceCleanedDoi) |
        !is.na(referenceCleanedPmcid) |
        !is.na(referenceCleanedPmid),
      yes = "reference_Yes",
      no = "reference_No"
    ),
    cleaned_organism = ifelse(
      test = !is.na(organismCleaned),
      yes = "organism_Yes",
      no = "organism_No"
    ),
    cleaned_structure = ifelse(
      test = !is.na(structureCleanedSmiles) |
        !is.na(structureCleanedInchi) |
        !is.na(structureCleanedInchikey),
      yes = "structure_Yes",
      no = "structure_No"
    ),
  ) %>%
  distinct(
    database,
    organismType,
    organismValue,
    structureType,
    structureValue,
    referenceType,
    referenceValue,
    organismCleaned,
    structureCleanedInchikey,
    referenceCleanedTitle,
    cleaned_structure,
    cleaned_organism,
    cleaned_reference,
    validation,
  ) %>%
  pivot_longer(
    cols = 12:14,
    values_drop_na = TRUE
  ) %>%
  select(
    database,
    structure_originalValue = structureValue,
    structure_originalType = structureType,
    organism_originalValue = organismValue,
    organism_originalType = organismType,
    reference_originalValue = referenceValue,
    reference_originalType = referenceType,
    cleanedType = name,
    cleanedValue = value,
    validation
  ) %>%
  distinct() %>%
  data.frame()

ready_1 <- full %>%
  pivot_wider(
    names_from = organism_originalType,
    names_prefix = "organism_",
    values_from = organism_originalValue,
    values_fn = first
  ) %>%
  pivot_wider(
    names_from = reference_originalType,
    names_prefix = "reference_",
    values_from = reference_originalValue,
    values_fn = first
  ) %>%
  pivot_wider(
    names_from = structure_originalType,
    names_prefix = "structure_",
    values_from = structure_originalValue,
    values_fn = first
  ) %>%
  select(
    database,
    organism_organismClean = organism_clean,
    organism_organismDirty = organism_dirty,
    structure_structureSmiles = structure_smiles,
    structure_structureInchi = structure_inchi,
    structure_structureNominal = structure_nominal,
    reference_referenceDoi = reference_doi,
    reference_referencePmid = reference_pubmed,
    reference_referencePublishingDetails = reference_publishingDetails,
    reference_referenceTitle = reference_title,
    reference_referenceOriginal = reference_original,
    reference_referenceSplit = reference_split,
    cleanedType,
    cleanedValue,
    validation
  )

ready_2 <- ready_1 %>%
  pivot_longer(
    cols = 2:12,
    names_to = c("origin", "originalType"),
    names_sep = "_",
    values_to = "originalValue"
  ) %>%
  select(
    database,
    originalType,
    originalValue,
    cleanedType,
    cleanedValue,
    validation
  ) %>%
  distinct() %>%
  filter(!is.na(originalValue))

ready_3 <- ready_2 %>%
  filter(validation == "validated") %>%
  group_by(
    database,
    validation
  ) %>%
  count(name = "count") %>%
  arrange(desc(count))

sunk <- ready_2 %>%
  filter(database %in% ready_3$database) %>%
  group_by(
    database,
    originalType,
    cleanedType,
    validation
  ) %>%
  count(name = "count") %>%
  filter(substr(
    x = gsub(
      pattern = ".*_",
      replacement = "",
      x = cleanedType
    ),
    start = 1,
    stop = 8
  ) %in% substr(
    x = originalType,
    start = 1,
    stop = 8
  )) %>%
  ungroup() %>%
  arrange(desc(
    count,
    validation,
    database
  )) %>%
  mutate(
    validation = gsub(
      pattern = "notValidated",
      replacement = "not_validated",
      x = validation
    ),
    originalType = gsub(
      pattern = "organismClean",
      replacement = "organism_clean",
      x = originalType
    ),
    originalType = gsub(
      pattern = "organismDirty",
      replacement = "organism_dirty",
      x = originalType
    ),
    originalType = gsub(
      pattern = "structureSmiles",
      replacement = "structure_smiles",
      x = originalType
    ),
    originalType = gsub(
      pattern = "structureInchi",
      replacement = "structure_inchi",
      x = originalType
    ),
    originalType = gsub(
      pattern = "structureNominal",
      replacement = "structure_name",
      x = originalType
    ),
    originalType = gsub(
      pattern = "referenceDoi",
      replacement = "reference_doi",
      x = originalType
    ),
    originalType = gsub(
      pattern = "referencePmid",
      replacement = "reference_pmid",
      x = originalType
    ),
    originalType = gsub(
      pattern = "referencePublishingDetails",
      replacement = "reference_publishing_details",
      x = originalType
    ),
    originalType = gsub(
      pattern = "referenceOriginal",
      replacement = "reference_original",
      x = originalType
    ),
    originalType = gsub(
      pattern = "referenceSplit",
      replacement = "reference_split",
      x = originalType
    ),
    originalType = gsub(
      pattern = "referenceTitle",
      replacement = "reference_title",
      x = originalType
    ),
    cleanedType = gsub(
      pattern = "cleaned_",
      replacement = "",
      x = cleanedType
    ),
  )

legend <- with(sunk, reorder(database, count))

legend_v <- with(sunk, reorder(validation, desc(count)))

cat("drawing alluvial \n")

if (mode == "full") {
  pdf(
    file = file.path("../res", "alluvial.pdf"),
    width = 96,
    height = 54
  )

  ggplot(
    as.data.frame(sunk),
    aes(
      y = count,
      axis1 = database,
      axis2 = originalType,
      axis3 = cleanedType
    )
  ) +
    geom_stratum(
      width = 1 / 2,
      aes(size = 1),
      decreasing = TRUE
    ) +
    geom_alluvium(
      width = 1 / 2,
      aes(fill = legend_v),
      aes.bind = "alluvia",
      lode.guidance = "forward",
      decreasing = TRUE
    ) +
    geom_flow(
      width = 1 / 2,
      aes(
        fill = legend_v,
        colour = legend_v
      ),
      aes.bind = "alluvia",
      aes.flow = "forward",
      stat = after_stat("alluvium"),
      decreasing = TRUE
    ) +
    geom_fit_text(
      stat = "stratum",
      min.size = 0,
      grow = TRUE,
      width = 1 / 2,
      aes(label = after_stat(stratum)),
      decreasing = TRUE
    ) +
    scale_x_discrete(limits = c("database", "original", "cleaned")) +
    # scale_y_continuous(trans = 'log10', name = "log10(count)") +
    scale_fill_manual(values = c("#D71D62", "#08589B")) +
    scale_colour_manual(values = c("#D71D62", "#08589B")) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      # axis.text.y = element_blank(),
      axis.text.y = element_text(size = rel(7)),
      axis.title.y = element_text(size = rel(1)),
      # axis.ticks = element_blank(),
      axis.text.x = element_text(size = rel(7)),
      axis.title.x = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = rel(5)),
    )
  dev.off()
}
## interactive alternative

# library(easyalluvial)
# library(parcats)
# sunk2 <- ready_2 %>%
#   filter(database %in% ready_3$database) %>%
#   arrange(desc(validation, database)) %>%
#   mutate(
#     validation = gsub(
#       pattern = "notValidated",
#       replacement = "not_validated",
#       x = validation
#     ),
#     originalType = gsub(
#       pattern = "organismOriginal",
#       replacement = "organism_original",
#       x = originalType
#     ),
#     originalType = gsub(
#       pattern = "structureSmiles",
#       replacement = "structure_smiles",
#       x = originalType
#     ),
#     originalType = gsub(
#       pattern = "structureInchi",
#       replacement = "structure_inchi",
#       x = originalType
#     ),
#     originalType = gsub(
#       pattern = "structureNominal",
#       replacement = "structure_name",
#       x = originalType
#     ),
#     originalType = gsub(
#       pattern = "referenceDoi",
#       replacement = "reference_doi",
#       x = originalType
#     ),
#     originalType = gsub(
#       pattern = "referencePmid",
#       replacement = "reference_pmid",
#       x = originalType
#     ),
#     originalType = gsub(
#       pattern = "referencePublishingDetails",
#       replacement = "reference_publishing_details",
#       x = originalType
#     ),
#     originalType = gsub(
#       pattern = "referenceOriginal",
#       replacement = "reference_original",
#       x = originalType
#     ),
#     originalType = gsub(
#       pattern = "referenceSplit",
#       replacement = "reference_split",
#       x = originalType
#     ),
#     originalType = gsub(
#       pattern = "referenceTitle",
#       replacement = "reference_title",
#       x = originalType
#     ),
#     cleanedType = gsub(
#       pattern = "cleaned_",
#       replacement = "",
#       x = cleanedType
#     ),
#   ) %>%
#   select(-originalValue,
#          -cleanedType,
#          -cleanedValue)
#
# sunk2$validation <- as.factor(sunk2$validation)
#
# test <- sunk %>%
#   arrange(desc(count)) %>%
#   distinct(database)
#
# p <- alluvial_wide(data = sunk2,
#                    order_levels = test$database,
#                    fill_by = "last_variable")
#
# parcats(
#   p,
#   marginal_histograms = TRUE,
#   sunk2,
#   width = 1600,
#   height = 900
# )
