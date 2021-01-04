# additional fig

# needs 3_metrics
fig1 <- plot_ly(type = "violin") %>%
  add_trace(
    x = 1,
    y = organismsPerStructure$n,
    legendgroup = "organisms per structure",
    scalegroup = "organisms per structure",
    name = "organisms per structure",
    side = "negative",
    box = list(visible = T),
    meanline = list(visible = T),
    color = I("#a6cee3")
  ) %>%
  add_trace(
    x = 1,
    y = structuresPerOrganism$n,
    legendgroup = "structures per organism",
    scalegroup = "structures per organism",
    name = "structures per organism",
    side = "positive",
    box = list(visible = T),
    meanline = list(visible = T),
    color = I("#1f78b4")
  ) %>%
  layout(yaxis = list(
    title = "",
    zeroline = F
  ))
fig1

orca(
  p = fig1,
  file = file.path(pathDataProcessedFigures, "violin.pdf")
)

fig2 <- plot_ly(type = "violin") %>%
  add_trace(
    x = 1,
    y = organismsPerStructure$n,
    legendgroup = "organisms per structure",
    scalegroup = "organisms per structure",
    name = "organisms per structure",
    side = "negative",
    box = list(visible = T),
    meanline = list(visible = T),
    color = I("#a6cee3")
  ) %>%
  add_trace(
    x = 1,
    y = structuresPerOrganism$n,
    legendgroup = "structures per organism",
    scalegroup = "structures per organism",
    name = "structures per organism",
    side = "positive",
    box = list(visible = T),
    meanline = list(visible = T),
    color = I("#1f78b4")
  ) %>%
  layout(yaxis = list(
    range = c(0, 100),
    title = "",
    zeroline = F
  ))
fig2

orca(
  p = fig2,
  file = file.path(pathDataProcessedFigures, "violin_zoomed.pdf")
)

## alluvial alternative

# library(easyalluvial)
# library(parcats)
#
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
# # sunk2$database <- as.factor(sunk2$database)
# # sunk2$originalType <- as.factor(sunk2$originalType)
# # sunk2$cleanedType <- as.factor(sunk2$cleanedType)
# # sunk2$cleanedValue <- as.factor(sunk2$cleanedValue)
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
