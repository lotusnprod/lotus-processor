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

orca(p = fig1, file = file.path(pathDataProcessedFigures, "violin.pdf"))

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

orca(p = fig2, file = file.path(pathDataProcessedFigures, "violin_zoomed.pdf"))
