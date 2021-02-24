cat("This script draws the violin plots \n")

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")

cat("... libraries \n")
library(data.table)
library(plotly)

organismsPerStructure_3D <-
  fread(file = file.path(
    pathDataProcessed,
    "organismsPerStructure.tsv.gz"
  ))

structuresPerOrganism_3D <-
  fread(file = file.path(
    pathDataProcessed,
    "structuresPerOrganism.tsv.gz"
  ))
# needs 3_metrics
fig1 <- plot_ly(type = "violin") %>%
  add_trace(
    x = 0,
    y = organismsPerStructure_3D$n,
    legendgroup = "organisms per structure",
    scalegroup = "organisms per structure",
    name = "organisms per structure",
    side = "negative",
    box = list(visible = TRUE),
    meanline = list(visible = TRUE),
    color = I("#a6cee3")
  ) %>%
  add_trace(
    x = 0,
    y = structuresPerOrganism_3D$n,
    legendgroup = "structures per organism",
    scalegroup = "structures per organism",
    name = "structures per organism",
    side = "positive",
    box = list(visible = TRUE),
    meanline = list(visible = TRUE),
    color = I("#1f78b4")
  ) %>%
  layout(
    yaxis = list(
      title = "",
      zeroline = FALSE
    ),
    xaxis = list(
      title = "",
      zeroline = FALSE,
      showline = FALSE,
      showticklabels = FALSE,
      showgrid = FALSE
    )
  )
fig1

orca(
  p = fig1,
  file = file.path(pathDataProcessedFigures, "violin.pdf")
)

fig2 <- plot_ly(type = "violin") %>%
  add_trace(
    x = 0,
    y = organismsPerStructure_3D$n,
    legendgroup = "organisms per structure",
    scalegroup = "organisms per structure",
    name = "organisms per structure",
    side = "negative",
    box = list(visible = TRUE),
    meanline = list(visible = TRUE),
    color = I("#a6cee3")
  ) %>%
  add_trace(
    x = 0,
    y = structuresPerOrganism_3D$n,
    legendgroup = "structures per organism",
    scalegroup = "structures per organism",
    name = "structures per organism",
    side = "positive",
    box = list(visible = TRUE),
    meanline = list(visible = TRUE),
    color = I("#1f78b4")
  ) %>%
  layout(
    yaxis = list(
      range = c(0, 200),
      title = "",
      zeroline = FALSE
    ),
    xaxis = list(
      title = "",
      zeroline = FALSE,
      showline = FALSE,
      showticklabels = FALSE,
      showgrid = FALSE
    )
  )
fig2

orca(
  p = fig2,
  file = file.path(pathDataProcessedFigures, "violin_zoomed.pdf")
)

htmlwidgets::saveWidget(
  widget = as_widget(fig2),
  file = file.path(pathDataProcessedFiguresHtml, "violin.html")
)

end <- Sys.time()

cat("Script finished in", format(end - start), "\n")
