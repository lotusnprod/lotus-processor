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
    "organismsPerStructure2D.tsv.gz"
  ))

structuresPerOrganism_3D <-
  fread(file = file.path(
    pathDataProcessed,
    "structures2DPerOrganism.tsv.gz"
  ))

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
    color = I("#2994D2")
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
    color = I("#08589B")
  ) %>%
  layout(
    yaxis = list(
      title = "Count",
      zeroline = FALSE
    ),
    xaxis = list(
      title = "Density",
      zeroline = FALSE,
      showline = FALSE,
      showticklabels = FALSE,
      showgrid = FALSE
    )
  )
fig1

orca(
  p = fig1,
  file = file.path("../res", "violin.pdf")
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
    color = I("#2994D2")
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
    color = I("#08589B")
  ) %>%
  layout(
    yaxis = list(
      range = c(0, 200),
      title = "Count",
      zeroline = FALSE
    ),
    xaxis = list(
      title = "Density",
      zeroline = FALSE,
      showline = FALSE,
      showticklabels = FALSE,
      showgrid = FALSE
    )
  )
fig2

orca(
  p = fig2,
  file = file.path("../res", "violin_zoomed.pdf")
)

setwd("../res/html")
htmlwidgets::saveWidget(
  widget = as_widget(fig2),
  file = "violin.html"
)

end <- Sys.time()

cat("Script finished in", format(end - start), "\n")
