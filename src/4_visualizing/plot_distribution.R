source("r/log_debug.R")
log_debug("This script draws the distribution of organisms per structure and vice versa")

start <- Sys.time()

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")

log_debug("... libraries")
library(dplyr)
library(plotly)
library(readr)

organismsPerStructure_2D <-
  read_delim(file = file.path(
    pathDataProcessed,
    "organismsPerStructure2D.tsv.gz"
  ))

structuresPerOrganism_2D <-
  read_delim(file = file.path(
    pathDataProcessed,
    "structures2DPerOrganism.tsv.gz"
  ))

organismsPerStructure_2D_cum <- organismsPerStructure_2D %>%
  ungroup() %>%
  arrange(n) %>%
  mutate(cum_sum = cumsum(n)) %>%
  mutate(
    per_x = n / max(n) * 100,
    per_y = cum_sum / max(cum_sum) * 100
  ) %>%
  distinct(per_x, .keep_all = TRUE)

structuresPerOrganism_2D_cum <- structuresPerOrganism_2D %>%
  ungroup() %>%
  arrange(n) %>%
  mutate(cum_sum = cumsum(n)) %>%
  mutate(
    per_x = n / max(n) * 100,
    per_y = cum_sum / max(cum_sum) * 100
  ) %>%
  distinct(per_x, .keep_all = TRUE)

top_structure <-
  organismsPerStructure_2D_cum[which.max(organismsPerStructure_2D_cum$n), ]

## actually third but better to ilustrate
top_organism <-
  structuresPerOrganism_2D_cum[structuresPerOrganism_2D_cum$organismCleaned == "Arabidopsis thaliana", ]

a <- list(
  x = top_organism$per_x,
  y = top_organism$per_y,
  text = paste("<i>", top_organism$organismCleaned, "</i> \n N =", top_organism$n),
  xref = "x",
  yref = "y",
  showarrow = TRUE,
  arrowhead = 1,
  ax = -200,
  ay = 40,
  font = list(color = "#2994D2"),
  arrowcolor = "#2994D2"
)

b <- list(
  x = top_structure$per_x,
  y = top_structure$per_y,
  text = paste(
    top_structure$structureCleaned_inchikey2D,
    "\n N =",
    top_structure$n
  ),
  xref = "x",
  yref = "y",
  showarrow = TRUE,
  arrowhead = 1,
  ax = -100,
  ay = 40,
  font = list(color = "#7CB13F"),
  arrowcolor = "#7CB13F"
)

c <- list(
  x = "per organism",
  y = mean(structuresPerOrganism_2D$n),
  text = paste("mean ", round(mean(
    structuresPerOrganism_2D$n
  ), 2)),
  xref = "x2",
  yref = "y2",
  showarrow = TRUE,
  arrowhead = 1,
  ax = 57,
  ay = -40,
  font = list(color = "#2994D2"),
  arrowcolor = "#2994D2"
)

d <- list(
  x = "per structure",
  y = mean(organismsPerStructure_2D$n),
  text = paste("mean ", round(mean(
    organismsPerStructure_2D$n
  ), 2)),
  xref = "x2",
  yref = "y2",
  showarrow = TRUE,
  arrowhead = 1,
  ax = -57,
  ay = -40,
  font = list(color = "#7CB13F"),
  arrowcolor = "#7CB13F"
)

e <- list(
  x = "per structure",
  y = median(organismsPerStructure_2D$n),
  text = paste("median ", median(organismsPerStructure_2D$n)),
  xref = "x2",
  yref = "y2",
  showarrow = TRUE,
  arrowhead = 1,
  ax = -90,
  ay = -5,
  font = list(color = "#7CB13F"),
  arrowcolor = "#7CB13F"
)

f <- list(
  x = "per organism",
  y = median(structuresPerOrganism_2D$n),
  text = paste("median ", median(structuresPerOrganism_2D$n)),
  xref = "x2",
  yref = "y2",
  showarrow = TRUE,
  arrowhead = 1,
  ax = -100,
  ay = 0,
  font = list(color = "#2994D2"),
  arrowcolor = "#2994D2"
)

# initialize plot
fig <- plot_ly() %>%
  add_lines(
    data = organismsPerStructure_2D_cum,
    name = "organisms per structure",
    x = ~per_x,
    y = ~per_y,
    mode = "line",
    color = I("#7CB13F"),
    showlegend = FALSE
  ) %>%
  add_lines(
    data = structuresPerOrganism_2D_cum,
    name = "structures per organism",
    x = ~per_x,
    y = ~per_y,
    mode = "line",
    color = I("#2994D2"),
    showlegend = FALSE
  ) %>%
  add_trace(
    data = structuresPerOrganism_2D,
    name = "structures per organism",
    x = "per organism",
    y = ~n,
    type = "box",
    quartilemethod = "inclusive",
    boxmean = TRUE,
    xaxis = "x2",
    yaxis = "y2",
    color = I("#2994D2")
  ) %>%
  add_trace(
    data = organismsPerStructure_2D,
    name = "organisms per structure",
    x = "per structure",
    y = ~n,
    type = "box",
    quartilemethod = "inclusive",
    boxmean = TRUE,
    xaxis = "x2",
    yaxis = "y2",
    color = I("#7CB13F")
  ) %>%
  layout(
    annotations = list(a, b, c, d, e, f),
    font = list(family = "helvetica neue", size = 18),
    yaxis = list(title = "Cumulative % of maximal contribution"),
    xaxis = list(title = "% of individuals"),
    yaxis2 = list(
      title = "Number of individuals",
      mirror = TRUE,
      showline = TRUE,
      zeroline = FALSE,
      showline = FALSE,
      showgrid = FALSE,
      range = c(0, 25),
      domain = c(0.3, 0.7),
      anchor = "x2"
    ),
    xaxis2 = list(
      mirror = TRUE,
      showline = TRUE,
      zeroline = FALSE,
      showline = FALSE,
      showticklabels = FALSE,
      showgrid = FALSE,
      domain = c(0.4, 0.99),
      anchor = "y2"
    ),
    legend = list(
      x = 0.4,
      y = 0.25,
      orientation = "h"
    )
  )

fig

if (mode == "full") {
  orca(
    p = fig,
    file = file.path("../res", "distribution.pdf"),
    width = 800,
    height = 450
  )

  setwd("../res/html")
  htmlwidgets::saveWidget(
    widget = as_widget(fig),
    file = "distribution.html"
  )
}

end <- Sys.time()

log_debug("Script finished in", format(end - start))
