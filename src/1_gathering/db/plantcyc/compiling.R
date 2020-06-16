#title: "PLANTCYC compileR"

# setting working directory
setwd("~/GitLab/opennaturalproductsdb/src/")

# loading paths
source("paths.R")

# loading functions
source("functions.R")

## paths
dirnames <-
  data.frame(list.dirs(pathDataExternalDbSourcePlantcycDir))

dirnames[, 1] <- as.character(dirnames[, 1])

dirnames <- dirnames[seq(4, nrow(dirnames), 6), ]

file <- "compounds.dat"

for (i in dirnames) {
  plantcycompiling(x = i)
}
