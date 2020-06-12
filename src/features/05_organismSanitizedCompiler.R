# title: "Organisms (sanitized) compileR"

# setting working directory
setwd("~/GitLab/opennaturalproductsdb/src/")

# loading paths
source("paths.R")

# loading functions
source("functions.R")

#writing path
dir <- pathTranslatedOrganismDistinct

length <- length(list.files(path = dir,
                            pattern = 'tsv'))

cut <- 10000

num <- as.integer(seq(
  from = 1 * cut,
  to = length * cut,
  by = cut
))

# cleaning GNFinder output
for (i in num) {
  gnfinder_cleaning(num = i)
}
