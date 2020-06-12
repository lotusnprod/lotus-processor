#title: "PLANTCYC compileR"

#loading
##functions
source("../../functions.R")

##db
db <- "PLANTCYC"

##paths
inpath <- "0_initial_files/0_data"

dirnames <- data.frame(list.dirs(inpath))

dirnames[, 1] <- as.character(dirnames[, 1])

dirnames <- dirnames[seq(4, nrow(dirnames), 6),]

file <- "compounds.dat"

for(i in dirnames) {
  plantcycleaning(i)
}
