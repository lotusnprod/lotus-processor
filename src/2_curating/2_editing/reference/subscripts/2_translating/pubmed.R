# title: "Ref translatoR"

# loading paths
source("paths.R")

# loading functions
source("functions/reference.R")

# getting references ##getting them with pubmed API and not crossRef because crossRef pubmed ID not working!!
## 2
# mc cores set to 2 because fails otherwise (entrez limitation probably)
reflistPubmed <- invisible(
  pbmclapply(
    FUN = getrefPubmed,
    X = as.character(dataReferenceFillPubmed$value),
    mc.preschedule = TRUE,
    mc.set.seed = TRUE,
    mc.silent = TRUE,
    mc.cores = 2,
    mc.cleanup = TRUE,
    mc.allow.recursive = TRUE, 
    ignore.interactive = TRUE
  )
)

reflistPubmedBound <- bind_rows(reflistPubmed)

# joining with original dataframe
for (i in 1:nrow(reflistPubmedBound)) {
  dataReferenceFillPubmed[i, "translatedDoi"] <-
    reflistPubmedBound[i, "translatedDoi"]
}

for (i in 1:nrow(reflistPubmedBound)) {
  dataReferenceFillPubmed[i, "translatedJournal"] <-
    reflistPubmedBound[i, "translatedJournal"]
}

for (i in 1:nrow(reflistPubmedBound)) {
  dataReferenceFillPubmed[i, "translatedTitle"] <-
    reflistPubmedBound[i, "translatedTitle"]
}

for (i in 1:nrow(reflistPubmedBound)) {
  dataReferenceFillPubmed[i, "translatedAuthor"] <-
    reflistPubmedBound[i, "translatedAuthor"]
}

for (i in 1:nrow(reflistPubmedBound)) {
  dataReferenceFillPubmed[i, "translatedDate"] <-
    reflistPubmedBound[i, "translatedDate"]
}

for (i in 1:nrow(reflistPubmedBound)) {
  dataReferenceFillPubmed[i, "translationScore"] <- 1
}
