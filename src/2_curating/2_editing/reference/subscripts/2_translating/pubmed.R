# title: "Ref translatoR"

# loading paths
source("paths.R")

# loading functions
source("functions.R")

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
    mc.allow.recursive = TRUE
  )
)

reflistPubmedBound <- bind_rows(reflistPubmed)

# joining with original dataframe
for (i in 1:nrow(reflistPubmedBound)) {
  dataReferenceFillPubmed[i, "translatedDoi"] <-
    reflistPubmedBound[i, "translatedDoi"]
  
  dataReferenceFillPubmed[i, "translatedJournal"] <-
    reflistPubmedBound[i, "translatedJournal"]
  
  dataReferenceFillPubmed[i, "translatedTitle"] <-
    reflistPubmedBound[i, "translatedTitle"]
  
  dataReferenceFillPubmed[i, "translatedAuthor"] <-
    reflistPubmedBound[i, "translatedAuthor"]
  
  dataReferenceFillPubmed[i, "translatedDate"] <-
    reflistPubmedBound[i, "translatedDate"]
  
  dataReferenceFillPubmed[i, "translationScore"] <- 1
}
