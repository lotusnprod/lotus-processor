# title: "treating bio"

# loading paths
source("paths.R")

# loading functions
#source("functions/bio.R")

log_debug("  Step 3")

system(command = paste("bash", pathTranslatedGnfinderScript))

length <- length(list.files(path = pathTranslatedOrganismDistinct,
                            pattern = 'tsv'))

num <- as.integer(seq(
  from = 1 * cut,
  to = length * cut,
  by = cut
))

dataCleanTranslatedOrganism <- list()

# cleaning GNFinder output
for (i in num) {
  j <- i / cut
  tryCatch({
    dataCleanTranslatedOrganism[[j]] <-
      gnfinder_cleaning(num = i,
                        organismCol = "organismInterim")
  }, error = function(e) {
    cat("ERROR :", conditionMessage(e), "\n")
  })
}

# selecting and reordering
dataCleanedTranslatedOrganism <-
  bind_rows(dataCleanTranslatedOrganism) %>%
  select(
    organismInterim,
    organismCleaned = canonicalname,
    organismDbTaxo = dbTaxo,
    everything()
  ) %>%
  select(-nchar, -sum)

dataCleanedTranslatedOrganism2join <-
  dataInterimOrganismToFill %>%
  select(organismOriginal, organismInterim) %>%
  mutate_all(as.character)

dataCleanedTranslatedOrganismFull <-
  left_join(dataCleanedTranslatedOrganism2join,
            dataCleanedTranslatedOrganism) %>%
  select(-organismInterim)

dataCleanedOrganism <-
  rbind(dataCleanedOriginalOrganism,
        dataCleanedTranslatedOrganismFull)

dataCleanedOrganism <- dataCleanedOrganism %>%
  distinct(organismOriginal,
           organismCleaned,
           .keep_all = TRUE) %>%
  group_by(organismOriginal) %>%
  add_count() %>%
  ungroup() %>%
  filter(!is.na(organismCleaned) |
           !n > 1) %>%
  select(-n)

dataCleanedOrganismManipulated <-
  manipulating_taxo(dfsel = dataCleanedOrganism,
                    dic = taxaRanksDictionary) %>%
  select(-rank, -taxonomy)
