# title: "Ref translatoR"

# loading paths
source("paths.R")

# loading functions
source("functions.R")

# loading files
## reference
dataReference <- read_delim(
  file = gzfile(pathOriginalRef),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

# splitting when MULTIPLE REFERENCES (|)
dataReferenceSplit <- dataReference %>%
  mutate(newreference = referenceOriginal) %>%
  cSplit(splitCols = "newreference",
         sep = "|") %>%
  mutate_all(as.character) %>%
  tibble()

# pivoting
dataReferenceLong <- dataReferenceSplit %>%
  rowwise() %>%
  pivot_longer(cols = 2:ncol(.)) %>%
  filter(!is.na(value)) %>%
  ungroup()

dataReferenceLong$value <- y_as_na(x = dataReferenceLong$value,
                                   y = "")

dataReferenceLong$value <- y_as_na(x = dataReferenceLong$value,
                                   y = "NA")

dataReferenceLong <- dataReferenceLong %>%
  filter(!is.na(value))

# splitting when multiple fields FOR A SINGLE REFERENCE (§)
dataReferenceLongSplit <- dataReferenceLong %>%
  cSplit(splitCols = "value",
         sep = "§") %>%
  rowwise() %>%
  pivot_longer(
    cols = (ncol(.) - 1):ncol(.),
    names_to = c(".value", "level"),
    names_sep = "_",
    values_to = "taxonomy",
    values_drop_na = TRUE
  ) %>%
  ungroup()

# selecting fields probably corresponding to a title or authors list (human)
dataReferenceFillAuto <- dataReferenceLongSplit %>%
  filter(
    !grepl(pattern = "^KNApSAcK Database$",
           x = value) &
      !grepl(pattern = "^DrDuke$",
             x = value) &
      !grepl(pattern = "^patent:",
             x = value) &
      !grepl(pattern = "^CHEBI:",
             x = value) &
      !grepl(pattern = "^doi:",
             x = value) &
      !grepl(pattern = "^10\\.\\d{4,9}",
             x = value)  &
      !grepl(pattern = "^pubmed:",
             x = value) &
      !grepl(pattern = "^PUBMED ",
             x = value) &
      !grepl(pattern = "^[0-9]*$",
             x = value) &
      !grepl(pattern = "^\\d+(;\\d+)*$",
             x = value)  &
      !grepl(pattern = "http://www.ncbi.nlm.nih.gov/pubmed/",
             x = value) &
      str_count(string = value) > 28
  )

dataReferenceFillManual <- dataReferenceLongSplit %>%
  filter(str_count(string = value) <= 28 &
           grepl(pattern = "et al",
                 x = value))

# getting references
## 1
reflist <- invisible(
  pbmclapply(
    FUN = getref,
    X = dataReferenceFillAuto$value,
    mc.preschedule = TRUE,
    mc.set.seed = TRUE,
    mc.silent = TRUE,
    mc.cores = (parallel::detectCores() - 2),
    mc.cleanup = TRUE,
    mc.allow.recursive = TRUE
  )
)

# joining with original dataframe
for (i in 1:length(reflist)) {
  dataReferenceFillAuto[i, "translatedDoi"] <-
    as.character(ifelse(
      test = !is.na(reflist[[i]]),
      yes = ifelse(
        test = !is.null(reflist[[i]][["data"]][["doi"]]),
        yes = reflist[[i]][["data"]][["doi"]],
        no = NA
      ),
      no = NA
    )[1])
  
  dataReferenceFillAuto[i, "translatedJournal"] <-
    as.character(ifelse(
      test = !is.na(reflist[[i]]),
      yes = ifelse(
        test = !is.null(reflist[[i]][["data"]][["container.title"]]),
        yes = reflist[[i]][["data"]][["container.title"]],
        no = NA
      ),
      no = NA
    )[1])
  
  dataReferenceFillAuto[i, "translatedTitle"] <-
    as.character(ifelse(
      test = !is.na(reflist[[i]]),
      yes = ifelse(
        test = !is.null(reflist[[i]][["data"]][["title"]]),
        yes = reflist[[i]][["data"]][["title"]],
        no = NA
      ),
      no = NA
    )[1])
  
  dataReferenceFillAuto[i, "translatedDate"] <-
    as.character(ifelse(
      test = !is.na(reflist[[i]]),
      yes = ifelse(
        test = !is.null(reflist[[i]][["data"]][["issued"]]),
        yes = reflist[[i]][["data"]][["issued"]],
        no = NA
      ),
      no = NA
    )[1])
  
  dataReferenceFillAuto[i, "translatedAuthor"] <-
    as.character(ifelse(
      test = !is.na(reflist[[i]][["data"]][["author"]][[1]][["family"]][1]),
      yes = ifelse(
        test = !is.null(reflist[[i]][["data"]][["author"]][[1]][["family"]][1]),
        yes = reflist[[i]][["data"]][["author"]][[1]][["family"]][1],
        no = NA
      ),
      no = NA
    )[1])
  
  dataReferenceFillAuto[i, "translationScore"] <-
    as.character(ifelse(
      test = !is.na(reflist[[i]][["data"]][["score"]]),
      yes = ifelse(
        test = !is.null(reflist[[i]][["data"]][["score"]]),
        yes = reflist[[i]][["data"]][["score"]],
        no = 0
      ),
      no = 0
    )[1])
}

# selecting fields probably corresponding to a title or authors list (human)
dataReferenceFillPubmed <- dataReferenceLongSplit %>%
  filter(
    grepl(pattern = "^pubmed:",
          x = value) |
      grepl(pattern = "^PUBMED ",
            x = value) |
      grepl(pattern = "^[0-9]*$",
            x = value) |
      grepl(pattern = "^\\d+(;\\d+)*$",
            x = value) |
      grepl(pattern = "http://www.ncbi.nlm.nih.gov/pubmed/",
            x = value)
  ) %>%
  cSplit(splitCols = "value",
         sep = ";") %>%
  select(-level) %>%
  mutate_all(as.character) %>%
  group_by(referenceOriginal) %>%
  pivot_longer(
    cols = 3:ncol(.),
    names_to = c(".value", "level"),
    names_sep = "_",
    values_to = "taxonomy",
    values_drop_na = TRUE
  ) %>%
  ungroup() %>%
  mutate(value = gsub(
    pattern = "pubmed:",
    replacement = "",
    x = value
  )) %>%
  mutate(value = gsub(
    pattern = "PUBMED ",
    replacement = "",
    x = value
  )) %>%
  mutate(value = gsub(
    pattern = "http://www.ncbi.nlm.nih.gov/pubmed/",
    replacement = "",
    x = value
  )) %>%
  mutate(value = gsub(
    pattern = "?term=",
    replacement = "",
    x = value
  ))

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

# selecting fields corresponding to other sources than articles
dataReferenceFillNoArticle <- dataReferenceLongSplit %>%
  filter(
    grepl(pattern = "^KNApSAcK Database$",
          x = value) |
      grepl(pattern = "^DrDuke$",
            x = value) |
      grepl(pattern = "^patent:",
            x = value) |
      grepl(pattern = "^CHEBI:",
            x = value)
  ) %>%
  mutate(
    translatedDoi = NA,
    translatedJournal = NA,
    translatedTitle = NA,
    translatedDate = NA,
    translatedAuthor = NA,
    translationScore = 0
  )

# selecting fields corresponding to DOIs
dataReferenceFillDoi <- dataReferenceLongSplit %>%
  filter(grepl("^doi:", value) |
           grepl("^10\\.\\d{4,9}", value))

# getting references
## 3
reflistDoi <- invisible(
  pbmclapply(
    FUN = getrefDoi,
    X = dataReferenceFillDoi$value,
    mc.preschedule = TRUE,
    mc.set.seed = TRUE,
    mc.silent = TRUE,
    mc.cores = (parallel::detectCores() - 2),
    mc.cleanup = TRUE,
    mc.allow.recursive = TRUE
  )
)

# joining with original dataframe
for (i in 1:length(reflistDoi)) {
  dataReferenceFillDoi[i, "translatedDoi"] <-
    as.character(ifelse(
      test = !is.na(reflistDoi[[i]]),
      yes = ifelse(
        test = !is.null(reflistDoi[[i]][["data"]][["doi"]]),
        yes = reflistDoi[[i]][["data"]][["doi"]],
        no = NA
      ),
      no = NA
    )[1])
  
  dataReferenceFillDoi[i, "translatedJournal"] <-
    as.character(ifelse(
      test = !is.na(reflistDoi[[i]]),
      yes = ifelse(
        test = !is.null(reflistDoi[[i]][["data"]][["container.title"]]),
        yes = reflistDoi[[i]][["data"]][["container.title"]],
        no = NA
      ),
      no = NA
    )[1])
  
  dataReferenceFillDoi[i, "translatedTitle"] <-
    as.character(ifelse(
      test = !is.na(reflistDoi[[i]]),
      yes = ifelse(
        test = !is.null(reflistDoi[[i]][["data"]][["title"]]),
        yes = reflistDoi[[i]][["data"]][["title"]],
        no = NA
      ),
      no = NA
    )[1])
  
  dataReferenceFillDoi[i, "translatedDate"] <-
    as.character(ifelse(
      test = !is.na(reflistDoi[[i]]),
      yes = ifelse(
        test = !is.null(reflistDoi[[i]][["data"]][["issued"]]),
        yes = reflistDoi[[i]][["data"]][["issued"]],
        no = NA
      ),
      no = NA
    )[1])
  
  dataReferenceFillDoi[i, "translatedAuthor"] <-
    as.character(ifelse(
      test = !is.na(reflistDoi[[i]]),
      yes = ifelse(
        test = !is.null(reflistDoi[[i]][["data"]][["author"]][[1]][["family"]][1]),
        yes = reflistDoi[[i]][["data"]][["author"]][[1]][["family"]][1],
        no = NA
      ),
      no = NA
    )[1])
  
  dataReferenceFillDoi[i, "translationScore"] <-
    as.character(ifelse(
      test = !is.na(reflistDoi[[i]][["data"]][["score"]]),
      yes = ifelse(
        test = !is.null(reflistDoi[[i]][["data"]][["score"]]),
        yes = reflistDoi[[i]][["data"]][["score"]],
        no = 0
      ),
      no = 0
    )[1])
}

# temporary
dataReferenceFillManual <- dataReferenceFillManual %>%
  mutate(
    translatedDoi = NA,
    translatedJournal = NA,
    translatedTitle = NA,
    translatedDate = NA,
    translatedAuthor = NA,
    translationScore = 0
  )

# joining all types together again
dataFull <- rbind(
  dataReferenceFillAuto,
  dataReferenceFillDoi,
  dataReferenceFillManual,
  dataReferenceFillNoArticle,
  dataReferenceFillPubmed
)

# joining with original df
dataReferenceMin <- dataReferenceLongSplit %>%
  select(referenceOriginal)

dataReferenced <- left_join(dataReferenceMin, dataFull)

### problematic reference field containing multiple references with no clear association ###
### I also think we should output a "nonarticleref" column for external DB links etc ###
### inconsistency of journal name depending on retrieval method (check JNP) ###

dataReferencedSelected <- dataReferenced %>%
  select(
    referenceOriginal,
    referenceSplit = value,
    translatedDoi,
    translatedJournal,
    translatedTitle,
    translatedDate,
    translatedAuthor,
    translationScore
  ) %>%
  mutate(translationScore = replace_na(translationScore, "0")) %>%
  distinct(
    referenceOriginal,
    referenceSplit,
    translatedTitle,
    translatedJournal,
    translatedDate,
    translatedAuthor,
    translatedDoi,
    translationScore,
    .keep_all = TRUE
  )

dataReferencedSelected$translationScore <-
  as.numeric(dataReferencedSelected$translationScore)

dataReferencedSelected$translationScore[dataReferencedSelected$translationScore == 1] <-
  100

# exporting
## ref
write.table(
  x = dataReferencedSelected,
  file = gzfile(
    description = pathTranslatedReference,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
