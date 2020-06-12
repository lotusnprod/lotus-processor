# title: "Reference sanitizeR"

# setting working directory
setwd("~/GitLab/opennaturalproductsdb/src/")

# loading paths
source("paths.R")

# loading functions
source("functions.R")

# loading files
dataTranslated <- read_delim(
  file = gzfile(pathTranslatedReference),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

### Find something appropriate
test <- dataTranslated %>%
  filter(!is.na(translatedTitle)) %>%
  group_by(translatedTitle) %>%
  distinct(referenceSplit, .keep_all = TRUE) %>%
  add_count() %>%
  arrange(desc(n))

test2 <- test %>%
  filter(grepl("Harborne, The Handbook of Natural Flavonoids", referenceSplit)) %>%
  mutate(
    referenceSplit = gsub(
      pattern = "Harborne, The Handbook of Natural Flavonoids, 2, (1999), 1,Anthocyanins",
      replacement = "",
      x = referenceSplit,
      fixed = TRUE
    ),
    referenceSplit = gsub(
      pattern = "Harborne, The Handbook of Natural Flavonoids, 1, (1999), 2.Flavones, John Wiley & Son",
      replacement = "",
      x = referenceSplit,
      fixed = TRUE
    ),
    referenceSplit = gsub(
      pattern = "Harborne, The Handbook of Natural Flavonoids, 1, (1999), 3.Flavone O-glycosides, John Wiley & Son",
      replacement = "",
      x = referenceSplit,
      fixed = TRUE
    ),
    referenceSplit = gsub(
      pattern = "Harborne, The Handbook of Natural Flavonoids, 1, (1999), 3.Flavone O-glycosides,",
      replacement = "",
      x = referenceSplit,
      fixed = TRUE
    ),
    referenceSplit = gsub(
      pattern = "Harborne, The Handbook of Natural Flavonoids, 1, (1999), 8.Flavone O-glycosides,",
      replacement = "",
      x = referenceSplit,
      fixed = TRUE
    ),
    referenceSplit = gsub(
      pattern = "Harborne, The Handbook of Natural Flavonoids, 1, (1999), 9.Flavone O-glycosides,",
      replacement = "",
      x = referenceSplit,
      fixed = TRUE
    ),
    referenceSplit = gsub(
      pattern = "Harborne, The Handbook of Natural Flavonoids, 1, (1999), 11.Flavone O-glycosides,",
      replacement = "",
      x = referenceSplit,
      fixed = TRUE
    ),
    referenceSplit = gsub(
      pattern = "Harborne, The Handbook of Natural Flavonoids, 1, (1999), 12.Flavone O-glycosides,",
      replacement = "",
      x = referenceSplit,
      fixed = TRUE
    ),
    referenceSplit = gsub(
      pattern = "Harborne, The Handbook of Natural Flavonoids, 1, (1999), 13.Flavone O-glycosides,",
      replacement = "",
      x = referenceSplit,
      fixed = TRUE
    ),
    referenceSplit = gsub(
      pattern = "Harborne, The Handbook of Natural Flavonoids, 1, (1999), 181.Flavonols",
      replacement = "",
      x = referenceSplit,
      fixed = TRUE
    ),
    referenceSplit = gsub(
      pattern = "Harborne, The Handbook of Natural Flavonoids, 1, (1999), 297.Flavonol O-glycosides",
      replacement = "",
      x = referenceSplit,
      fixed = TRUE
    ),
    referenceSplit = gsub(
      pattern = "Harborne, The Handbook of Natural Flavonoids, 2, (1999), 355,Flavans and proanthocyanidins",
      replacement = "",
      x = referenceSplit,
      fixed = TRUE
    ),
    referenceSplit = gsub(
      pattern = "Harborne, The Handbook of Natural Flavonoids, 2, (1999), 517,Isoflavonoids and neoflavonoids",
      replacement = "",
      x = referenceSplit,
      fixed = TRUE
    ),
    referenceSplit = gsub(
      pattern = "Harborne, The Handbook of Natural Flavonoids, 1, (1999), 549,C-glycosylflavones",
      replacement = "",
      x = referenceSplit,
      fixed = TRUE
    ),
    referenceSplit = gsub(
      pattern = "Harborne, The Handbook of Natural Flavonoids, 1, (1999), 645,Biflavonyls",
      replacement = "",
      x = referenceSplit,
      fixed = TRUE
    ),
    referenceSplit = gsub(
      pattern = "Harbone, Comparative Biochemisty of the Flavonoids,(1965),39, Academic Press",
      replacement = "",
      x = referenceSplit,
      fixed = TRUE
    ),
    referenceSplit = gsub(
      pattern = "Harborne, Comparative Biochemstry of the Flavonoides,(1967),41, Academic Press",
      replacement = "",
      x = referenceSplit,
      fixed = TRUE
    ),
    referenceSplit = gsub(
      pattern = "Harborne,Comparative Biochemistry of the Flavonoids, (1967) Academic Press",
      replacement = "",
      x = referenceSplit,
      fixed = TRUE
    ),
    referenceSplit = gsub(
      pattern = "Chalcones,dihydrochalcones and aurones",
      replacement = "",
      x = referenceSplit,
      fixed = FALSE
    ),
    referenceSplit = gsub(
      pattern = "Harborne, The Handbook of Natural Flavonoids, 2, (1999)",
      replacement = "",
      x = referenceSplit,
      fixed = TRUE
    ),
    referenceSplit = gsub(
      pattern = ", [0-9]{3},",
      replacement = "",
      x = referenceSplit,
      fixed = FALSE
    ),
    referenceSplit = gsub(
      pattern = "Harborne, The Handbook of Natural Flavonoids, 2, (1999)Flavanones and dihydroflavonols",
      replacement = "",
      x = referenceSplit,
      fixed = TRUE
    )
  )



test3 <- test %>%
  filter(n == 631) %>%
  mutate(
    referenceSplit = gsub(
      pattern = "; Khimiya .*",
      replacement = "",
      x = referenceSplit,
      fixed = FALSE
    )
  )


test4 <- test %>%
  filter(n != 631 &
           !grepl("Harborne, The Handbook of Natural Flavonoids", referenceSplit) &
           n > 10)

test5 <- test4 %>%
  filter(n == 274) # problematic issue from PlantaMed

test6 <- test4 %>%
  filter(n == 194) # problematic too... DNP

test6 <- test4 %>%
  filter(n == 194)

RefShouldBeOk <- test %>%
  filter(n < 5 & translationScore > 80 & translationScore <= 100)

### analyzing how low we can go
low70 <- test %>%
  filter(n < 5 & translationScore > 70 & translationScore <= 80)

## example
### Luesch, Hendrik; Yoshida, Wesley Y.; Moore, Richard E.; Paul, Valerie J.; Journal of Natural Products; vol. 63; 10; (2000); p. 1437 - 1439.
#### returning erroneous result
##### Journal of Natural Products; vol. 63; 10; (2000); p. 1437 - 1439
###### test <- cr_works(query = "Journal of Natural Products; vol. 63; 10; (2000); p. 1437 - 1439") #returning right result (but very low score)


### score above 100 seems to indicate multiple references in 1 row

dataSanitized <- dataTranslated %>%
  mutate(
    sanitizedTitle = translatedTitle,
    sanitizedJournal = translatedJournal,
    sanitizedDoi = translatedDoi,
    sanitizedAuthor = translatedAuthor,
    sanitizedDate = translatedDate,
    sanitizedTranslationScore = translationScore
  )

# export
write.table(
  x = dataSanitized,
  file = gzfile(
    description = pathSanitizedReference,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
