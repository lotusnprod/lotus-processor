# title: "Phytohub cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

groundhog.library(splitstackshape, date = groundhog.day)
groundhog.library(tidyverse, date = groundhog.day)

# get paths
database <- databases$get("phytohub")

## files
data_original <- read_delim(
  file = gzfile(database$sourceFiles$tsv),
  delim = "\t",
  quote = ""
) %>%
  mutate_all(as.character)

data_manipulated <- data_original %>%
  mutate(
    reference = ifelse(
      test = reference == "Belitz, H.-D., Grosch, W., Schieberle, P. (2009). Food Chemistry. Springer. [ISBN:3540408185  ] " |
        reference == "21. (2009). In Food Chemistry (pp. 938). Springer. [ISBN:3540408185  ] " |
        reference == "21. (2009). In Food Chemistry (pp. 943). Springer. [ISBN:3540408185  ] FooDB [Link] " |
        reference == "Speer, K., Hruschka, A., Kurzrock, T. and Kölling-Speer, I. (2000). 25. In Caffeinated Beverages. Health benefits, physiological effects, and chemistry (pp. 241-251). ACS Symposium series 754. [ISBN:0-841-3654-2  ] " |
        reference == "Speer, K., Hruschka, A., Kurzrock, T. and Kölling-Speer, I.  (2000). 25. In Caffeinated Beverages. Health benefits, physiological effects, and chemistry (pp. 241-251). ACS Symposium series 754. [ISBN:0-841-3654-2  ] ",
      yes = reference,
      no = sub(
        pattern = ":",
        replacement = "§",
        x = reference
      )
    )
  ) %>%
  cSplit("reference",
    sep = "§",
    fixed = TRUE,
    drop = FALSE
  ) %>%
  cSplit("reference_1",
    sep = "[ISBN:",
    fixed = TRUE,
    stripWhite = FALSE
  ) %>%
  cSplit("reference_2",
    sep = "[ISBN:",
    fixed = TRUE,
    stripWhite = FALSE
  ) %>%
  mutate(reference_2_1 = sub(
    pattern = "doi: ",
    replacement = "doi§ ",
    x = reference_2_1
  )) %>%
  cSplit(
    "reference_2_1",
    sep = "doi§ ",
    fixed = TRUE,
    stripWhite = FALSE,
  ) %>%
  mutate(
    reference_2_1_1 = sub(
      pattern = "[PubMed:",
      replacement = "[PubMed§",
      fixed = TRUE,
      x = reference_2_1_1
    ),
    reference_2_1_2 = sub(
      pattern = "[PubMed:",
      replacement = "[PubMed§",
      fixed = TRUE,
      x = reference_2_1_2
    )
  ) %>%
  cSplit(
    "reference_2_1_1",
    sep = "[PubMed§",
    fixed = TRUE,
    stripWhite = FALSE
  ) %>%
  cSplit(
    "reference_2_1_2",
    sep = "[PubMed§",
    fixed = TRUE,
    stripWhite = FALSE
  ) %>%
  mutate_all(as.character) %>%
  mutate(
    referenceIsbn_1 = gsub(
      pattern = "  ]",
      replacement = "",
      x = reference_1_2
    ),
    referenceIsbn_2 = gsub(
      pattern = "  ]",
      replacement = "",
      x = reference_2_2
    )
  ) %>%
  mutate(
    referenceAuthors_1 = ifelse(
      test = !is.na(referenceIsbn_1),
      yes = NA,
      no = reference_1_1
    ),
    referenceUnsplittable_1 = reference_2_1_1_1,
    referenceDoi_1_temp = reference_2_1_2_1,
    referencePubmed_1_temp = ifelse(
      test = !is.na(reference_2_1_2_2),
      yes = reference_2_1_2_2,
      no = reference_2_1_1_2
    )
  ) %>%
  mutate(
    referenceDoi_1_temp = sub(
      pattern = " ",
      replacement = "§",
      x = referenceDoi_1_temp
    ),
    referencePubmed_1_temp = sub(
      pattern = "  ]",
      replacement = "§",
      x = referencePubmed_1_temp
    )
  ) %>%
  cSplit("referenceDoi_1_temp", sep = "§") %>%
  cSplit("referencePubmed_1_temp", sep = "§") %>%
  mutate_all(as.character) %>%
  select(
    name,
    inchi,
    smiles,
    biologicalsource,
    reference,
    referenceIsbn_1,
    referenceAuthors_1,
    referenceUnsplittable_1,
    referenceDoi_1 = referenceDoi_1_temp_1,
    referencePubmed_1 = referencePubmed_1_temp_1,
    referenceNew = referencePubmed_1_temp_2,
    referenceIsbn_2
  ) %>%
  mutate(
    referenceNew = ifelse(
      test = referenceNew == "Belitz, H.-D., Grosch, W., Schieberle, P. (2009). Food Chemistry. Springer. [ISBN:3540408185  ] " |
        referenceNew == "21. (2009). In Food Chemistry (pp. 938). Springer. [ISBN:3540408185  ] " |
        referenceNew == "21. (2009). In Food Chemistry (pp. 943). Springer. [ISBN:3540408185  ] FooDB [Link] " |
        referenceNew == "Speer, K., Hruschka, A., Kurzrock, T. and Kölling-Speer, I. (2000). 25. In Caffeinated Beverages. Health benefits, physiological effects, and chemistry (pp. 241-251). ACS Symposium series 754. [ISBN:0-841-3654-2  ] " |
        referenceNew == "Speer, K., Hruschka, A., Kurzrock, T. and Kölling-Speer, I.  (2000). 25. In Caffeinated Beverages. Health benefits, physiological effects, and chemistry (pp. 241-251). ACS Symposium series 754. [ISBN:0-841-3654-2  ] ",
      yes = referenceNew,
      no = sub(
        pattern = ":",
        replacement = "§",
        x = referenceNew
      )
    )
  ) %>%
  cSplit("referenceNew",
    sep = "§",
    fixed = TRUE,
    drop = FALSE
  ) %>%
  cSplit(
    "referenceNew_1",
    sep = "[ISBN:",
    fixed = TRUE,
    stripWhite = FALSE
  ) %>%
  cSplit(
    "referenceNew_2",
    sep = "[ISBN:",
    fixed = TRUE,
    stripWhite = FALSE
  ) %>%
  mutate(referenceNew_2_1 = sub(
    pattern = "doi: ",
    replacement = "doi§ ",
    x = referenceNew_2_1
  )) %>%
  cSplit(
    "referenceNew_2_1",
    sep = "doi§ ",
    fixed = TRUE,
    stripWhite = FALSE,
  ) %>%
  mutate(
    referenceNew_2_1_1 = sub(
      pattern = "[PubMed:",
      replacement = "[PubMed§",
      fixed = TRUE,
      x = referenceNew_2_1_1
    ),
    referenceNew_2_1_2 = sub(
      pattern = "[PubMed:",
      replacement = "[PubMed§",
      fixed = TRUE,
      x = referenceNew_2_1_2
    )
  ) %>%
  cSplit(
    "referenceNew_2_1_1",
    sep = "[PubMed§",
    fixed = TRUE,
    stripWhite = FALSE
  ) %>%
  cSplit(
    "referenceNew_2_1_2",
    sep = "[PubMed§",
    fixed = TRUE,
    stripWhite = FALSE
  ) %>%
  mutate_all(as.character) %>%
  mutate(
    referenceAuthors_2 = ifelse(
      test = !is.na(referenceIsbn_2),
      yes = NA,
      no = referenceNew_1_1
    ),
    referenceUnsplittable_2 = referenceNew_2_1_1_1,
    referenceDoi_2_temp = referenceNew_2_1_2_1,
    referencePubmed_2_temp = ifelse(
      test = !is.na(referenceNew_2_1_2_2),
      yes = referenceNew_2_1_2_2,
      no = referenceNew_2_1_1_2
    )
  ) %>%
  mutate(
    referenceDoi_2_temp = sub(
      pattern = " ",
      replacement = "§",
      x = referenceDoi_2_temp
    ),
    referencePubmed_2_temp = sub(
      pattern = "  ]",
      replacement = "§",
      x = referencePubmed_2_temp
    )
  ) %>%
  cSplit("referenceDoi_2_temp", sep = "§") %>%
  cSplit("referencePubmed_2_temp", sep = "§") %>%
  mutate_all(as.character) %>%
  select(
    name,
    inchi,
    smiles,
    biologicalsource,
    reference,
    referenceIsbn_1,
    referenceAuthors_1,
    referenceUnsplittable_1,
    referenceDoi_1,
    referencePubmed_1,
    referenceNew,
    referenceIsbn_2,
    referenceAuthors_2,
    referenceUnsplittable_2,
    referenceDoi_2 = referenceDoi_2_temp_1,
    referencePubmed_2 = referencePubmed_2_temp_1,
    referenceNewNew = referencePubmed_2_temp_2,
  ) %>%
  mutate(
    referenceNewNew = ifelse(
      test = referenceNewNew == "Belitz, H.-D., Grosch, W., Schieberle, P. (2009). Food Chemistry. Springer. [ISBN:3540408185  ] " |
        referenceNewNew == "21. (2009). In Food Chemistry (pp. 938). Springer. [ISBN:3540408185  ] " |
        referenceNewNew == "21. (2009). In Food Chemistry (pp. 943). Springer. [ISBN:3540408185  ] FooDB [Link] " |
        referenceNewNew == "Speer, K., Hruschka, A., Kurzrock, T. and Kölling-Speer, I. (2000). 25. In Caffeinated Beverages. Health benefits, physiological effects, and chemistry (pp. 241-251). ACS Symposium series 754. [ISBN:0-841-3654-2  ] " |
        referenceNewNew == "Speer, K., Hruschka, A., Kurzrock, T. and Kölling-Speer, I.  (2000). 25. In Caffeinated Beverages. Health benefits, physiological effects, and chemistry (pp. 241-251). ACS Symposium series 754. [ISBN:0-841-3654-2  ] ",
      yes = referenceNewNew,
      no = sub(
        pattern = ":",
        replacement = "§",
        x = referenceNewNew
      )
    )
  ) %>%
  cSplit("referenceNewNew",
    sep = "§",
    fixed = TRUE,
    drop = FALSE
  ) %>%
  cSplit(
    "referenceNewNew_1",
    sep = "[ISBN:",
    fixed = TRUE,
    stripWhite = FALSE
  ) %>%
  cSplit(
    "referenceNewNew_2",
    sep = "[ISBN:",
    fixed = TRUE,
    stripWhite = FALSE
  ) %>%
  mutate(referenceNewNew_2_1 = sub(
    pattern = "doi: ",
    replacement = "doi§ ",
    x = referenceNewNew_2_1
  )) %>%
  cSplit(
    "referenceNewNew_2_1",
    sep = "doi§ ",
    fixed = TRUE,
    stripWhite = FALSE,
  ) %>%
  mutate(
    referenceNewNew_2_1_1 = sub(
      pattern = "[PubMed:",
      replacement = "[PubMed§",
      fixed = TRUE,
      x = referenceNewNew_2_1_1
    ),
    referenceNewNew_2_1_2 = sub(
      pattern = "[PubMed:",
      replacement = "[PubMed§",
      fixed = TRUE,
      x = referenceNewNew_2_1_2
    )
  ) %>%
  cSplit(
    "referenceNewNew_2_1_1",
    sep = "[PubMed§",
    fixed = TRUE,
    stripWhite = FALSE
  ) %>%
  cSplit(
    "referenceNewNew_2_1_2",
    sep = "[PubMed§",
    fixed = TRUE,
    stripWhite = FALSE
  ) %>%
  mutate_all(as.character) %>%
  mutate(
    referenceAuthors_3 = referenceNewNew_1_1,
    referenceUnsplittable_3 = referenceNewNew_2_1_1_1,
    referenceDoi_3_temp = referenceNewNew_2_1_2_1,
    referencePubmed_3_temp = ifelse(
      test = !is.na(referenceNewNew_2_1_2_2),
      yes = referenceNewNew_2_1_2_2,
      no = referenceNewNew_2_1_1_2
    )
  ) %>%
  mutate(
    referenceDoi_3_temp = sub(
      pattern = " ",
      replacement = "§",
      x = referenceDoi_3_temp
    ),
    referencePubmed_3_temp = sub(
      pattern = "  ]",
      replacement = "§",
      x = referencePubmed_3_temp
    )
  ) %>%
  cSplit("referenceDoi_3_temp", sep = "§") %>%
  cSplit("referencePubmed_3_temp", sep = "§") %>%
  mutate_all(as.character) %>%
  select(
    name,
    inchi,
    smiles,
    biologicalsource,
    reference,
    referenceIsbn_1,
    referenceAuthors_1,
    referenceUnsplittable_1,
    referenceDoi_1,
    referencePubmed_1,
    referenceIsbn_2,
    referenceAuthors_2,
    referenceUnsplittable_2,
    referenceDoi_2,
    referencePubmed_2,
    referenceAuthors_3,
    referenceUnsplittable_3,
    referenceDoi_3 = referenceDoi_3_temp_1,
    referencePubmed_3 = referencePubmed_3_temp_1
  )

data_pivoted <- data_manipulated %>%
  pivot_longer(
    6:ncol(.),
    names_to = c(".value", "level"),
    names_sep = "_",
    values_to = "reference",
    values_drop_na = TRUE
  ) %>%
  select(
    name,
    inchi,
    smiles,
    biologicalsource,
    reference_isbn = referenceIsbn,
    reference_authors = referenceAuthors,
    reference_original = referenceUnsplittable,
    reference_doi = referenceDoi,
    reference_pubmed = referencePubmed
  ) %>%
  mutate(reference_external = ifelse(
    test = grepl(x = reference_authors, pattern = "[Link]", fixed = TRUE),
    yes = reference_authors,
    no = NA
  )) %>%
  mutate(reference_external = ifelse(
    test = grepl(x = reference_isbn, pattern = "[Link]", fixed = TRUE),
    yes = gsub(
      pattern = "[0-9]",
      replacement = "",
      x = reference_isbn
    ),
    no = reference_external
  )) %>%
  mutate(reference_authors = ifelse(
    test = !is.na(reference_external),
    yes = NA,
    no = reference_authors
  )) %>%
  mutate(reference_isbn = gsub(
    pattern = " .*",
    replacement = "",
    x = reference_isbn
  )) %>%
  data.frame()

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_pivoted,
    db = "phy_1",
    structure_field = c("name", "inchi", "smiles"),
    reference_field = c(
      "reference_isbn",
      "reference_authors",
      "reference_original",
      "reference_doi",
      "reference_pubmed",
      "reference_external"
    )
  )

# exporting
database$writeInterim(data_standard)
