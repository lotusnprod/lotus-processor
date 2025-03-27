# title: "SANCDB cleaneR"

# loading paths
source("paths.R")
source("r/y_as_na.R")
source("r/standardizing_original.R")

library(dplyr)
library(readr)
library(stringr)
library(tidyr)

# get paths
database <- databases$get("sancdb")

## files
data_original <- readr::read_delim(
  file = gzfile(description = database$sourceFiles$tsv),
  col_types = cols(.default = "c")
)

#  cleaning
## function
internal_clean <- function(df, num, var) {
  colname_1 <-
    paste(
      "V1",
      stringr::str_pad(
        string = num,
        width = 3,
        pad = 0
      ),
      sep = "_"
    )
  colname_2 <-
    paste(
      "V1",
      stringr::str_pad(
        string = num + 1,
        width = 3,
        pad = 0
      ),
      sep = "_"
    )
  df1 <- df |>
    dplyr::filter(!!as.name(colname_1) != var)
  df1[, (num + 2):ncol(df)] <- df1[, num:(ncol(df) - 2)]
  df1 <- df1 |>
    dplyr::mutate(
      !!as.name(colname_1) := NA_character_,
      !!as.name(colname_2) := NA_character_,
    )
  df2 <- df |>
    dplyr::filter(!!as.name(colname_1) == var)
  df_clean <- df1 |>
    dplyr::bind_rows(df2)
  return(df_clean)
}

int_clean_min <- function(df, num, var) {
  colname_1 <-
    paste(
      "V1",
      stringr::str_pad(
        string = num,
        width = 3,
        pad = 0
      ),
      sep = "_"
    )
  colname_2 <-
    paste(
      "V1",
      stringr::str_pad(
        string = num + 1,
        width = 3,
        pad = 0
      ),
      sep = "_"
    )
  df1 <- df |>
    dplyr::filter(!!as.name(colname_1) == var)
  df1[, (num + 1):ncol(df)] <- df1[, num:(ncol(df) - 1)]
  df1 <- df1 |>
    dplyr::mutate(!!as.name(colname_1) := NA_character_)
  df2 <- df |>
    dplyr::filter(
      !!as.name(colname_1) != var |
        !!as.name(colname_2) == var
    )
  df_clean <- df1 |>
    dplyr::bind_rows(df2)
  return(df_clean)
}

int_clean_rev <- function(df, num, var) {
  colname_1 <-
    paste(
      "V1",
      stringr::str_pad(
        string = num,
        width = 3,
        pad = 0
      ),
      sep = "_"
    )
  df1 <- df |>
    dplyr::filter(!!as.name(colname_1) != var)
  df1[, (num):(ncol(df) - 1)] <- df1[, (num + 1):(ncol(df))]
  df2 <- df |>
    dplyr::filter(!!as.name(colname_1) == var)
  df_clean <- df1 |>
    dplyr::bind_rows(df2)
  return(df_clean)
}

SANCDB_clean <- function(dfsel) {
  df1 <- dfsel |>
    dplyr::filter(V1_002 == "Entry name:")

  df2 <- internal_clean(df = df1, num = 6, var = "Molecular mass:")
  df3 <- internal_clean(df = df2, num = 8, var = "ChEMBL ID:")
  df4 <- internal_clean(df = df3, num = 10, var = "ChemSpider ID:")
  df5 <- internal_clean(df = df4, num = 12, var = "ZINC ID:")
  df6 <- internal_clean(df = df5, num = 14, var = "PubChem ID:")
  df7 <- internal_clean(df = df6, num = 16, var = "DrugBank ID:")
  df8 <- int_clean_min(df = df7, num = 17, var = "CAS No.")
  df9 <- int_clean_rev(df = df8, num = 33, var = "References")
  df10 <- int_clean_rev(df = df9, num = 33, var = "References")
  df11 <- int_clean_rev(df = df10, num = 33, var = "References")
  df12 <- int_clean_rev(df = df11, num = 33, var = "References")
  df13 <- int_clean_rev(df = df12, num = 33, var = "References")
  df14 <- int_clean_rev(df = df13, num = 33, var = "References")
  df15 <- int_clean_rev(df = df14, num = 33, var = "References")
  df16 <- int_clean_rev(df = df15, num = 33, var = "References")
  df17 <- int_clean_rev(df = df16, num = 33, var = "References")
  df18 <- int_clean_rev(df = df17, num = 33, var = "References")
  df19 <- int_clean_rev(df = df18, num = 33, var = "References")
  df20 <- int_clean_rev(df = df19, num = 33, var = "References")
  df21 <- int_clean_rev(df = df20, num = 33, var = "References")
  df22 <- int_clean_rev(df = df21, num = 33, var = "References")
  df23 <- int_clean_rev(df = df22, num = 33, var = "References")
  df24 <- int_clean_rev(df = df23, num = 33, var = "References")
  df25 <- int_clean_rev(df = df24, num = 33, var = "References")
  df26 <- int_clean_rev(df = df25, num = 33, var = "References")
  df27 <- int_clean_rev(df = df26, num = 33, var = "References")
  df28 <- int_clean_rev(df = df27, num = 33, var = "References")
  df29 <-
    int_clean_min(df = df28, num = 36, var = "Classifications")
  df30 <-
    int_clean_min(df = df29, num = 37, var = "Classifications")
  df31 <-
    int_clean_min(df = df30, num = 38, var = "Classifications")
  df32 <-
    int_clean_min(df = df31, num = 39, var = "Classifications")
  df33 <-
    int_clean_min(df = df32, num = 40, var = "Classifications")
  df34 <-
    int_clean_min(df = df33, num = 41, var = "Classifications")
  df35 <-
    int_clean_min(df = df34, num = 42, var = "Classifications")
  df36 <-
    int_clean_min(df = df35, num = 43, var = "Classifications")
  df37 <-
    int_clean_min(df = df36, num = 44, var = "Classifications")
  df38 <-
    int_clean_min(df = df37, num = 45, var = "Classifications")
  df39 <- int_clean_min(df = df38, num = 48, var = "Other Names")
  df40 <- int_clean_min(df = df39, num = 49, var = "Other Names")
  df41 <- int_clean_min(df = df40, num = 50, var = "Other Names")
  df42 <- int_clean_min(df = df41, num = 51, var = "Other Names")
  df43 <-
    int_clean_min(df = df42, num = 54, var = "Source Organisms")
  df44 <-
    int_clean_min(df = df43, num = 55, var = "Source Organisms")
  df45 <-
    int_clean_min(df = df44, num = 56, var = "Source Organisms")
  df46 <-
    int_clean_min(df = df45, num = 57, var = "Source Organisms")
  df47 <-
    int_clean_min(df = df46, num = 58, var = "Source Organisms")
  df48 <-
    int_clean_min(df = df47, num = 59, var = "Source Organisms")
  df49 <-
    int_clean_min(df = df48, num = 60, var = "Source Organisms")
  df50 <-
    int_clean_min(df = df49, num = 61, var = "Source Organisms")
  df51 <-
    int_clean_min(df = df50, num = 62, var = "Source Organisms")
  df52 <-
    int_clean_min(df = df51, num = 63, var = "Source Organisms")
  df53 <-
    int_clean_min(df = df52, num = 64, var = "Source Organisms")
  df54 <-
    int_clean_min(df = df53, num = 65, var = "Source Organisms")
  df55 <-
    int_clean_min(df = df54, num = 66, var = "Source Organisms")
  df56 <-
    int_clean_min(df = df55, num = 67, var = "Source Organisms")
  df57 <-
    int_clean_min(df = df56, num = 68, var = "Source Organisms")
  df58 <-
    int_clean_min(df = df57, num = 69, var = "Source Organisms")
  df59 <-
    int_clean_min(df = df58, num = 70, var = "Source Organisms")
  df60 <-
    int_clean_min(df = df59, num = 71, var = "Source Organisms")
  df61 <-
    int_clean_min(df = df60, num = 72, var = "Source Organisms")
  df62 <-
    int_clean_min(df = df61, num = 73, var = "Source Organisms")
  df63 <-
    int_clean_min(df = df62, num = 74, var = "Source Organisms")
  df64 <-
    int_clean_min(df = df63, num = 75, var = "Source Organisms")
  df65 <-
    int_clean_min(df = df64, num = 76, var = "Source Organisms")
  df66 <-
    int_clean_min(df = df65, num = 77, var = "Source Organisms")
  df67 <-
    int_clean_min(df = df66, num = 78, var = "Source Organisms")
  df68 <-
    int_clean_min(df = df67, num = 79, var = "Source Organisms")
  df69 <-
    int_clean_min(df = df68, num = 80, var = "Source Organisms")
  df70 <-
    int_clean_min(df = df69, num = 81, var = "Source Organisms")
  df71 <-
    int_clean_min(df = df70, num = 82, var = "Source Organisms")
  df72 <-
    int_clean_min(df = df71, num = 83, var = "Source Organisms")
  df73 <-
    int_clean_min(df = df72, num = 84, var = "Source Organisms")
  df74 <-
    int_clean_min(df = df73, num = 85, var = "Source Organisms")
  df75 <-
    int_clean_min(df = df74, num = 86, var = "Source Organisms")
  df76 <-
    int_clean_min(df = df75, num = 87, var = "Source Organisms")
  df77 <-
    int_clean_min(df = df76, num = 88, var = "Source Organisms")
  df78 <-
    int_clean_min(df = df77, num = 89, var = "Source Organisms")
  df79 <-
    int_clean_min(df = df78, num = 90, var = "Source Organisms")
  df80 <-
    int_clean_min(df = df79, num = 91, var = "Source Organisms")
  df81 <-
    int_clean_min(df = df80, num = 92, var = "Source Organisms")
  df82 <-
    int_clean_min(df = df81, num = 93, var = "Source Organisms")
  df83 <-
    int_clean_min(df = df82, num = 94, var = "Source Organisms")
  df84 <-
    int_clean_min(df = df83, num = 95, var = "Source Organisms")
  df85 <-
    int_clean_min(df = df84, num = 96, var = "Source Organisms")
  df86 <-
    int_clean_min(df = df85, num = 97, var = "Source Organisms")
  df87 <-
    int_clean_min(df = df86, num = 98, var = "Source Organisms")
  df88 <-
    int_clean_min(df = df87, num = 99, var = "Source Organisms")
  df89 <-
    int_clean_min(df = df88, num = 100, var = "Source Organisms")
  df90 <-
    int_clean_min(df = df89, num = 101, var = "Source Organisms")
  df91 <-
    int_clean_min(df = df90, num = 102, var = "Source Organisms")
  df92 <-
    int_clean_min(df = df91, num = 103, var = "Source Organisms")
  df93 <-
    int_clean_min(df = df92, num = 104, var = "Source Organisms")
  df94 <-
    int_clean_min(df = df93, num = 105, var = "Source Organisms")
  df95 <-
    int_clean_min(df = df94, num = 106, var = "Source Organisms")
  df96 <-
    int_clean_min(df = df95, num = 109, var = "Compound Uses")
  df97 <-
    int_clean_min(df = df96, num = 110, var = "Compound Uses")
  df98 <-
    int_clean_min(df = df97, num = 111, var = "Compound Uses")
  df99 <- df98 |>
    dplyr::mutate(
      V1_112 = ifelse(
        test = V1_112 == "Compound Uses",
        yes = NA_character_,
        no = V1_112
      )
    )

  df_selected <- df99 |>
    dplyr::select(
      uniqueid = 1,
      name = 3,
      formula = 5,
      mass = 7,
      smiles = 32,
      reference_1 = 34,
      reference_2 = 35,
      reference_3 = 36,
      reference_4 = 37,
      reference_5 = 38,
      reference_6 = 39,
      reference_7 = 40,
      reference_8 = 41,
      reference_9 = 42,
      reference_10 = 43,
      reference_11 = 44,
      reference_12 = 45,
      biologicalsource_1 = 108,
      biologicalsource_2 = 109,
      biologicalsource_3 = 110,
      biologicalsource_4 = 111,
      biologicalsource_5 = 112
    )

  df_pivoted <- df_selected |>
    tidyr::pivot_longer(
      cols = 18:22,
      names_prefix = "biologicalsource_",
      names_repair = "minimal"
    ) |>
    data.frame() |>
    dplyr::filter(!is.na(value)) |>
    dplyr::select(biologicalsource = value, dplyr::everything()) |>
    tidyr::pivot_longer(
      cols = 7:18,
      names_prefix = "reference_",
      names_repair = "minimal"
    ) |>
    data.frame() |>
    dplyr::filter(!is.na(value)) |>
    dplyr::filter(grepl(pattern = "^\\([0-9]{4}\\)", x = value)) |>
    dplyr::mutate(
      reference = gsub(
        pattern = "^\\([0-9]{4}\\) ",
        replacement = "",
        x = value
      )
    )

  return(df_pivoted)
}

## applying
sancdb_clean <- SANCDB_clean(dfsel = data_original)

# selecting
data_manipulated <- sancdb_clean |>
  dplyr::select(
    uniqueid,
    structure_name = name,
    structure_smiles = smiles,
    organism_clean = biologicalsource,
    reference_title = reference
  ) |>
  data.frame()

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_manipulated,
    db = "sancdb",
    structure_field = "structure_smiles",
    organism_field = "organism_clean",
    reference_field = c("reference_title")
  )

# exporting
database$writeInterim(data_standard)
