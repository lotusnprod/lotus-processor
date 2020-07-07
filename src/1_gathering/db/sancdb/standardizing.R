# title: "SANCDB cleaneR"

# loading paths
source("paths.R")
source("functions/helpers.R")
source("functions/standardizing.R")

library(dplyr)
library(readr)
library(splitstackshape)
library(stringr)
library(tidyr)

# get paths
database <- databases$get("sancdb")

## files
data_original <- read_delim(
  file = gzfile(database$sourceFiles$tsv),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

#  cleaning
## function
SANCDB_clean <- function(dfsel)
{
  df1 <- dfsel %>% filter(V1_02 == "Entry name:")
  
  df11 <- df1 %>%
    filter(V1_06 != "ChEMBL ID:")
  
  for (i in (ncol(df11):7))
  {
    df11[, i] <- df11[, i - 2]
  }
  
  df11$V1_06 <- ""
  df11$V1_07 <- ""
  df11$V1_06 <- y_as_na(x = df11$V1_06, y = "")
  df11$V1_07 <- y_as_na(x = df11$V1_07, y = "")
  df11$V1_06 <- as.character(df11$V1_06)
  df11$V1_07 <- as.character(df11$V1_07)
  
  df12 <- df1 %>%
    filter(V1_06 == "ChEMBL ID:")
  
  df2 <- full_join(df12, df11)
  
  df21 <- df2 %>%
    filter(V1_08 != "ChemSpider ID:")
  
  for (i in (ncol(df21):9))
  {
    df21[, i] <- df21[, i - 2]
  }
  
  df21$V1_08 <- ""
  df21$V1_09 <- ""
  df21$V1_08 <- y_as_na(x = df21$V1_08, y = "")
  df21$V1_09 <- y_as_na(x = df21$V1_09, y = "")
  df21$V1_08 <- as.character(df21$V1_08)
  df21$V1_09 <- as.character(df21$V1_09)
  
  df22 <- df2 %>%
    filter(V1_08 == "ChemSpider ID:")
  
  df3 <- full_join(df22, df21)
  
  df31 <- df3 %>%
    filter(V1_10 != "ZINC ID:")
  
  for (i in (ncol(df31):11))
  {
    df31[, i] <- df31[, i - 2]
  }
  
  df31$V1_10 <- ""
  df31$V1_11 <- ""
  df31$V1_10 <- y_as_na(x = df31$V1_10, y = "")
  df31$V1_11 <- y_as_na(x = df31$V1_11, y = "")
  df31$V1_10 <- as.character(df31$V1_10)
  df31$V1_11 <- as.character(df31$V1_11)
  
  df32 <- df3 %>%
    filter(V1_10 == "ZINC ID:")
  
  df4 <- full_join(df32, df31)
  
  df41 <- df4 %>%
    filter(V1_12 != "PubChem ID:")
  
  for (i in (ncol(df41):13))
  {
    df41[, i] <- df41[, i - 2]
  }
  
  df41$V1_12 <- ""
  df41$V1_13 <- ""
  df41$V1_12 <- y_as_na(x = df41$V1_12, y = "")
  df41$V1_13 <- y_as_na(x = df41$V1_13, y = "")
  df41$V1_12 <- as.character(df41$V1_12)
  df41$V1_13 <- as.character(df41$V1_13)
  
  df42 <- df4 %>%
    filter(V1_12 == "PubChem ID:")
  
  df5 <- full_join(df42, df41)
  
  df51 <- df5 %>%
    filter(V1_14 != "DrukBank ID:")
  
  for (i in (ncol(df51):15))
  {
    df51[, i] <- df51[, i - 2]
  }
  
  df51$V1_14 <- ""
  df51$V1_15 <- ""
  df51$V1_14 <- y_as_na(x = df51$V1_14, y = "")
  df51$V1_15 <- y_as_na(x = df51$V1_15, y = "")
  df51$V1_14 <- as.character(df51$V1_14)
  df51$V1_15 <- as.character(df51$V1_15)
  
  df52 <- df5 %>%
    filter(V1_14 == "DrukBank ID:")
  
  df6 <- full_join(df52, df51)
  
  df61 <- df6 %>%
    filter(V1_31 == "Classifications")
  
  for (i in (ncol(df61):32))
  {
    df61[, i] <- df61[, i - 1]
  }
  
  df61$V1_31 <- ""
  df61$V1_31 <- y_as_na(x = df61$V1_31, y = "")
  df61$V1_31 <- as.character(df61$V1_31)
  
  df62 <- df6 %>%
    filter(V1_31 != "Classifications")
  
  df7 <- full_join(df62, df61)
  
  df71 <- df7 %>%
    filter(V1_32 == "Classifications")
  
  for (i in (ncol(df71):33))
  {
    df71[, i] <- df71[, i - 1]
  }
  
  df71$V1_32 <- ""
  df71$V1_32 <- y_as_na(x = df71$V1_32, y = "")
  df71$V1_32 <- as.character(df71$V1_32)
  
  df72 <- df7 %>%
    filter(V1_32 != "Classifications")
  
  df8 <- full_join(df72, df71)
  
  df81 <- df8 %>%
    filter(V1_33 == "Classifications")
  
  for (i in (ncol(df81):34))
  {
    df81[, i] <- df81[, i - 1]
  }
  
  df81$V1_33 <- ""
  df81$V1_33 <- y_as_na(x = df81$V1_33, y = "")
  df81$V1_33 <- as.character(df81$V1_33)
  
  df82 <- df8 %>%
    filter(V1_33 != "Classifications")
  
  df9 <- full_join(df82, df81)
  
  df91 <- df9 %>%
    filter(V1_34 == "Classifications")
  
  for (i in (ncol(df91):35))
  {
    df91[, i] <- df91[, i - 1]
  }
  
  df91$V1_34 <- ""
  df91$V1_34 <- y_as_na(x = df91$V1_34, y = "")
  df91$V1_34 <- as.character(df91$V1_34)
  
  df92 <- df9 %>%
    filter(V1_34 != "Classifications")
  
  df10 <- full_join(df92, df91)
  
  dfa1 <- df10 %>%
    filter(V1_35 == "Classifications")
  
  for (i in (ncol(dfa1):36))
  {
    dfa1[, i] <- dfa1[, i - 1]
  }
  
  dfa1$V1_35 <- ""
  dfa1$V1_35 <- y_as_na(x = dfa1$V1_35, y = "")
  dfa1$V1_35 <- as.character(dfa1$V1_35)
  
  dfa2 <- df10 %>%
    filter(V1_35 != "Classifications")
  
  dfa <- full_join(dfa2, dfa1)
  
  dfb1 <- dfa %>%
    filter(V1_36 == "Classifications")
  
  for (i in (ncol(dfb1):37))
  {
    dfb1[, i] <- dfb1[, i - 1]
  }
  
  dfb1$V1_36 <- ""
  dfb1$V1_36 <- y_as_na(x = dfb1$V1_36, y = "")
  dfb1$V1_36 <- as.character(dfb1$V1_36)
  
  dfb2 <- dfa %>%
    filter(V1_36 != "Classifications")
  
  dfb <- full_join(dfb2, dfb1)
  
  dfc1 <- dfb %>%
    filter(V1_37 == "Classifications")
  
  for (i in (ncol(dfc1):38))
  {
    dfc1[, i] <- dfc1[, i - 1]
  }
  
  dfc1$V1_37 <- ""
  dfc1$V1_37 <- y_as_na(x = dfc1$V1_37, y = "")
  dfc1$V1_37 <- as.character(dfc1$V1_37)
  
  dfc2 <- dfb %>%
    filter(V1_37 != "Classifications")
  
  dfc <- full_join(dfc2, dfc1)
  
  dfd1 <- dfc %>%
    filter(V1_38 == "Classifications")
  
  for (i in (ncol(dfd1):39))
  {
    dfd1[, i] <- dfd1[, i - 1]
  }
  
  dfd1$V1_38 <- ""
  dfd1$V1_38 <- y_as_na(x = dfd1$V1_38, y = "")
  dfd1$V1_38 <- as.character(dfd1$V1_38)
  
  dfd2 <- dfc %>%
    filter(V1_38 != "Classifications")
  
  dfd <- full_join(dfd2, dfd1)
  
  dfe1 <- dfd %>%
    filter(V1_39 == "Classifications")
  
  for (i in (ncol(dfe1):40))
  {
    dfe1[, i] <- dfe1[, i - 1]
  }
  
  dfe1$V1_39 <- ""
  dfe1$V1_39 <- y_as_na(x = dfe1$V1_39, y = "")
  dfe1$V1_39 <- as.character(dfe1$V1_39)
  
  dfe2 <- dfd %>%
    filter(V1_39 != "Classifications")
  
  dfe <- full_join(dfe2, dfe1)
  
  dff1 <- dfe %>%
    filter(V1_40 == "Classifications")
  
  for (i in (ncol(dff1):41))
  {
    dff1[, i] <- dff1[, i - 1]
  }
  
  dff1$V1_40 <- ""
  dff1$V1_40 <- y_as_na(x = dff1$V1_40, y = "")
  dff1$V1_40 <- as.character(dff1$V1_40)
  
  dff2 <- dfe %>%
    filter(V1_40 != "Classifications")
  
  dff <- full_join(dff2, dff1)
  
  dfg1 <- dff %>%
    filter(V1_43 == "Other Names")
  
  for (i in (ncol(dfg1):44))
  {
    dfg1[, i] <- dfg1[, i - 1]
  }
  
  dfg1$V1_43 <- y_as_na(x = dfg1$V1_43, y = "Other Names")
  dfg1$V1_43 <- as.character(dfg1$V1_43)
  
  dfg2 <- dff %>%
    filter(V1_43 != "Other Names")
  
  dfg <- full_join(dfg2, dfg1)
  
  dfh1 <- dfg %>%
    filter(V1_44 == "Other Names")
  
  for (i in (ncol(dfh1):45))
  {
    dfh1[, i] <- dfh1[, i - 1]
  }
  
  dfh1$V1_44 <- y_as_na(x = dfh1$V1_44, y = "Other Names")
  dfh1$V1_44 <- as.character(dfh1$V1_44)
  
  dfh2 <- dfg %>%
    filter(V1_44 != "Other Names")
  
  dfh <- full_join(dfh2, dfh1)
  
  dfi1 <- dfh %>%
    filter(V1_45 == "Other Names")
  
  for (i in (ncol(dfi1):46))
  {
    dfi1[, i] <- dfi1[, i - 1]
  }
  
  dfi1$V1_45 <- y_as_na(x = dfi1$V1_45, y = "Other Names")
  dfi1$V1_45 <- as.character(dfi1$V1_45)
  
  dfi2 <- dfh %>%
    filter(V1_45 != "Other Names")
  
  dfi <- full_join(dfi2, dfi1)
  
  dfj1 <- dfi %>%
    filter(V1_47 == "Source Organisms")
  
  for (i in (ncol(dfj1):48))
  {
    dfj1[, i] <- dfj1[, i - 1]
  }
  
  dfj1$V1_47 <- y_as_na(x = dfj1$V1_47, y = "Source Organisms")
  dfj1$V1_47 <- as.character(dfj1$V1_47)
  
  dfj2 <- dfi %>%
    filter(V1_47 != "Source Organisms")
  
  dfj <- full_join(dfj2, dfj1)
  
  dfk1 <- dfj %>%
    filter(V1_48 == "Source Organisms")
  
  for (i in (ncol(dfk1):49))
  {
    dfk1[, i] <- dfk1[, i - 1]
  }
  
  dfk1$V1_48 <- y_as_na(x = dfk1$V1_48, y = "Source Organisms")
  dfk1$V1_48 <- as.character(dfk1$V1_48)
  
  dfk2 <- dfj %>%
    filter(V1_48 != "Source Organisms")
  
  dfk <- full_join(dfk2, dfk1)
  
  dfl1 <- dfk %>%
    filter(V1_49 == "Source Organisms")
  
  for (i in (ncol(dfl1):50))
  {
    dfl1[, i] <- dfl1[, i - 1]
  }
  
  dfl1$V1_49 <- y_as_na(x = dfl1$V1_49, y = "Source Organisms")
  dfl1$V1_49 <- as.character(dfl1$V1_49)
  
  dfl2 <- dfk %>%
    filter(V1_49 != "Source Organisms")
  
  dfl <- full_join(dfl2, dfl1)
  
  dfm1 <- dfl %>%
    filter(V1_50 == "Source Organisms")
  
  for (i in (ncol(dfm1):51))
  {
    dfm1[, i] <- dfm1[, i - 1]
  }
  
  dfm1$V1_50 <- y_as_na(x = dfm1$V1_50, y = "Source Organisms")
  dfm1$V1_50 <- as.character(dfm1$V1_50)
  
  dfm2 <- dfl %>%
    filter(V1_50 != "Source Organisms")
  
  dfm <- full_join(dfm2, dfm1)
  
  dfn1 <- dfm %>%
    filter(V1_51 == "Source Organisms")
  
  for (i in (ncol(dfn1):52))
  {
    dfn1[, i] <- dfn1[, i - 1]
  }
  
  dfn1$V1_51 <- y_as_na(x = dfn1$V1_51, y = "Source Organisms")
  dfn1$V1_51 <- as.character(dfn1$V1_51)
  
  dfn2 <- dfm %>%
    filter(V1_51 != "Source Organisms")
  
  dfn <- full_join(dfn2, dfn1)
  
  dfo1 <- dfn %>%
    filter(V1_52 == "Source Organisms")
  
  for (i in (ncol(dfo1):53))
  {
    dfo1[, i] <- dfo1[, i - 1]
  }
  
  dfo1$V1_52 <- y_as_na(x = dfo1$V1_52, y = "Source Organisms")
  dfo1$V1_52 <- as.character(dfo1$V1_52)
  
  dfo2 <- dfn %>%
    filter(V1_52 != "Source Organisms")
  
  dfo <- full_join(dfo2, dfo1)
  
  dfp1 <- dfo %>%
    filter(V1_53 == "Source Organisms")
  
  for (i in (ncol(dfp1):54))
  {
    dfp1[, i] <- dfp1[, i - 1]
  }
  
  dfp1$V1_53 <- y_as_na(x = dfp1$V1_53, y = "Source Organisms")
  dfp1$V1_53 <- as.character(dfp1$V1_53)
  
  dfp2 <- dfo %>%
    filter(V1_53 != "Source Organisms")
  
  dfp <- full_join(dfp2, dfp1)
  
  dfq1 <- dfp %>%
    filter(V1_54 == "Source Organisms")
  
  for (i in (ncol(dfq1):55))
  {
    dfq1[, i] <- dfq1[, i - 1]
  }
  
  dfq1$V1_54 <- y_as_na(x = dfq1$V1_54, y = "Source Organisms")
  dfq1$V1_54 <- as.character(dfq1$V1_54)
  
  dfq2 <- dfp %>%
    filter(V1_54 != "Source Organisms")
  
  dfq <- full_join(dfq2, dfq1)
  
  dfr1 <- dfq %>%
    filter(V1_55 == "Source Organisms")
  
  for (i in (ncol(dfr1):56))
  {
    dfr1[, i] <- dfr1[, i - 1]
  }
  
  dfr1$V1_55 <- y_as_na(x = dfr1$V1_55, y = "Source Organisms")
  dfr1$V1_55 <- as.character(dfr1$V1_55)
  
  dfr2 <- dfq %>%
    filter(V1_55 != "Source Organisms")
  
  dfr <- full_join(dfr2, dfr1)
  
  dfs1 <- dfr %>%
    filter(V1_56 == "Source Organisms")
  
  for (i in (ncol(dfs1):57))
  {
    dfs1[, i] <- dfs1[, i - 1]
  }
  
  dfs1$V1_56 <- y_as_na(x = dfs1$V1_56, y = "Source Organisms")
  dfs1$V1_56 <- as.character(dfs1$V1_56)
  
  dfs2 <- dfr %>%
    filter(V1_56 != "Source Organisms")
  
  dfs <- full_join(dfs2, dfs1)
  
  dft1 <- dfs %>%
    filter(V1_57 == "Source Organisms")
  
  for (i in (ncol(dft1):58))
  {
    dft1[, i] <- dft1[, i - 1]
  }
  
  dft1$V1_57 <- y_as_na(x = dft1$V1_57, y = "Source Organisms")
  dft1$V1_57 <- as.character(dft1$V1_57)
  
  dft2 <- dfs %>%
    filter(V1_57 != "Source Organisms")
  
  dft <- full_join(dft2, dft1)
  
  dfu1 <- dft %>%
    filter(V1_58 == "Source Organisms")
  
  for (i in (ncol(dfu1):59))
  {
    dfu1[, i] <- dfu1[, i - 1]
  }
  
  dfu1$V1_58 <- y_as_na(x = dfu1$V1_58, y = "Source Organisms")
  dfu1$V1_58 <- as.character(dfu1$V1_58)
  
  dfu2 <- dft %>%
    filter(V1_58 != "Source Organisms")
  
  dfu <- full_join(dfu2, dfu1)
  
  dfv1 <- dfu %>%
    filter(V1_59 == "Source Organisms")
  
  for (i in (ncol(dfv1):60))
  {
    dfv1[, i] <- dfv1[, i - 1]
  }
  
  dfv1$V1_59 <- y_as_na(x = dfv1$V1_59, y = "Source Organisms")
  dfv1$V1_59 <- as.character(dfv1$V1_59)
  
  dfv2 <- dfu %>%
    filter(V1_59 != "Source Organisms")
  
  dfv <- full_join(dfv2, dfv1)
  
  dfw1 <- dfv %>%
    filter(V1_60 == "Source Organisms")
  
  for (i in (ncol(dfw1):61))
  {
    dfw1[, i] <- dfw1[, i - 1]
  }
  
  dfw1$V1_60 <- y_as_na(x = dfw1$V1_60, y = "Source Organisms")
  dfw1$V1_60 <- as.character(dfw1$V1_60)
  
  dfw2 <- dfv %>%
    filter(V1_60 != "Source Organisms")
  
  dfw <- full_join(dfw2, dfw1)
  
  dfx1 <- dfw %>%
    filter(V1_61 == "Source Organisms")
  
  for (i in (ncol(dfx1):62))
  {
    dfx1[, i] <- dfx1[, i - 1]
  }
  
  dfx1$V1_61 <- y_as_na(x = dfx1$V1_61, y = "Source Organisms")
  dfx1$V1_61 <- as.character(dfx1$V1_61)
  
  dfx2 <- dfw %>%
    filter(V1_61 != "Source Organisms")
  
  dfx <- full_join(dfx2, dfx1)
  
  dfy1 <- dfx %>%
    filter(V1_62 == "Source Organisms")
  
  for (i in (ncol(dfy1):63))
  {
    dfy1[, i] <- dfy1[, i - 1]
  }
  
  dfy1$V1_62 <- y_as_na(x = dfy1$V1_62, y = "Source Organisms")
  dfy1$V1_62 <- as.character(dfy1$V1_62)
  
  dfy2 <- dfx %>%
    filter(V1_62 != "Source Organisms")
  
  dfy <- full_join(dfy2, dfy1)
  
  dfz1 <- dfy %>%
    filter(V1_63 == "Source Organisms")
  
  for (i in (ncol(dfz1):64))
  {
    dfz1[, i] <- dfz1[, i - 1]
  }
  
  dfz1$V1_63 <- y_as_na(x = dfz1$V1_63, y = "Source Organisms")
  dfz1$V1_63 <- as.character(dfz1$V1_63)
  
  dfz2 <- dfy %>%
    filter(V1_63 != "Source Organisms")
  
  dfz2[1, 64:67] <- dfz2[1, 68:71]
  dfz2[2, 64:70] <- dfz2[2, 71:77]
  dfz2[3, 64:67] <- dfz2[3, 70:73]
  dfz2[4, 64:67] <- dfz2[4, 68:71]
  dfz2[5, 64:67] <- dfz2[5, 70:73]
  dfz2[6, 64:67] <- dfz2[6, 68:71]
  dfz2[7, 64:67] <- NA
  dfz2[8, 64:67] <- dfz2[8, 65:68]
  dfz2[9, 64:67] <- dfz2[9, 65:68]
  dfz2[10, 64:67] <- dfz2[10, 67:70]
  dfz2[11, 64:67] <- dfz2[11, 82:85]
  
  dfz <- full_join(dfz2, dfz1)
  return(dfz)
}

##applying
sancdb_clean <- SANCDB_clean(dfsel = data_original)

#selecting
data_selected <- sancdb_clean %>%
  select(
    uniqueid = 1,
    name = 3,
    smiles = 27,
    biologicalsource = 65,
    cas = 17,
    pubchem = 13,
    V1_29,
    V1_30,
    V1_31,
    V1_32,
    V1_33,
    V1_34,
    V1_35,
    V1_36,
    V1_37,
    V1_38,
    V1_39,
    V1_40
  )

data_selected_21 <- data_selected %>%
  filter(grepl(pattern = "^\\(", x = V1_30))

for (i in (ncol(data_selected_21):8))
{
  data_selected_21[, i] <- data_selected_21[, i - 1]
}

data_selected_21$V1_30 <- NA
data_selected_21$V1_30 <- as.character(data_selected_21$V1_30)

data_selected_22 <- data_selected %>%
  filter(!grepl(pattern = "^\\(", x = V1_30))

data_selected_3 <- full_join(data_selected_22, data_selected_21)

data_selected_31 <- data_selected_3 %>%
  filter(grepl(pattern = "^\\(", x = V1_33))

for (i in (ncol(data_selected_31):11))
{
  data_selected_31[, i] <- data_selected_31[, i - 1]
}

data_selected_31$V1_33 <- NA
data_selected_31$V1_33 <- as.character(data_selected_31$V1_33)

data_selected_32 <- data_selected_3 %>%
  filter(!grepl(pattern = "^\\(", x = V1_33))

data_manipulated <-
  full_join(data_selected_32, data_selected_31) %>%
  mutate(
    referenceAuthors_1 = paste(V1_29, V1_30, sep = " "),
    referenceAuthors_2 = paste(V1_32, V1_33, sep = " "),
    referenceAuthors_3 = paste(V1_35, V1_36, sep = " "),
    referenceAuthors_4 = paste(V1_38, V1_39, sep = " ")
  ) %>%
  select(
    uniqueid,
    name,
    smiles,
    biologicalsource,
    cas,
    pubchem,
    referenceAuthors_1,
    referenceAuthors_2,
    referenceAuthors_3,
    referenceAuthors_4,
    referenceTitle_1 = V1_31,
    referenceTitle_2 = V1_34,
    referenceTitle_3 = V1_37,
    referenceTitle_4 = V1_40
  ) %>%
  pivot_longer(
    7:ncol(.),
    names_to = c(".value", "level"),
    names_sep = "_",
    values_to = "reference",
    values_drop_na = TRUE
  ) %>%
  mutate(
    referenceAuthors = paste(
      gsub(
        pattern = "NA",
        replacement = "",
        x = referenceAuthors,
        fixed = TRUE
      ),
      str_extract(string = referenceTitle, pattern = "\\([0-9]{4}\\)")
    ),
    referenceTitle = gsub(
      pattern = "\\([0-9]{4}\\)",
      replacement = "",
      x = referenceTitle
    )
  )

data_manipulated$referenceAuthors <-
  gsub(
    pattern = " NA",
    replacement = "",
    x = data_manipulated$referenceAuthors
  )

data_manipulated$referenceAuthors <-
  trimws(data_manipulated$referenceAuthors)

data_manipulated$referenceAuthors <-
  y_as_na(x = data_manipulated$referenceAuthors, y = "")

data_manipulated$referenceTitle <-
  trimws(data_manipulated$referenceTitle)

data_manipulated <- data_manipulated %>% 
  select(uniqueid,
        name,
        smiles,
        biologicalsource,
        cas,
        pubchem,
        reference_authors = referenceAuthors,
        reference_title = referenceTitle) %>% 
  data.frame()

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_manipulated,
    db = "san_1",
    structure_field = c("name", "smiles"),
    reference_field = c("reference_authors", "reference_title")
  )

# exporting
database$writeInterim(data_standard)
