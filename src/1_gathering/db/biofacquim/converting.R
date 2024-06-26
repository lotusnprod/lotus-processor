# title: "biofacquim converteR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

if (!require(BiocManager)) {
  install.packages("BiocManager")
  library(BiocManager)
}
if (!require(ChemmineR)) {
  BiocManager::install("BiocManager")
  library(ChemmineR)
}


# get paths
database <- databases$get("biofacquim")

data <- ChemmineR::read.SDFset(sdfstr = database$sourceFiles$sdf)

df <- data.frame()
df[1, 1] <- ""
df[1, 2] <- ""
df[1, 3] <- ""
df[1, 4] <- ""
df[1, 5] <- ""
df[1, 6] <- ""
df[1, 7] <- ""
df[1, 8] <- ""
df[1, 9] <- ""
df[1, 10] <- ""
colnames(df)[1] <- "ID"
colnames(df)[2] <- "Name"
colnames(df)[3] <- "SMILES"
colnames(df)[4] <- "Reference"
colnames(df)[5] <- "Year"
colnames(df)[6] <- "Genus"
colnames(df)[7] <- "Specie"
colnames(df)[8] <- "DOI"
colnames(df)[9] <- "Journal"
colnames(df)[10] <- "Site"

for (i in seq_along(data@SDF)) {
  df[i, 1] <- data@SDF[[i]]@datablock[["ID"]]
  df[i, 2] <- data@SDF[[i]]@datablock[["Name"]]
  df[i, 3] <- data@SDF[[i]]@datablock[["SMILES"]]
  df[i, 4] <- data@SDF[[i]]@datablock[["Reference"]]
  df[i, 5] <- data@SDF[[i]]@datablock[["Year"]]
  df[i, 6] <- data@SDF[[i]]@datablock[["Genus"]]
  df[i, 7] <- data@SDF[[i]]@datablock[["Specie"]]
  df[i, 8] <- data@SDF[[i]]@datablock[["DOI"]]
  df[i, 9] <- data@SDF[[i]]@datablock[["Journal"]]
  df[i, 10] <- data@SDF[[i]]@datablock[["Site"]]
}

# exporting
database$writeFile(database$sourceFiles$tsv, df)
