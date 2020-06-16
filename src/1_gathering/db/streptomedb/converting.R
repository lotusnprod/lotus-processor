# title: "STREPTOMEDB compileR"

# loading paths
source("paths.R")

# loading functions
source("functions.R")

data <- read.SDFset(pathDataExternalDbSourceStreptomedbOriginal)

df <- data.frame()
df[1, 1] <- ""
df[1, 2] <- ""
df[1, 3] <- ""
df[1, 4] <- ""
df[1, 5] <- ""
df[1, 6] <- ""
colnames(df)[1] <- "name"
colnames(df)[2] <- "smiles"
colnames(df)[3] <- "uniqueid"
colnames(df)[4] <- "pubchem"
colnames(df)[5] <- "biologicalsource"
colnames(df)[6] <- "pubmedid"

for (i in 1:length(data@SDF)) {
  df[i, 1] <- data@SDF[[i]]@header[["Molecule_Name"]]
  df[i, 2] <- data@SDF[[i]]@datablock[["canonical_smiles"]]
  df[i, 3] <- data@SDF[[i]]@datablock[["compound_id"]]
  df[i, 4] <- data@SDF[[i]]@datablock[["pubchem cid"]]
  df[i, 5] <- data@SDF[[i]]@datablock[["organisms"]]
  df[i, 6] <- data@SDF[[i]]@datablock[["pmids"]]
}

# exporting
write.table(
  x = df,
  file = gzfile(
    description = pathDataExternalDbSourceStreptomedbCompiled,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
