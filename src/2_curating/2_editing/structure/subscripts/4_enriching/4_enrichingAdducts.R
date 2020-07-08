#title: "DB adducteR"

#loading functions
source("functions.R")

#writing paths
##inpath
inpath <-
  paste("../Outputs/tables/2_sanitized/sanitizedStructure.tsv")

##outpath pos
outpathPos <- paste("../Outputs/tables/DBAdductedPos.tsv.zip")

##outpath neg
outpathNeg <- paste("../Outputs/tables/DBAdductedNeg.tsv.zip")

##loading db
db <- read_delim(
  file = inpath,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  select(structureTranslated,
         validatorLog,
         smiles = smiles_sanitized,
         inchi = inchi_sanitized,
         inchikey = inchikeySanitized,
         inchikey2D = shortikSanitized,
         molecularFormula = formulaSanitized,
         exactMass = exactmassSanitized,
         xlogP = xlogpSanitized
  ) %>% 
  filter(!is.na(molecularFormula)) %>%
  distinct(molecularFormula, .keep_all = TRUE)

#writing exact masses of adducts
proton <- 1.0072765

#pos
ammonium  <- 17.026547
sodium <- 22.98922
methanol <- 32.02621
potassium <- 38.96316
acetonitrile <- 41.02655
ethylamine <- 46.06513
isopropanol <- 61.06534
dmso <- 79.02122

#neg
H2O <- 18.01056
chlorine <- 34.96940
formic <- 46.00548
acetic <- 60.02113
bromine <- 78.91834
TFA <- 113.99286

#adducting
##pos
dataAdductedPos <- db %>%
  mutate(Pos3H = (exactMass + 3 * proton) / 3) %>%
  mutate(Pos2HNa = (exactMass + 2 * proton + sodium) / 3) %>%
  mutate(PosH2Na = (exactMass + proton + 2 * sodium) / 3) %>%
  mutate(Pos3Na = (exactMass + 3 * sodium) / 3) %>%
  mutate(Pos2H = ((exactMass + 2 * proton) / 2)) %>%
  mutate(PosHNH3 = ((exactMass + 2 * proton + ammonium) / 2)) %>%
  mutate(PosHNa = ((exactMass + proton + sodium) / 2)) %>%
  mutate(PosHK = ((exactMass + proton + potassium) / 2)) %>%
  mutate(Pos2HCH3CN = ((exactMass + 2 * proton + acetonitrile) / 2)) %>%
  mutate(Pos2Na = ((exactMass + 2 * sodium) / 2)) %>%
  mutate(Pos2H2CH3CN = ((exactMass + 2 * proton + 2 * acetonitrile) / 2)) %>%
  mutate(Pos2H3CH3CN = ((exactMass + 2 * proton + 3 * acetonitrile) / 2)) %>%
  mutate(PosH = exactMass + proton) %>%
  mutate(PosHNH3 = exactMass + proton + ammonium) %>%
  mutate(PosNa = exactMass + sodium) %>%
  mutate(PosHCH3OH = exactMass + proton + methanol) %>%
  mutate(PosK = exactMass + potassium) %>%
  mutate(PosHCH3CN = exactMass + proton + acetonitrile) %>%
  mutate(Pos2NaH = exactMass - proton + 2 * sodium) %>%
  mutate(PosHC2H7N = exactMass + proton + ethylamine) %>%
  mutate(PosHIsoP = exactMass + proton + isopropanol) %>%
  mutate(PosCH3CNNa = exactMass + acetonitrile + sodium) %>%
  mutate(Pos2KH = exactMass - proton + 2 * potassium) %>%
  mutate(PosHDMSO = exactMass + proton + dmso) %>%
  mutate(PosH2ACN = exactMass + proton + 2 * acetonitrile) %>%
  #mutate(PosIsoPNa-H = exactMass - proton + isopropanol + sodium)
  mutate(Pos2MH = 2 * exactMass + proton) %>%
  mutate(Pos2MHNH3 = 2 * exactMass + proton + ammonium) %>%
  mutate(Pos2MNa = 2 * exactMass + sodium) %>%
  mutate(Pos2MK = 2 * exactMass + potassium) %>%
  mutate(Pos2MHCH3CN = 2 * exactMass + proton + acetonitrile) %>%
  mutate(Pos2MCH3CNNa = 2 * exactMass + acetonitrile + sodium)

##neg
dataAdductedNeg <- db %>%
  mutate(Neg3H = (exactMass - 3 * proton) / 3) %>%
  mutate(Neg2H = ((exactMass - 2 * proton) / 2)) %>%
  #mutate(adductNegH2OH = exactMass - H2O - proton) %>%
  mutate(NegH = exactMass - proton) %>%
  mutate(NegNa2H = exactMass + sodium - 2 * proton) %>%
  mutate(NegCl = exactMass + chlorine) %>%
  mutate(NegK2H = exactMass + potassium - 2 * proton) %>%
  mutate(NegFAH = exactMass + formic - proton) %>%
  mutate(NegACH = exactMass + acetic - proton) %>%
  mutate(NegFANa2H = exactMass + formic + sodium - 2 * proton) %>%
  mutate(NegBr = exactMass + bromine) %>%
  mutate(NegTFAH = exactMass + TFA - proton) %>%
  mutate(Neg2MH = 2 * exactMass - proton) %>%
  mutate(Neg2MFAH = 2 * exactMass + formic - proton) %>%
  mutate(Neg2MACH = 2 * exactMass + acetic - proton) %>%
  mutate(Neg3MH = 3 * exactMass - proton)

#changing colnames
colnames(dataAdductedPos)[(ncol(db)+1):ncol(dataAdductedPos)] <-
  paste("adduct", colnames(dataAdductedPos)[(ncol(db)+1):ncol(dataAdductedPos)], sep = "_")

colnames(dataAdductedNeg)[(ncol(db)+1):ncol(dataAdductedNeg)] <-
  paste("adduct", colnames(dataAdductedNeg)[(ncol(db)+1):ncol(dataAdductedNeg)], sep = "_")

#pivoting long
##pos
dataAdductedLongPos <- dataAdductedPos %>%
  pivot_longer((ncol(db)+1):ncol(.)) %>%
  select(everything(),
         adduct = name,
         adductMass = value)

##neg
dataAdductedLongNeg <- dataAdductedNeg %>%
  pivot_longer((ncol(db)+1):ncol(.)) %>%
  select(everything(),
         adduct = name,
         adductMass = value)

#exporting
##pos
write.table(
  x = dataAdductedLongPos,
  file = gzfile(
    description = outpathPos,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

##neg
write.table(
  x = dataAdductedLongNeg,
  file = gzfile(
    description = outpathNeg,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
