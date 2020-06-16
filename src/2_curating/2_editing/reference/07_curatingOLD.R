#title: "Reference curatoR"

#writing paths
##inputs
inpath <- "outputs/tables/2_sanitized/sanitizedReference.tsv.zip"

##outputs
outpath <- "outputs/tables/3_curated/curatedReference.tsv.zip"

#loading files
dataSanitized <- read_delim(
  file = gzfile(inpath),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)


###Find something appropriate


dataCurated <- dataSanitized %>%
  mutate(
    curatedTitle = sanitizedTitle,
    curatedJournal = sanitizedJournal,
    curatedDoi = sanitizedDoi,
    curatedAuthor = sanitizedAuthor,
    curatedDate = sanitizedDate,
    curatedTranslationScore = sanitizedTranslationScore,
  )

#export
write.table(
  x = dataCurated,
  file = gzfile(
    description = outpath,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
