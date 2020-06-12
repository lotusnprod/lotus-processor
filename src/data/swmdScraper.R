#title: "SWMD scrapeR"

##functions
source("../../functions.R")

outpath <- "0_initial_files/SWMD_scraped.tsv.zip"

ids <- gsub(".mol", "", list.files(path = "0_initial_files/Mol/"))

url <- 'http://www.swmd.co.in/search.php?No='

X <- ids

getswmd <- function(X)
{
  tryCatch({
    cd_id <- X
    url_id <- paste(url, cd_id, sep = "")
    
    df0 <- read_html(url_id) %>%
      html_node(xpath = "body")
    
    df1 <- df0 %>%
      html_node(xpath = "div[1]/table[3]") %>%
      html_table(fill = TRUE) %>% select(1, 2)
    
    df2 <- df0 %>%
      html_node(xpath = "div[2]/table") %>%
      html_table(fill = TRUE)
    
    df3 <- df0 %>%
      html_node(xpath = "div[4]/div/table") %>%
      html_table(fill = TRUE)
    
    df4 <- rbind(df1, df2, df3)
    
    return(df4)
  })
}

SWMD <- invisible(
  pbmclapply(
    FUN = getswmd,
    X = X,
    mc.preschedule = TRUE,
    mc.set.seed = TRUE,
    mc.silent = TRUE,
    mc.cores = (parallel::detectCores() - 2),
    mc.cleanup = TRUE,
    mc.allow.recursive = TRUE
  )
)

SWMD_2 <- bind_rows(SWMD)

SWMD_2$level <- as.numeric(gl(nrow(SWMD_2) / 21, 21))

colnames(SWMD_2) <- c("name", "value", "level")

SWMD_2$name <- y_as_na(SWMD_2$name, "")
SWMD_2$value <- y_as_na(SWMD_2$value, "")

SWMD_3 <- SWMD_2 %>%
  filter(!str_detect(name, "Accession Number\r\n")) %>%
  filter(!is.na(name)) %>%
  group_by(level) %>%
  pivot_wider(names_from = name,
              values_from = value) %>%
  ungroup() %>%
  select(-level)

SWMD_3[] <- lapply(SWMD_3, function(x)
  gsub("\r\n", " ", x))
SWMD_3[] <- lapply(SWMD_3, function(x)
  gsub("\r", " ", x))
SWMD_3[] <- lapply(SWMD_3, function(x)
  gsub("\n", " ", x))

#exporting
write.table(
  x = SWMD_3,
  file = gzfile(description = outpath,
                compression = 9,
                encoding = "UTF-8"),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
