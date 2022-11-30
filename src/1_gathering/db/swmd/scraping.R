# title: "SWMD scrapeR"

# loading paths
source("paths.R")
source("r/y_as_na.R")

library(dplyr)
library(data.table)
library(splitstackshape) # provides cSplit
library(rvest) # provides read_html

# get paths
database <- databases$get("swmd")

ids <-
  gsub(
    ".mol",
    "",
    list.files(path = pathDataExternalDbSourceSwmdDirectory)
  )

url <- "http://www.swmd.co.in/search.php?No="

X <- ids

getswmd <- function(X) {
  tryCatch({
    cd_id <- X
    url_id <- paste0(url, cd_id)

    df0 <- read_html(url_id) %>%
      html_element(xpath = "body")

    df1 <- df0 %>%
      html_element(xpath = "div[1]/table[3]") %>%
      html_table(fill = TRUE) %>%
      select(1, 2)

    df2 <- df0 %>%
      html_element(xpath = "div[2]/table") %>%
      html_table(fill = TRUE)

    df3 <- df0 %>%
      html_element(xpath = "div[4]/div/table") %>%
      html_table(fill = TRUE)

    df4 <- rbind(df1, df2, df3)

    return(df4)
  })
}

SWMD <- invisible(
  lapply(
    FUN = getswmd,
    X = X
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
  pivot_wider(
    names_from = name,
    values_from = value
  ) %>%
  ungroup() %>%
  select(-level)

SWMD_3[] <- lapply(SWMD_3, function(x) {
  gsub("\r\n", " ", x)
})
SWMD_3[] <- lapply(SWMD_3, function(x) {
  gsub("\r", " ", x)
})
SWMD_3[] <- lapply(SWMD_3, function(x) {
  gsub("\n", " ", x)
})

# exporting
ifelse(
  test = !dir.exists(dirname(database$sourceFiles$tsv)),
  yes = dir.create(dirname(database$sourceFiles$tsv)),
  no = paste(dirname(database$sourceFiles$tsv), "exists")
)

database$writeFile(database$sourceFiles$tsv, SWMD_3)
