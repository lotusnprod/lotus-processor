# title: "Knapsack scrapeR"

# loading paths
source("paths.R")

library(dplyr)
library(data.table)
library(readr)
library(splitstackshape) # provides cSplit
library(stringr) # provides str_pad
library(rvest) # provides read_html
library(xml2)

# get paths
database <- databases$get("knapsack")

url <-
  "http://www.knapsackfamily.com/knapsack_core/information.jsp?word=C00"

X <- 1:63574

GetKnapSackTable <- function(X) {
  tryCatch(
    {
      cd_id <- stringr::str_pad(string = X, width = 6, pad = "0")
      url_id <- paste(url, cd_id) |>
        gsub(pattern = "\\s", replacement = "")
      df0 <- rvest::read_html(x = url_id) |>
        rvest::html_element(xpath = "/html/body/div/div[2]/div[2]/table") |>
        xml2::xml_child(2) |>
        xml2::xml_child(1) |>
        xml2::xml_child(1)

      df1 <- df0 |>
        rvest::html_table(fill = TRUE) |>
        head(8) |>
        dplyr::select(1, 2) |>
        t()

      colnames(df1) <- df1[1, ]

      df1 <- data.frame(df1) %>%
        dplyr::filter(rownames(.) == "X2") %>%
        dplyr::mutate(joincol = url_id)

      df2 <- df0 |>
        xml2::xml_child(10) |>
        xml2::xml_child(2) |>
        xml2::xml_child(1) |>
        rvest::html_table(fill = TRUE)

      i <- 0:(nrow(df2) - 1) %>%
        unlist() %>%
        data.frame() %>%
        dplyr::mutate(joincol = url_id)

      df_for_ref <- dplyr::full_join(df1, i) %>%
        dplyr::mutate(link = paste0(joincol, "&key=", .)) %>%
        dplyr::select(-joincol, -.)

      return(df_for_ref)
    },
    error = function(e) {
      NA
    }
  )
}

GetKnapSackRef <- function(X) {
  tryCatch(
    {
      df3 <- rvest::read_html(X) |>
        rvest::html_element(xpath = "/html/body/div/div[2]/div[2]/table") |>
        xml2::xml_child(2) |>
        xml2::xml_child(2) |>
        xml2::xml_child(7) |>
        rvest::html_table(fill = TRUE) |>
        t()

      colnames(df3) <- df3[1, ]

      df3 <- data.frame(df3) %>%
        dplyr::filter(rownames(.) == "X2")

      df4 <- cbind(c("link" = X, df3)) |>
        t() |>
        data.frame() |>
        dplyr::mutate_all(as.character)

      return(df4)
    },
    error = function(e) {
      NA
    }
  )
}

invisible(lapply(
  FUN = GetKnapSackTable,
  X = head(X)
))
df1 <- invisible(pbmcapply::pbmclapply(
  FUN = GetKnapSackTable,
  X = X,
  mc.cores = parallel::detectCores(logical = TRUE)
))

KnapSackTable <- dplyr::bind_rows(df1[!is.na(df1)])

X <- KnapSackTable$link

invisible(lapply(
  FUN = GetKnapSackRef,
  X = head(X)
))
df3 <- invisible(pbmcapply::pbmclapply(
  FUN = GetKnapSackRef,
  X = X,
  mc.cores = parallel::detectCores(logical = TRUE)
))

df4 <- bind_rows(df3[!is.na(df3)])

KNApSAcK_db <- dplyr::full_join(KnapSackTable, df4)

# exporting
ifelse(
  test = !dir.exists(dirname(database$sourceFiles$tsv)),
  yes = dir.create(dirname(database$sourceFiles$tsv)),
  no = paste(dirname(database$sourceFiles$tsv), "exists")
)

database$writeFile(database$sourceFiles$tsv, KNApSAcK_db)
