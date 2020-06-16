#title: "Knapsack scrapeR"

# setting working directory
setwd("~/GitLab/opennaturalproductsdb/src/")

# loading paths
source("paths.R")

# loading functions
source("functions.R")
  
url <-
  'http://www.knapsackfamily.com/knapsack_core/information.jsp?word=C00'

X <- 1:52045

GetKnapSackTable <- function(X)
{
  tryCatch({
    cd_id <- str_pad (X, 6, pad = "0")
    url_id <- paste(url, cd_id)
    url_id <- gsub("\\s", "", url_id)
    
    df0 <- read_html(url_id) %>%
      html_node(xpath = "/html/body/div/div[2]/div[2]/table") %>%
      xml_child(2) %>%
      xml_child(1) %>%
      xml_child(1)
    
    df1 <- df0 %>%
      html_table(fill = TRUE) %>%
      head(8) %>%
      select(1, 2) %>%
      t()
    
    colnames(df1) <- df1[1,]
    
    df1 <- data.frame(df1) %>%
      filter(rownames(.) == "X2") %>%
      mutate(joincol = url_id)
    
    df2 <- df0 %>%
      xml_child(10) %>%
      xml_child(2) %>%
      xml_child(1) %>%
      html_table(fill = TRUE)
    
    i <- 0:(nrow(df2) - 1) %>%
      unlist() %>%
      data.frame() %>%
      mutate(joincol = url_id)
    
    df_for_ref <- full_join(df1, i) %>%
      mutate(link = paste(joincol,
                          "&key=",
                          .,
                          sep = "")) %>%
      select(-joincol, -.)
    
    return(df_for_ref)
  },
  error = function(e) {
    NA
  })
}

GetKnapSackRef <- function(X)
{
  tryCatch({
    df3 <- read_html(X) %>%
      html_node(xpath = "/html/body/div/div[2]/div[2]/table") %>%
      xml_child(2) %>%
      xml_child(2) %>%
      xml_child(7) %>%
      html_table(fill = TRUE) %>%
      t()
    colnames(df3) <- df3[1, ]
    df3 <- data.frame(df3) %>%
      filter(rownames(.) == "X2")
    return(df3)
  },
  error = function(e) {
    NA
  })
}

df1 <- invisible(
  pbmclapply(
    FUN = GetKnapSackTable,
    X = X,
    mc.preschedule = FALSE,
    mc.set.seed = TRUE,
    mc.silent = TRUE,
    mc.cores = (parallel::detectCores()-2),
    mc.cleanup = FALSE,
    mc.allow.recursive = FALSE
  )
) 

KnapSackTable <- bind_rows(df1[!is.na(df1)])

X <- KnapSackTable$link

df3 <- invisible(
  pbmclapply(
    FUN = GetKnapSackRef,
    X = X,
    mc.preschedule = TRUE,
    mc.set.seed = TRUE,
    mc.silent = TRUE,
    mc.cores = (parallel::detectCores() - 2),
    mc.cleanup = TRUE,
    mc.allow.recursive = TRUE
  )
)

df4 <- bind_rows(df3[!is.na(df3)])

KNApSAcK_db <- cbind(KnapSackTable,df4)

# exporting
write.table(
  x = data_standard,
  file = gzfile(description = pathDataExternalDbSourceKnapsackOriginal,
                compression = 9,
                encoding = "UTF-8"),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
