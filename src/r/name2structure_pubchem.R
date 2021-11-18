#' Title
#'
#' @param i
#'
#' @return
#' @export
#'
#' @examples
name2smiles_pubchem <- function(i) {
  tryCatch(
    {
      cpd <- dataForPubchem[i, "nameCleaned"]
      url <-
        paste0(
          "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/",
          cpd,
          "/property/IsomericSMILES/TXT"
        )
      url <- gsub(
        pattern = "\\s",
        replacement = "%20",
        x = url
      )
      read_html(url) %>%
        html_text() %>%
        gsub(
          pattern = "\n$",
          replacement = ""
        )
    },
    error = function(e) {
      NA
    }
  )
}
