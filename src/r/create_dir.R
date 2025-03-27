#' @title Create directory
#'
#' @noRd
#'
#' @param export TODO
#'
#' @return TODO
#'
#' @export
#'
#' @examples TODO
create_dir <- function(export) {
  if (
    grepl(
      pattern = "[[:alnum:]]\\.",
      x = export,
    )
  ) {
    ifelse(
      test = !dir.exists(paths = dirname(path = export)),
      yes = dir.create(path = dirname(path = export), recursive = TRUE),
      no = paste(
        dirname(path = export),
        "exists"
      )
    )
  } else {
    ifelse(
      test = !dir.exists(paths = export),
      yes = dir.create(path = export, recursive = TRUE),
      no = paste(
        export,
        "exists"
      )
    )
  }
}

#' @title Create directory and remove files
#'
#' @noRd
#'
#' @param export TODO
#'
#' @return TODO
#'
#' @export
#'
#' @examples TODO
create_dir_with_rm <- function(export) {
  ifelse(
    test = !dir.exists(export),
    yes = dir.create(export),
    no = file.remove(
      list.files(
        path = export,
        full.names = TRUE
      )
    ) &
      dir.create(
        export,
        showWarnings = FALSE
      )
  )
}
