# These are helpers functions and classes to handle the databases

#' Store a database files paths
#'
#' @field name Name of the database
#' @field sourceFile The source file of the database
#' @field interimFile The transformed file of the database
Database <-
  setRefClass(
    "Database",
    fields = list(
      name = "character",
      sourceFiles = "list",
      interimFile = "character"
    ),
    methods = list(
      writeFile = function(file, data) {
        write.table(
          x = data,
          file = gzfile(
            description = file,
            compression = 9,
            encoding = "UTF-8"
          ),
          row.names = FALSE,
          quote = FALSE,
          sep = "\t",
          fileEncoding = "UTF-8"
        )
      },
      writeInterim = function(data) {
        writeFile(interimFile, data)
      }
    )
  )

#' A container for databases
#'
#' @field dbSource path for the sources, any added database will use this path as starting path
#' @field dbInterim path for the interims, any added database will use this path as starting path
#' @field paths this is used internally to store the databases
#'
#' The database name will always be added to every source path as well (not interim yet)
Databases <-
  setRefClass(
    "Databases",
    fields = list(
      pathDbSource = "character",
      pathDbInterim = "character",
      paths = "environment"
    ),
    methods = list(
      initialize = function(pathDbSource, pathDbInterim) {
        .self$pathDbSource <- pathDbSource
        .self$pathDbInterim <- pathDbInterim
        .self$paths <- new.env(hash = T)
      },
      #' @param name Name of the database
      #' @param sourceFile Name of the source file (without path)
      #' @param interimFile Name of the interim file (without path)
      add = function(name, sourceFiles, interimFile) {
        .self$paths[[name]] <- Database$new(
          name = name,
          sourceFiles = lapply(
            sourceFiles,
            function(path) {
              file.path(pathDbSource, file.path(name, path))
            }
          ),
          interimFile = file.path(pathDbInterim, interimFile)
        )
      },
      get = function(name) {
        return(.self$paths[[name]])
      }
    )
  )