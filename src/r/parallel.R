numCores <-
  if (.Platform$OS.type == "unix") {
    min(max(1, parallel::detectCores() - 1), 10)
  } else {
    1
  }
