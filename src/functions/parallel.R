numCores <- Sys.getenv("NUM_CORES", unset=(parallel::detectCores() - 1))
