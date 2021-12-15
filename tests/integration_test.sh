#!/usr/bin/env bash

Rscript -e "
if (!require(remotes)) {
  install.packages(remotes)
  require(
    package = remotes,
    quietly = TRUE
  )
}
remotes::install_github('ropensci/rcrossref')
"
cd src
Rscript 1_gathering/db/manual/standardizing.R ../data/example/manual_example.tsv
cd ..
make MODE=manual manual-entry
make MODE=test lotus-quick
make MODE=test tests