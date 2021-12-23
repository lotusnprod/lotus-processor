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
Rscript 1_gathering/db/custom/standardizing.R ../data/example/custom_example.tsv
cd ..
make MODE=custom lotus-bloom
make MODE=test lotus-bloom
make MODE=test lotus-check