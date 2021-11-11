#!/usr/bin/env bash

Rscript -e "remotes::install_github('ropensci/rcrossref')" &&
make MODE=test lotus-quick &&
cd src &&
Rscript 1_gathering/db/manual/standardizing.R ../data/example/manual_example.tsv &&
cd .. &&
make MODE=manual manual-entry
