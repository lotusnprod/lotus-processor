#!/usr/bin/env bash

cd src
Rscript 1_gathering/db/custom/standardizing.R ../data/example/custom_example.tsv
cd ..
make MODE=custom lotus-bloom
make MODE=test lotus-bloom
make MODE=test lotus-check