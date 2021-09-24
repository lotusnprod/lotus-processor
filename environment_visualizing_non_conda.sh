#!/usr/bin/env bash

echo "Installing packages missing from conda"

# R --slave -e "install.packages('ggstar', repo = \"https://cran.rstudio.com/\")"

conda skeleton cran ggstar
