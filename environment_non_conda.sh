#!/usr/bin/env bash

echo "Installing packages missing from conda"

# R --slave -e "install.packages('rcrossref', repo = \"https://cran.rstudio.com/\")"

# R --slave -e "install.packages('readr', repo = \"https://cran.rstudio.com/\")"

# R --slave -e "install.packages('tidyjson', repo = \"https://cran.rstudio.com/\")"

# R --slave -e "install.packages('https://cran.r-project.org/src/contrib/Archive/classyfireR/classyfireR_0.3.6.tar.gz', repo = NULL, type = 'source')"

conda skeleton cran rcrossref
conda skeleton cran readr
conda skeleton cran tidyjson
conda skeleton cran classyfirer --allow-archived