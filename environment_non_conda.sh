#!/usr/bin/env bash

echo "Installing packages missing from conda"

# R --slave -e "install.packages('tidyjson', repo = \"https://cran.rstudio.com/\")"

# R --slave -e "install.packages('https://cran.r-project.org/src/contrib/Archive/classyfireR/classyfireR_0.3.6.tar.gz', repo = NULL, type = 'source')"

conda run -n lotus_env skeleton cran tidyjson
conda run -n lotus_env skeleton cran classyfirer --allow-archived
conda run -n lotus_env skeleton cran https://github.com/ropensci/rcrossref