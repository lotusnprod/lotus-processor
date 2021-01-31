#!/usr/bin/env bash

echo "Installing packages missing from conda"

R --slave -e "install.packages('rcrossref', repo=\"https://cran.rstudio.com/\")"

R --slave -e "install.packages('classyfireR', repo=\"https://cran.rstudio.com/\")"