#!/usr/bin/env bash

echo "Installing packages missing from conda"

R --slave -e "install.packages('rcrossref', repo=\"https://cran.rstudio.com/\")"