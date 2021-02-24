#!/usr/bin/env bash

echo "Installing packages missing from conda"

R --slave -e "install.packages('rcrossref', repo=\"https://cran.rstudio.com/\")"

R --slave -e "install.packages('classyfireR', repo=\"https://cran.rstudio.com/\")"

R --slave -e "install.packages('ggstar', repo=\"https://cran.rstudio.com/\")"

R --slave -e "install.packages('ggtree', repo=\"https://bioconductor.org/\")"

R --slave -e "install.packages('ggtreeExtra', repo=\"https://bioconductor.org/\")"

R --slave -e "install.packages('treeio', repo=\"https://bioconductor.org/\")"