#!/usr/bin/env bash

echo "Installing packages missing from conda"

R --slave -e "install.packages('rcrossref', repo=\"https://cran.rstudio.com/\")"

R --slave -e "install.packages('classyfireR', "https://cran.r-project.org/src/contrib/Archive/classyfireR/classyfireR_0.3.6.tar.gz", repo=NULL, type="source")"

R --slave -e "install.packages('ggstar', repo=\"https://cran.rstudio.com/\")"

R --slave -e "install.packages('ggtree', repo=\"https://bioconductor.org/packages/release/bioc/\")"

R --slave -e "install.packages('ggtreeExtra', repo=\"https://bioconductor.org/packages/release/bioc/\")"
