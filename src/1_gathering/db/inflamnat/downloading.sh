#!/usr/bin/env bash
# -*- coding: utf-8 -*-

if [ ! -f ../data/external/dbSource/inflamnat/ci8b00560_si_001.xlsx ]; then

mkdir -p ../data/external/dbSource/inflamnat/

## need ACS access... TODO
# wget "https://pubs.acs.org/doi/suppl/10.1021/acs.jcim.8b00560/suppl_file/ci8b00560_si_001.xlsx" -O ../data/external/dbSource/inflamnat/ci8b00560_si_001.xlsx

fi

echo "Done"