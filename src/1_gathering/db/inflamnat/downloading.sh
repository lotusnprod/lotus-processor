#!/usr/bin/env bash
# -*- coding: utf-8 -*-

# if [ ! -f ../data/external/dbSource/inflamnat/ci8b00560_si_001.xlsx ]; then
if [ ! -f ../data/external/dbSource/inflamnat/13321_2022_608_MOESM2_ESM.xlsx ]; then

mkdir -p ../data/external/dbSource/inflamnat/

## need ACS access... TODO
# wget "https://pubs.acs.org/doi/suppl/10.1021/acs.jcim.8b00560/suppl_file/ci8b00560_si_001.xlsx" -O ../data/external/dbSource/inflamnat/ci8b00560_si_001.xlsx
# wget "https://static-content.springer.com/esm/art%3A10.1186%2Fs13321-022-00608-5/MediaObjects/13321_2022_608_MOESM2_ESM.xlsx" -O ../data/external/dbSource/inflamnat/13321_2022_608_MOESM2_ESM.xlsx

fi

echo "Done"