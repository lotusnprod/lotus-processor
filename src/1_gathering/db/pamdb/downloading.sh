#!/usr/bin/env bash
# -*- coding: utf-8 -*-

if [ ! -f ../data/external/dbSource/pamdb/PaMet.xlsx ]; then

mkdir -p ../data/external/dbSource/pamdb/

wget "http://pseudomonas.umaryland.edu/PaDl/PaMet.xlsx" -O ../data/external/dbSource/pamdb/PaMet.xlsx

fi

echo "Done"