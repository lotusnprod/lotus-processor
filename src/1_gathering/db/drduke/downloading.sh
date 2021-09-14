#!/usr/bin/env bash
# -*- coding: utf-8 -*-

mkdir -p ../data/external/dbSource/drduke/Duke-Source-CSV

wget "https://data.nal.usda.gov/system/files/Duke-Source-CSV.zip" -O ../data/external/dbSource/drduke/Duke-Source-CSV.zip

tar -xf ../data/external/dbSource/drduke/Duke-Source-CSV.zip -C ../data/external/dbSource/drduke/Duke-Source-CSV

rm ../data/external/dbSource/drduke/Duke-Source-CSV.zip