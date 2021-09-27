#!/usr/bin/env bash
# -*- coding: utf-8 -*-

if [ ! -d ../data/external/dbSource/drduke/Duke-Source-CSV ]; then

mkdir -p ../data/external/dbSource/drduke/Duke-Source-CSV

wget -nv "https://data.nal.usda.gov/system/files/Duke-Source-CSV.zip" -O ../data/external/dbSource/drduke/Duke-Source-CSV.zip

unzip ../data/external/dbSource/drduke/Duke-Source-CSV.zip -d ../data/external/dbSource/drduke/Duke-Source-CSV

rm ../data/external/dbSource/drduke/Duke-Source-CSV.zip

fi

echo "Done"