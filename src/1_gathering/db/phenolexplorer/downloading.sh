#!/usr/bin/env bash
# -*- coding: utf-8 -*-

if [ ! -d ../data/external/dbSource/phenolexplorer/ ]; then

mkdir -p ../data/external/dbSource/phenolexplorer/

wget "http://phenol-explorer.eu/system/downloads/current/composition-data.xlsx.zip" -O ../data/external/dbSource/phenolexplorer/composition-data.xlsx.zip

unzip ../data/external/dbSource/phenolexplorer/composition-data.xlsx.zip -d ../data/external/dbSource/phenolexplorer/

rm ../data/external/dbSource/phenolexplorer/composition-data.xlsx.zip

wget "http://phenol-explorer.eu/system/downloads/current/compounds.csv.zip" -O ../data/external/dbSource/phenolexplorer/compounds.csv.zip

unzip ../data/external/dbSource/phenolexplorer/compounds.csv.zip -d ../data/external/dbSource/phenolexplorer/

rm ../data/external/dbSource/phenolexplorer/compounds.csv.zip

wget "http://phenol-explorer.eu/system/downloads/current/compounds-structures.csv.zip" -O ../data/external/dbSource/phenolexplorer/compounds-structures.csv.zip

unzip ../data/external/dbSource/phenolexplorer/compounds-structures.csv.zip -d ../data/external/dbSource/phenolexplorer/

rm ../data/external/dbSource/phenolexplorer/compounds-structures.csv.zip

wget "http://phenol-explorer.eu/system/downloads/current/foods.csv.zip" -O ../data/external/dbSource/phenolexplorer/foods.csv.zip

unzip ../data/external/dbSource/phenolexplorer/foods.csv.zip -d ../data/external/dbSource/phenolexplorer/

rm ../data/external/dbSource/phenolexplorer/foods.csv.zip

wget "http://phenol-explorer.eu/system/downloads/current/publications.csv.zip" -O ../data/external/dbSource/phenolexplorer/publications.csv.zip

unzip ../data/external/dbSource/phenolexplorer/publications.csv.zip -d ../data/external/dbSource/phenolexplorer/

rm ../data/external/dbSource/phenolexplorer/publications.csv.zip

fi

echo "Done"