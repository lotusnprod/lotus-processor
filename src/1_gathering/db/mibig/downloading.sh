#!/usr/bin/env bash
# -*- coding: utf-8 -*-

if [ ! -d ../data/external/dbSource/mibig/ ]; then

mkdir -p ../data/external/dbSource/mibig/

wget -nv "https://dl.secondarymetabolites.org/mibig/mibig_json_2.0.tar.gz" -O ../data/external/dbSource/mibig/mibig_json_2.0.tar.gz

tar -xzf ../data/external/dbSource/mibig/mibig_json_2.0.tar.gz -C ../data/external/dbSource/mibig/

zip -r ../data/external/dbSource/mibig/mibig_json_2.0.zip ../data/external/dbSource/mibig/mibig_json_2.0

rm ../data/external/dbSource/mibig/mibig_json_2.0.tar.gz

rm -r ../data/external/dbSource/mibig/mibig_json_2.0

fi

echo "Done"