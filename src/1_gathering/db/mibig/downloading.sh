#!/usr/bin/env bash
# -*- coding: utf-8 -*-

if [ ! -d ../data/external/dbSource/mibig/ ]; then

mkdir -p ../data/external/dbSource/mibig/

wget -nv "https://dl.secondarymetabolites.org/mibig/mibig_json_3.1.tar.gz" -O ../data/external/dbSource/mibig/mibig_json_3.1.tar.gz

tar -xzf ../data/external/dbSource/mibig/mibig_json_3.1.tar.gz -C ../data/external/dbSource/mibig/

zip -r ../data/external/dbSource/mibig/mibig_json_3.1.zip ../data/external/dbSource/mibig/mibig_json_3.1

rm ../data/external/dbSource/mibig/mibig_json_3.1.tar.gz

rm -r ../data/external/dbSource/mibig/mibig_json_3.1

fi

echo "Done"