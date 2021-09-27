#!/usr/bin/env bash
# -*- coding: utf-8 -*-

if [ ! -f ../data/external/dbSource/tmmc/compound.xlsx ]; then

mkdir -p ../data/external/dbSource/tmmc/

wget -nv "http://informatics.kiom.re.kr/compound/download/compound.xlsx" -O ../data/external/dbSource/tmmc/compound.xlsx

fi

echo "Done"