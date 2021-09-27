#!/usr/bin/env bash
# -*- coding: utf-8 -*-

if [ ! -f ../data/external/dbSource/tppt/TPPT_database.xlsx ]; then

mkdir -p ../data/external/dbSource/tppt/

wget -nv "https://www.agroscope.admin.ch/dam/agroscope/de/dokumente/publikationen/tppt-xls.xlsx.download.xlsx/TPPT_database.xlsx" -O ../data/external/dbSource/tppt/TPPT_database.xlsx

fi

echo "Done"