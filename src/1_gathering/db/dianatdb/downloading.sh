#!/usr/bin/env bash
# -*- coding: utf-8 -*-

if [ ! -f ../data/external/dbSource/dianatdb/2020_DiaNatDB_336.xlsx ]; then

mkdir -p ../data/external/dbSource/dianatdb/

wget "http://rdu.iquimica.unam.mx/bitstream/20.500.12214/1186/3/2020_DiaNatDB_336.xlsx" -O ../data/external/dbSource/dianatdb/2020_DiaNatDB_336.xlsx

fi

echo "Done"