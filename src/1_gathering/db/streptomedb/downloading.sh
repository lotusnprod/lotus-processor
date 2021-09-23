#!/usr/bin/env bash
# -*- coding: utf-8 -*-

if [ ! -f ../data/external/dbSource/streptomedb/streptomedb.sdf ]; then

mkdir -p ../data/external/dbSource/streptomedb/

wget "http://132.230.56.4/streptomedb/get_download_file/" -O ../data/external/dbSource/streptomedb/streptomedb.sdf

fi

echo "Done"