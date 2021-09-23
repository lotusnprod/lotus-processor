#!/usr/bin/env bash
# -*- coding: utf-8 -*-

if [ ! -f ../data/external/dbSource/biofacquim/BIOFACQUIM_V2.sdf ]; then

mkdir -p ../data/external/dbSource/biofacquim/

wget "https://ndownloader.figshare.com/files/20050244" -O ../data/external/dbSource/biofacquim/BIOFACQUIM_V2.sdf

fi

echo "Done"