#!/usr/bin/env bash
# -*- coding: utf-8 -*-

if [ ! -f ../data/external/dbSource/npcare/npcare.zip ]; then

mkdir -p ../data/external/dbSource/npcare/

wget "http://silver.sejong.ac.kr/npcare/csv/npcare.zip" -O ../data/external/dbSource/npcare/npcare.zip

fi

echo "Done"

