#!/usr/bin/env bash
# -*- coding: utf-8 -*-

mkdir -p ../data/external/dbSource/foodb/

wget "https://foodb.ca/public/system/downloads/foodb_2020_4_7_csv.tar.gz" -P ../data/external/dbSource/foodb/

tar -xzf ../data/external/dbSource/foodb/foodb_2020_4_7_csv.tar.gz -C ../data/external/dbSource/foodb/

rm ../data/external/dbSource/foodb/foodb_2020_4_7_csv.tar.gz