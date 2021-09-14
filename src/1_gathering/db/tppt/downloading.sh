#!/usr/bin/env bash
# -*- coding: utf-8 -*-

mkdir -p ../data/external/dbSource/tppt/

wget "https://www.agroscope.admin.ch/dam/agroscope/de/dokumente/publikationen/tppt-xls.xlsx.download.xlsx/TPPT_database.xlsx" -O ../data/external/dbSource/tppt/TPPT_database.xlsx
