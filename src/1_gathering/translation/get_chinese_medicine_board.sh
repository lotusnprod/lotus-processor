#!/usr/bin/env bash
# -*- coding: utf-8 -*-

if [ ! -f ../data/external/translationSource/tcm/Chinese-Medicine-Board---List---Nomenclature-compendium-of-commonly-used-Chinese-herbal-medicines.XLSX ]; then

mkdir -p ../data/external/translationSource/tcm/

wget --no-check-certificate "https://www.chinesemedicineboard.gov.au/documents/default.aspx?record=WD15%2f18746&dbid=AP&chksum=Cs%2baCFhYVrzbzL2%2bYnhtRA%3d%3d" -O ../data/external/translationSource/tcm/Chinese-Medicine-Board---List---Nomenclature-compendium-of-commonly-used-Chinese-herbal-medicines.XLSX

fi

echo "Done"