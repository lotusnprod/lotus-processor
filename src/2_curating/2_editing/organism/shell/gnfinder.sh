#!/usr/bin/env bash
# -*- coding: utf-8 -*-

# A simple script with a function to launch gnfinder

gnfinder() {
  INPUT=$1
  OUTPUT=$2

  for file in $INPUT/*tsv; do

    echo "Submitting $INPUT to gnfinder servers"

    filename=$(basename $file .tsv | sed s/translated/sanitized/g) # file has extension, it return only filename, and here we add a sed line to chenge a given string in the filename

    ../bin/gnfinder -f "pretty" -s 3,4,5,6,8,9,11,12,118,128,132,147,148,150,155,158,163,164,165,167,169,174,175,180,187 <$file >$OUTPUT/$filename.json

    echo "Result is saved in $OUTPUT"
  done
}
