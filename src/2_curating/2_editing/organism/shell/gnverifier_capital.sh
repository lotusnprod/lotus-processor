#!/bin/sh

# A simple script with a function to launch gnverifier

gnverifier_capital() {
  INPUT=$1
  OUTPUT=$2

  echo "Submitting $INPUT to gnverifier servers"

  gzip -d $INPUT -k

  ../bin/gnverifier ${INPUT%.*} -s 3,4,5,6,8,9,11,12,118,128,132,147,148,150,155,158,163,164,165,167,169,174,175,179,180,187 -j 200 -f compact -c >$OUTPUT

  rm ${INPUT%.*}

  echo "Result is saved in $OUTPUT"
}