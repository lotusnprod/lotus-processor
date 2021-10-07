#!/bin/bash

timestamps="../data/interim/timestamps.txt"
> "$timestamps"

for filename in ../data/external/dbSource/**/*; do
    echo "$filename" >> "$timestamps"
    date -r "$filename" -u +"%Y-%m-%dT%H:%M:%SZ" >> "$timestamps"
done
