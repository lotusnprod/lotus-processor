#!/usr/bin/env bash
export DB_DIR=./data/external/dbSource/nubbe
export DB_DIR_OUT="$DB_DIR"/folder
mkdir -p "$DB_DIR" "$DB_DIR_OUT"

curl 'https://nubbe.iq.unesp.br/portal/do/Query' --compressed -H 'Content-Type: application/x-www-form-urlencoded' -H 'X-Requested-With: XMLHttpRequest' -H 'Origin: https://nubbe.iq.unesp.br' -H 'Referer: https://nubbe.iq.unesp.br/portal/nubbe-search.html' --data-raw 'service=17&tipo_1=' > "$DB_DIR"/nubbe_ids.xml
rm "$DB_DIR_OUT"/download.txt
for i in $(grep -o '<id>.*</id>' "$DB_DIR"/nubbe_ids.xml | sed 's/\(<id>\|<\/id>\)//g' | sort -u); do
    echo "url=\"https://nubbe.iq.unesp.br/portal/do/Query?service=21&id=$i\"" >> "$DB_DIR_OUT"/download.txt
    echo "output=\"$i.xml\"" >> "$DB_DIR_OUT"/download.txt
done

cd "$DB_DIR_OUT"
curl -K "$DB_DIR_OUT"/download.txt
rm "$DB_DIR_OUT"/download.txt