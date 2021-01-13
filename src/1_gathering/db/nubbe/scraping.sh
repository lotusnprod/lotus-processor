curl 'https://nubbe.iq.unesp.br/portal/do/Query' --compressed -H 'Content-Type: application/x-www-form-urlencoded' -H 'X-Requested-With: XMLHttpRequest' -H 'Origin: https://nubbe.iq.unesp.br' -H 'Referer: https://nubbe.iq.unesp.br/portal/nubbe-search.html' --data-raw 'service=17&tipo_1=' > ../data/external/dbSource/nubbe/nubbe_ids.xml
echo "" > ../data/external/dbSource/nubbe/folder/download.txt
for i in $(grep -o '<id>.*</id>' ../data/external/dbSource/nubbe/nubbe_ids.xml | sed 's/\(<id>\|<\/id>\)//g' | sort -u); do

    echo "url=\"https://nubbe.iq.unesp.br/portal/do/Query?service=21&id=$i\"" >> ../data/external/dbSource/nubbe/folder/download.txt
    echo "output=\"$i.xml\"" >> ../data/external/dbSource/nubbe/folder/download.txt
done
curl -K ../data/external/dbSource/nubbe/folder/download.txt
rm ../data/external/dbSource/nubbe/folder/download.txt