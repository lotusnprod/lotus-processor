NPCLASSIFIER_VERSION="1.5"
INDEX_VERSION="1"

if [ ! -f ../data/external/taxonomySource/structure/index_v"$INDEX_VERSION".json ]; then
  echo "Downloading"
  curl -o ../data/external/taxonomySource/structure/index_v"$INDEX_VERSION".json https://raw.githubusercontent.com/mwang87/NP-Classifier/master/Classifier/dict/index_v"$INDEX_VERSION".json
fi

echo "You can now import in sqlite: "
echo " it is a json"
echo " next steps to do"