if [ ! -f ../data/external/taxonomySource/structure/classyfire/tax_nodes.json ]; then
  echo "Downloading"
  curl -o ../data/external/taxonomySource/structure/classyfire/tax_nodes.json http://classyfire.wishartlab.com/tax_nodes.json
fi

echo "You can now import in sql: "
echo " it is a json"
echo " next steps to do"
