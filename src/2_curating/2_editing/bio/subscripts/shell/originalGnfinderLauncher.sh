for file in ../data/interim/tables_min/0_original/gnfinder/*.tsv

do
	echo "Submitting File $file to gnfinder servers"
	filename=$(basename $file .tsv | sed s/translated/sanitized/g) #file has extension, it return only filename, and here we add a sed line to chenge a given string in the filename
  gnfinder find -c -l detect -s $(seq -s, 1 1 182 | sed s/,\$//) -t "0" < $file > ../data/interim/tables/2_cleaned/gnfinder/original/json/$filename.json
	echo "The result is saved as $filename.json"
done


