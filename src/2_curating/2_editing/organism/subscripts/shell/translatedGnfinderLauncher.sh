for file in ../data/interim/tables/1_translated/organism/*.tsv

do
	echo "Submitting File $file to gnfinder servers"
	f=$(echo "${file##*/}")
  	filename=$(echo $f| cut  -d'.' -f 1| sed  s/translated/translated/g) #file has extension, it return only filename, and here we add a sed line to chenge a given string in the filename
	cat $file | gnfinder find -c -l detect -s $(seq -s, 1 1 182 | sed s/,\$//) -t "0" > ../data/interim/tables/2_cleaned/organism/translated/$filename.json
	echo "The result is saved as $filename.json"
done


