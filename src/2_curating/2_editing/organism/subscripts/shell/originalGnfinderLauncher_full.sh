for file in ../data/interim/tables/0_original/organism/*.tsv

do
	echo "Submitting File $file to gnfinder servers"
	filename=$(basename $file .tsv | sed s/translated/sanitized/g) #file has extension, it return only filename, and here we add a sed line to chenge a given string in the filename
  gnfinder find -c -l detect -s 3,4,6,11,146,158,165,167,180 -t "0" < $file > ../data/interim/tables/2_cleaned/organism/original/$filename.json
	echo "The result is saved as $filename.json"
done


