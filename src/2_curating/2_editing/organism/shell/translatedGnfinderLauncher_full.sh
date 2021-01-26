for file in ../data/interim/tables/1_translated/organism/*.tsv; do
  echo "Submitting File $file to gnfinder servers"
  f=$(echo "${file##*/}")
  filename=$(echo $f | cut -d'.' -f 1 | sed s/translated/translated/g) # file has extension, it return only filename, and here we add a sed line to chenge a given string in the filename
  cat $file | gnfinder find -c -l detect -s 3,4,5,6,8,9,11,12,118,128,132,147,148,150,155,158,163,164,165,167,169,174,175,180,187 -t "0" >../data/interim/tables/2_cleaned/organism/translated/$filename.json
  echo "The result is saved as $filename.json"
done
