gzip -d ../data/interim/tables_test/2_cleaned/organism/verify.tsv.gz -k

../bin/gnverifier ../data/interim/tables_test/2_cleaned/organism/verify.tsv -s 3,4,5,6,8,9,11,12,118,128,132,147,148,150,155,158,163,164,165,167,169,174,175,179,180,187 -j 200 -f compact >../data/interim/tables_test/2_cleaned/organism/verified.json

rm ../data/interim/tables_test/2_cleaned/organism/verify.tsv