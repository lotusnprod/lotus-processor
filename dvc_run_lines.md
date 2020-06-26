
Following this example trying to run a dvc pipeline <https://dvc.org/doc/tutorials/get-started/data-pipelines#pipeline-stages>

dvc run -n prepare \
          -p prepare.seed,prepare.split \
          -d src/prepare.py -d data/data.xml \
          -o data/prepared \
          python src/prepare.py data/data.xml


Let's traduce this for the smiles.py file

dvc run -n smiler \
          -d src/2_curating/2_editing/chemo/subscripts/1_translating/smiles.py -d data/interim/tables/0_original/originalStructureSmiles.tsv.zip \
          -o data/interim/tables/1_translated/translatedStructureSmiles_dvc.tsv.zip \
          python src/2_curating/2_editing/chemo/subscripts/1_translating/smiles.py

This outputs a collision with dvc outputs 

ERROR: failed to run command - Paths for outs:                          
'data/interim'('data/interim.dvc')
'data/interim/tables/1_translated/translatedStructureSmiles_dvc.tsv.zip'('pipelines.yaml:smiler')
overlap. To avoid unpredictable behaviour, rerun command with non overlapping outs paths.

See here for details https://github.com/iterative/dvc/issues/2984


dvc run will systematically remove the -o directory so beware. Its a surprising and somewhat dangerous ? behavior.
See
https://github.com/iterative/dvc/issues/2027


dvc run -n smiler \
          -d src/2_curating/2_editing/chemo/subscripts/1_translating/smiles.py -d data/interim/tables/0_original/originalStructureSmiles.tsv.zip \
          -o data/dvc_outputs \
          -f data_smiled.dvc \
          python src/2_curating/2_editing/chemo/subscripts/1_translating/smiles.py


dvc run -n smiler \
          -d src/2_curating/2_editing/chemo/subscripts/1_translating/smiles_min.py -d data/interim/tables/0_original/originalStructureSmiles.tsv.zip \
          -o data/prepared/ \
          cd src && python 2_curating/2_editing/chemo/subscripts/1_translating/smiles_min.py


dvc run -n smiler \
          -d src/2_curating/2_editing/chemo/subscripts/1_translating/smiles_min.py -d data/interim/tables/0_original/originalStructureSmiles.tsv.zip \
          python src/2_curating/2_editing/chemo/subscripts/1_translating/smiles_min.py

dvc run -n smiler \
          -d src/2_curating/2_editing/chemo/subscripts/1_translating/smiles_min.py -d data/interim/tables/0_original/originalStructureSmiles.tsv.zip \
          -o data/out \
          python src/2_curating/2_editing/chemo/subscripts/1_translating/smiles_min.py
