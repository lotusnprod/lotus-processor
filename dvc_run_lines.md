
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


Works now !

so the "dvc traduction of 

curating-editing-chemo-smiles: ${INTERIM_TABLE_TRANSLATED_PATH}/translatedStructureSmiles_min.tsv.zip
${INTERIM_TABLE_TRANSLATED_PATH}/translatedStructureSmiles_min.tsv.zip: ${SRC_CURATING_EDITING_CHEMO_SUBSCRIPTS_TRANSLATING_PATH}/smiles_min.py
${SRC_CURATING_EDITING_CHEMO_SUBSCRIPTS_TRANSLATING_PATH}/smiles_min.py: ${INTERIM_TABLE_ORIGINAL_PATH}/originalStructureSmiles.tsv.zip
	cd src && python 2_curating/2_editing/chemo/subscripts/1_translating/smiles_min.py

is 


dvc run -n smiler \
          -d src/2_curating/2_editing/chemo/subscripts/1_translating/smiles_min.py -d data/interim/tables/0_original/originalStructureSmiles.tsv.zip \
          -o data/dvc_pipeline_outputs/translatedStructureSmiles_min.tsv.zip \
          python src/2_curating/2_editing/chemo/subscripts/1_translating/smiles_min.py


Let's try to traduce this one 

curating-editing-chemo-name: ${INTERIM_TABLE_TRANSLATED_PATH}/translatedStructureNominal.tsv.zip

${INTERIM_TABLE_TRANSLATED_PATH}/translatedStructureNominal.tsv.zip: ${SRC_CURATING_EDITING_CHEMO_SUBSCRIPTS_TRANSLATING_PATH}/names.R
${SRC_CURATING_EDITING_CHEMO_SUBSCRIPTS_TRANSLATING_PATH}/names.R: ${INTERIM_TABLE_ORIGINAL_PATH}/originalStructureNominal.tsv.zip
	cd src && Rscript 2_curating/2_editing/chemo/subscripts/1_translating/names.R

dvc run -n namer \
          -d src/2_curating/2_editing/chemo/subscripts/1_translating/names.R -d data/interim/tables/0_original/originalStructureNames.tsv.zip \
          -o data/dvc_pipeline_outputs/translatedStructureSmiles_min.tsv.zip \
          R src/2_curating/2_editing/chemo/subscripts/1_translating/smiles_min.py