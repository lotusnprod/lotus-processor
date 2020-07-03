
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


dvc run -n smiler \
          -d src/2_curating/2_editing/chemo/subscripts/1_translating/smiles_dvc.py -d data/interim/tables_min/0_original/smiles.tsv.zip \
          -o data/dvc_pipeline_outputs/translated_smiles.tsv.zip \
          python src/2_curating/2_editing/chemo/subscripts/1_translating/smiles_dvc.py data/interim/tables_min/0_original/smiles.tsv.zip data/dvc_pipeline_outputs/translated_smiles.tsv.zip structureOriginalSmiles


##Let's try to traduce this one 

curating-editing-chemo-name: ${INTERIM_TABLE_TRANSLATED_PATH}/translatedStructureNominal.tsv.zip

${INTERIM_TABLE_TRANSLATED_PATH}/translatedStructureNominal.tsv.zip: ${SRC_CURATING_EDITING_CHEMO_SUBSCRIPTS_TRANSLATING_PATH}/names.R
${SRC_CURATING_EDITING_CHEMO_SUBSCRIPTS_TRANSLATING_PATH}/names.R: ${INTERIM_TABLE_ORIGINAL_PATH}/originalStructureNominal.tsv.zip
	cd src && Rscript 2_curating/2_editing/chemo/subscripts/1_translating/names.R

dvc run -n namer \
          -d src/2_curating/2_editing/chemo/subscripts/1_translating/names_dvc.R -d data/interim/tables/0_original/originalStructureNominal.tsv.zip \
          -o data/dvc_pipeline_outputs/translatedStructureNominal.tsv.zip \
          Rscript src/2_curating/2_editing/chemo/subscripts/1_translating/names_dvc.R

Mooving to src to launch this one (too long and path mess)
dvc run -n namer \
          -d 2_curating/2_editing/chemo/subscripts/1_translating/names_dvc.R -d ../data/interim/tables/0_original/originalStructureNominal.tsv.zip \
          -o ../data/dvc_pipeline_outputs/translatedStructureNominal.tsv.zip \
          Rscript 2_curating/2_editing/chemo/subscripts/1_translating/names_dvc.R

## Working on the sanitizing script

dvc run -n sanitizer \
          -d src/2_curating/2_editing/chemo/subscripts/2_curatingAndEnriching/04_sanitizing.py -d data/interim/tables/1_translated/translatedTable.tsv.zip \
          -o data/dvc_pipeline_outputs/sanitized_structures.tsv.zip \
          python src/2_curating/2_editing/chemo/subscripts/2_curatingAndEnriching/04_sanitizing.py data/interim/tables/1_translated/translatedTable.tsv.zip data/dvc_pipeline_outputs/sanitized_structures.tsv.zip structureTranslated


Now adding parameters (these are kept in the params.yaml file) using the -p option


dvc run -n sanitizer \
          -d src/2_curating/2_editing/chemo/subscripts/2_curatingAndEnriching/04_sanitizing.py -d data/interim/tables/1_translated/translatedTable.tsv.zip \
          -p head_lenght \
          -o data/dvc_pipeline_outputs/sanitized_structures.tsv.zip \
          python src/2_curating/2_editing/chemo/subscripts/2_curatingAndEnriching/04_sanitizing.py data/interim/tables/1_translated/translatedTable.tsv.zip data/dvc_pipeline_outputs/sanitized_structures.tsv.zip structureTranslated


dvc run -n sanitizer \
          -d src/2_curating/2_editing/chemo/subscripts/3_cleaningAndEnriching/sanitizing.py -d data/dvc_pipeline_outputs/translated_smiles.tsv.zip \
          -p head_lenght \
          -o data/dvc_pipeline_outputs/sanitized_structures.tsv.zip \
          python src/2_curating/2_editing/chemo/subscripts/3_cleaningAndEnriching/sanitizing.py data/dvc_pipeline_outputs/translated_smiles.tsv.zip data/dvc_pipeline_outputs/sanitized_structures.tsv.zip structureTranslatedSmiles