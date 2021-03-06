```mermaid

graph TD

subgraph legend
style legend fill:#FFFFFF
A([file in ../data/..])
B[[script]]
end

subgraph src
style src fill:#FFFFFF

subgraph 1_gathering
style 1_gathering fill:#FFFFFF

subgraph db
style db fill:#FFFFFF
010([external/dbSource/..]) -- x times -->
020[[db/../standardizing.R]] -- x times -->
030([interim/db/..])
end

subgraph translation
style translation fill:#A6CEE3
010([external/dbSource/..])
011([external/translation/..]) -- y times -->
040[[common.R]] -->
050([interim/dictionaries/common/..])
010([external/dbSource/..])
011([external/translation/..]) -- z times -->
060[[tcm.R]] -->
070([interim/dictionaries/tcm/..])
end
end

050([interim/dictionaries/common/..]) -->
105[[2_translating.kt]]
070([interim/dictionaries/tcm/..]) -->
105[[2_translating.kt]]

subgraph 2_curating
style 2_curating fill:#FFFFFF
030([interim/db/...]) -->
080[[1_integrating.R]]
080[[1_integrating.R]] -->
100([interim/tables/0_original/organism/*.tsv])
080[[1_integrating.R]] -->
210([interim/tables/0_original/structure/inchi.tsv.gz])
080[[1_integrating.R]] -->
220([interim/tables/0_original/structure/smiles.tsv.gz])
080[[1_integrating.R]] -->
230([interim/tables/0_original/structure/nominal.tsv.gz])
080[[1_integrating.R]] -->
320([interim/tables/0_original/reference/doi.tsv.gz])
080[[1_integrating.R]] -->
330([interim/tables/0_original/reference/pubmed.tsv.gz])
080[[1_integrating.R]] -->
340([interim/tables/0_original/reference/title.tsv.gz])
080[[1_integrating.R]] -->
350([interim/tables/0_original/reference/split.tsv.gz])
080[[1_integrating.R]] -->
361([interim/tables/0_original/reference/pubDetails.tsv.gz])
080[[1_integrating.R]] -->
370([interim/tables/0_original/reference/original.tsv.gz])
080[[1_integrating.R]] -->
400([interim/tables/0_original/table.tsv.gz])

subgraph 2_editing
style 2_editing fill:#FFFFFF

subgraph organism
style organism fill:#A6CEE3
100([interim/tables/0_original/organism/*.tsv]) -->
101[[1_cleaningOriginal.R]] -->
102([interim/tables/2_cleaned/organism/original/*.json]) -->
103[[1_cleaningOriginal.R]] -->
104([interim/tables/2_cleaned/organism/original.tsv.gz]) -->
111[[4_cleaningTaxonomy.R]]
104([interim/tables/2_cleaned/organism/original.tsv.gz]) -->
105[[2_translating_organism/../main.kt]] -->
106([interim/1_translated/organism/*.tsv]) -->
107[[3_cleaningTranslated.R]] -->
108([interim/tables/2_cleaned/organism/translated/*.json]) -->
109[[3_cleaningTranslated.R]] -->
110([interim/tables/2_cleaned/organism/translated.tsv.gz]) -->
111[[4_cleaningTaxonomy.R]] -->
120([interim/tables/2_cleaned/organism/verify.tsv.gz]) -->
112[[4_cleaningTaxonomy.R]] -->
121([interim/tables/2_cleaned/organism/verified.json]) -->
113[[4_cleaningTaxonomy.R]] -->
122([interim/tables/2_cleaned/organism/cleaned.tsv.gz]) -->
123[[5_addingOTL.R]] -->
124([interim/dictionaries/organism/otl.sqlite]) -->
123[[5_addingOTL.R]]
end

subgraph structure
style structure fill:#B2DF8A
210([interim/tables/0_original/structure/inchi.tsv.gz]) -->
240[[2_integrating.R]]
220([interim/tables/0_original/structure/smiles.tsv.gz]) -->
221[[smiles.py]] -->
222([interim/tables/1_translated/structure/smiles.tsv.gz]) -->
240[[2_integrating.R]]
230([interim/tables/0_original/structure/nominal.tsv.gz]) -->
231[[names.R]] -->
232([interim/tables/1_translated/structure/names.tsv.gz]) -->
240[[2_integrating.R]]-->
251([interim/tables/1_translated/structure/final.tsv.gz])
240[[2_integrating.R]]-->
250([interim/tables/1_translated/structure/unique.tsv.gz]) -->
260[[3_cleaningAndEnriching/sanitizing.py]] -->
270([interim/tables/2_cleaned/structure/cleaned.tsv.gz]) -->
280[[3_cleaningAndEnriching/stereocounting.py]] -->
281([interim/tables/2_cleaned/structure/counted.tsv.gz]) -->
290[[4_enriching/naming.R]] -->
291([interim/tables/2_cleaned/structure/named.tsv.gz])
281([interim/tables/2_cleaned/structure/counted.tsv.gz]) -->
292[[4_enriching/np-classifier.R]] -->
293([interim/dictionaries/structure/npclassifier/smiles_np_classified.tsv.gz])
end

subgraph reference
style reference fill:#FB9A99
320([interim/tables/0_original/reference/doi.tsv.gz]) -->
321[[1_translating/doi.R]] -->
322([interim/tables/1_translated/reference/doi.tsv.gz]) -->
360[[2_integrating.R]]
330([interim/tables/0_original/reference/pubmed.tsv.gz]) -->
331[[1_translating/pubmed.R]] -->
332([interim/tables/1_translated/reference/pubmed.tsv.gz]) -->
360[[2_integrating.R]]
340([interim/tables/0_original/reference/title.tsv.gz]) -->
341[[1_translating/title.R]] -->
342([interim/tables/1_translated/reference/title.tsv.gz]) -->
360[[2_integrating.R]]
350([interim/tables/0_original/reference/split.tsv.gz]) -->
351[[1_translating/split.R]] -->
352([interim/tables/1_translated/reference/split.tsv.gz]) -->
360[[2_integrating.R]]
361([interim/tables/0_original/reference/pubDetails.tsv.gz]) -->
363[[1_translating/pubDetails.R]] -->
362([interim/tables/1_translated/reference/pubDetails.tsv.gz]) -->
360[[2_integrating.R]]
370([interim/tables/0_original/reference/original.tsv.gz]) -->
371[[1_translating/original.R]] -->
372([interim/tables/1_translated/reference/original.tsv.gz]) -->
360[[2_integrating.R]] -->
379([interim/dictionaries/reference/dictionary.tsv.gz])
360[[2_integrating.R]] -->
380([interim/tables/1_translated/reference/integrated.tsv.gz]) -->
385[[3_cleaning.R]] -->
390([interim/tables/2_cleaned/reference/cleaned.tsv.gz])
end

122([interim/tables/2_cleaned/organism/cleaned.tsv.gz]) -->
360[[2_integrating.R]]
122([interim/tables/2_cleaned/organism/cleaned.tsv.gz]) -->
998[[3_integrating.R]]
251([interim/tables/1_translated/structure/final.tsv.gz])-->
998[[3_integrating.R]]
291([interim/tables/2_cleaned/structure/named.tsv.gz]) -->
998[[3_integrating.R]]
293([interim/dictionaries/structure/npclassifier/smiles_np_classified.tsv.gz]) -->
998[[3_integrating.R]]
390([interim/tables/2_cleaned/reference/cleaned.tsv.gz]) -->
998[[3_integrating.R]] 
400([interim/tables/0_original/table.tsv.gz]) -->
998[[3_integrating.R]] -->
999([interim/tables/3_curated/table.tsv.gz])
998[[3_integrating.R]] -->
1001([interim/dictionary/organism/dictionary.tsv.gz])
998[[3_integrating.R]] -->
1002([interim/dictionary/organism/metadata.tsv.gz])
998[[3_integrating.R]] -->
1003([interim/dictionary/structure/dictionary.tsv.gz])
998[[3_integrating.R]] -->
1004([interim/dictionary/structure/metadata.tsv.gz])
998[[3_integrating.R]] -->
1005([interim/dictionary/reference/dictionaryOrganism.tsv.gz])
998[[3_integrating.R]] -->
1006([interim/dictionary/reference/metadata.tsv.gz])
998[[3_integrating.R]] -->
1007([interim/tables/3_curated/table.tsv.gz])
end
end

subgraph 3_analyzing
style 3_analyzing fill:#FFFFFF
1007([interim/tables/3_curated/table.tsv.gz]) -->
1010[[1_sampling.R]] -->
1015([interim/tables/validation/..]) -->
1020[[2_validating.R]]
1001([interim/dictionary/organism/dictionary.tsv.gz]) -->
1020[[2_validating.R]]
1002([interim/dictionary/organism/metadata.tsv.gz]) -->
1020[[2_validating.R]]
1003([interim/dictionary/structure/dictionary.tsv.gz]) -->
1020[[2_validating.R]]
1004([interim/dictionary/structure/metadata.tsv.gz]) -->
1020[[2_validating.R]]
1005([interim/dictionary/reference/dictionaryOrganism.tsv.gz]) -->
1020[[2_validating.R]]
1006([interim/dictionary/reference/metadata.tsv.gz]) -->
1020[[2_validating.R]]
1007([interim/tables/3_curated/table.tsv.gz]) -->
1020[[2_validating.R]] -->
1030([interim/tables/4_analysed/platinum.tsv.gz])
end

subgraph 4_visualizing
style 4_visualizing fill:#FFFFFF
1030([interim/tables/4_analysed/platinum.tsv.gz]) -->
1040[[visualizing.R]] -->
1050([figures])
end
1030([interim/tables/4_analysed/platinum.tsv.gz]) -->
1060{"EXPORT TO WIKIDATA - npwdimporter #9829;"}
style 1060 fill:#f9f,stroke:#333,stroke-width:4px
end

```
