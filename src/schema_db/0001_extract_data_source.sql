SELECT *
FROM data_source
LEFT JOIN databases_source 
  ON data_source.databaseSourceId = databases_source.id
LEFT JOIN organisms_source 
    ON data_source.organismSourceId = organisms_source.id
LEFT JOIN references_source 
    ON data_source.referenceSourceId = references_source.id
LEFT JOIN structures_source 
    ON data_source.structureSourceId = structures_source.id
    ; 
