SELECT data_source.id, 
database_source.name AS database, 
organism_source.value AS organism, 
reference_source.value AS reference, 
structure_source.value as structure
FROM data_source
LEFT JOIN database_source 
  ON data_source.databaseSourceId = database_source.id
LEFT JOIN organism_source 
    ON data_source.organismSourceId = organism_source.id
LEFT JOIN reference_source 
    ON data_source.referenceSourceId = reference_source.id
LEFT JOIN structure_source 
    ON data_source.structureSourceId = structure_source.id
;
