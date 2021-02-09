SELECT data_source.id,
       database_source.name   AS database,
       organism_source.value  AS \"organismValue\",
       organism_type.name     AS \"organismType\",
       reference_source.value AS \"referenceValue\",
       reference_type.name    AS \"referenceType\",
       structure_source.value AS \"structureValue\",
       structure_type.name    AS \"structureType\"
FROM data_source
         LEFT JOIN database_source
                   ON data_source.\"databaseSourceId\" = database_source.id
         LEFT JOIN organism_source
                   ON data_source.\"organismSourceId\" = organism_source.id
         LEFT JOIN organism_type
                   ON organism_source.\"organismTypeId\" = organism_type.id
         LEFT JOIN reference_source
                   ON data_source.\"referenceSourceId\" = reference_source.id
         LEFT JOIN reference_type
                   ON reference_source.\"referenceTypeId\" = reference_type.id
         LEFT JOIN structure_source
                   ON data_source.\"structureSourceId\" = structure_source.id
         LEFT JOIN structure_type
                   ON structure_source.\"structureTypeId\" = structure_type.id;