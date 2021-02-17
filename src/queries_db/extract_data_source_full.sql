SELECT data_source.id,
       database_source.id     AS database_source_id,
       database_source.name   AS database,
       organism_source.id     AS organism_source_id,
       organism_source.value  AS organism_value,
       organism_type.name     AS organism_type,
       reference_source.id    AS reference_source_id,
       reference_source.value AS reference_value,
       reference_type.name    AS reference_type,
       structure_source.id    AS structure_source_id,
       structure_source.value AS structure_value,
       structure_type.name    AS structure_type
FROM data_source
         LEFT JOIN database_source
                   ON data_source.database_source_id = database_source.id
         LEFT JOIN organism_source
                   ON data_source.organism_source_id = organism_source.id
         LEFT JOIN organism_type
                   ON organism_source.organism_type_id = organism_type.id
         LEFT JOIN reference_source
                   ON data_source.reference_source_id = reference_source.id
         LEFT JOIN reference_type
                   ON reference_source.reference_type_id = reference_type.id
         LEFT JOIN structure_source
                   ON data_source.structure_source_id = structure_source.id
         LEFT JOIN structure_type
                   ON structure_source.structure_type_id = structure_type.id;