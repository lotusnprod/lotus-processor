ALTER TABLE data_cleaned RENAME COLUMN structureCleanedId TO structure_cleaned_id;
ALTER TABLE data_cleaned RENAME COLUMN organismCleanedId TO organism_cleaned_id;
ALTER TABLE data_cleaned RENAME COLUMN referenceCleanedId TO reference_cleaned_id;
ALTER TABLE data_cleaned RENAME COLUMN curationTypeId TO curation_type_id;

ALTER TABLE data_source RENAME COLUMN databaseSourceId TO database_source_id;
ALTER TABLE data_source RENAME COLUMN organismSourceId TO organism_source_id;
ALTER TABLE data_source RENAME COLUMN structureSourceId TO structure_source_id;
ALTER TABLE data_source RENAME COLUMN referenceSourceId TO reference_source_id;

ALTER TABLE data_source__data_cleaned RENAME COLUMN dataSourceId TO data_source_id;
ALTER TABLE data_source__data_cleaned RENAME COLUMN dataCleanedId TO data_cleaned_id;

ALTER TABLE database_source RENAME COLUMN databaseTypeId TO database_type_id;

ALTER TABLE organism_information RENAME COLUMN organismCleanedId TO organism_cleaned_id;
ALTER TABLE organism_information RENAME COLUMN organismDatabaseId TO organism_database_id;
ALTER TABLE organism_information RENAME COLUMN taxonId TO taxon_id;

ALTER TABLE organism_source RENAME COLUMN organismTypeId TO organism_type_id;

ALTER TABLE organism_synonym RENAME COLUMN organismCleanedId TO organism_cleaned_id;

ALTER TABLE reference_information RENAME COLUMN referenceCleanedId TO reference_cleaned_id;
ALTER TABLE reference_information RENAME COLUMN referenceDatabaseId TO reference_database_id;
ALTER TABLE reference_information RENAME COLUMN referenceId TO reference_id;
2
ALTER TABLE reference_source RENAME COLUMN referenceTypeId TO reference_type_id;

ALTER TABLE structure_cleaned RENAME COLUMN traditionalName TO traditional_name;
ALTER TABLE structure_cleaned RENAME COLUMN iupacName TO iupac_name;
ALTER TABLE structure_cleaned RENAME COLUMN inchikey2D TO inchikey2d;
ALTER TABLE structure_cleaned RENAME COLUMN inchi2D TO inchi2d;
ALTER TABLE structure_cleaned RENAME COLUMN smiles2D TO smiles2d;
ALTER TABLE structure_cleaned RENAME COLUMN stereocentersTotal TO stereocenters_total;
ALTER TABLE structure_cleaned RENAME COLUMN stereocentersUnspecified TO stereocenters_unspecified;
ALTER TABLE structure_cleaned RENAME COLUMN molecularFormula TO molecular_formula;
ALTER TABLE structure_cleaned RENAME COLUMN exactMass TO exact_mass;

ALTER TABLE structure_cleaned__structure_information RENAME COLUMN structureCleanedId TO structure_cleaned_id;
ALTER TABLE structure_cleaned__structure_information RENAME COLUMN structureInformationId TO structure_information_id;

ALTER TABLE structure_information RENAME COLUMN structureInformationTypeId TO structure_information_type_id;
ALTER TABLE structure_information RENAME COLUMN valueString TO value_string;
ALTER TABLE structure_information RENAME COLUMN valueInt TO value_int;
ALTER TABLE structure_information RENAME COLUMN valueReal TO value_real;

ALTER TABLE structure_source RENAME COLUMN structureTypeId TO structure_type_id;

