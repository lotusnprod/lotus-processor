create table curation_type
(
    id                              serial primary key,
    name                            text
);

create table structure_cleaned
(
    id                              serial primary key,
    "traditionalName"               text,
    "iupacName"                     text,
    inchikey                        text,
    "inchikey2D"                    text,
    inchi                           text,
    "inchi2D"                       text,
    smiles                          text,
    "smiles2D"                      text,
    "stereocentersTotal"            integer,
    "stereocentersUnspecified"      integer,
    "molecularFormula"              text,
    "exactMass"                     numeric,
    xlogp                           numeric
);

create table structure_information_type
(
    id                              serial primary key,
    name                            text
);

create table structure_information
(
    id                              serial primary key,
    "structureInformationTypeId"    integer
                                    constraint structure_information_structure_information_type_id_fk
                                                references structure_information_type,
    "valueString"                   text,
    "valueInt"                      integer,
    "valueReal"                     numeric
);

create table structure_type
(
    id                              serial primary key,
    name                            text
);

create table structure_source
(
    id                              serial primary key,
    value                           text,
    "structureTypeId"               integer
                                    constraint structure_source_structure_type_id_fk
                                                references structure_type
);

create table structure_cleaned__structure_information
(
    id                              serial primary key,
    "structureCleanedId"            integer
                                    constraint structure_cleaned__structure_information_structure_cleaned_id_f
                                                references structure_cleaned,
    "structureInformationId"        integer
                                    constraint structure_cleaned__structure_information_structure_information_
                                                references structure_information
);

create table organism_cleaned
(
    id                              serial primary key,
    name                            text
);

create table organism_synonym
(
    id                              serial primary key,
    name                            text,
    "organismCleanedId"             integer
                                    constraint organism_synonym_organism_cleaned_id_fk
                                                references organism_cleaned
);

create table organism_database
(
    id                              serial primary key,
    name                            text
);

create table organism_information
(
    id                              serial primary key,
    "organismCleanedId"             integer
                                    constraint organism_information_organism_cleaned_id_fk
                                                references organism_cleaned,
    "organismDatabaseId"            integer
                                    constraint organism_information_organism_database_id_fk
                                                references organism_database,
    "taxonId"                       text,
    ranks                           text,
    taxonomy                        text,
    rank                            text
);

create table reference_cleaned
(
    id                              serial primary key,
    doi                             text,
    pmcid                           text,
    pmid                            text,
    title                           text
);

create table reference_database
(
    id                              serial primary key,
    name                            text
);

create table reference_information
(
    id                              serial primary key,
    "referenceCleanedId"            integer
                                    constraint reference_information_reference_cleaned_id_fk
                                                references reference_cleaned,
    "referenceDatabaseId"           integer
                                    constraint reference_information_reference_database_id_fk
                                                references reference_database,
    "referenceId"                   text,
    data                            text
);

create table reference_type
(
    id                              serial primary key,
    name                            text
);

create table reference_source
(
    id                              serial primary key,
    value                           text,
    "referenceTypeId"               integer
                                    constraint reference_source_reference_type_id_fk
                                                references reference_type
);

create table data_cleaned
(
    id                              serial primary key,
    "structureCleanedId"            integer
                                    constraint data_cleaned_structure_cleaned_id_fk
                                                references structure_cleaned,
    "organismCleanedId"             integer
                                    constraint data_cleaned_organism_cleaned_id_fk
                                                references organism_cleaned,
    "referenceCleanedId"            integer
                                    constraint data_cleaned_reference_cleaned_id_fk
                                                references reference_cleaned,
    "curationTypeId"                integer
                                    constraint data_cleaned_curation_type_id_fk
                                                references curation_type
);

create table database_type
(
    id                              serial primary key,
    name                            text
);

create table database_source
(
    id                              serial primary key,
    name                            text,
    "databaseTypeId"                integer
                                    constraint database_source_database_type_id_fk
                                                references database_type
);

create table organism_type
(
    id                              serial primary key,
    name                            text
);

create table organism_source
(
    id                              serial primary key,
    value                           text,
    "organismTypeId"                integer
                                    constraint organism_source_organism_type_id_fk
                                                references organism_type
);

create table data_source
(
    id                              serial primary key,
    "databaseSourceId"              integer
                                    constraint data_source_database_source_id_fk
                                                references database_source,
    "organismSourceId"              integer
                                    constraint data_source_organism_source_id_fk
                                                references organism_source,
    "structureSourceId"             integer
                                    constraint data_source_structure_source_id_fk
                                                references structure_source,
    "referenceSourceId"             integer
                                    constraint data_source_reference_source_id_fk
                                                references reference_source
);

create table data_source__data_cleaned
(
    id                              serial primary key,
    "dataSourceId"                  integer
                                    constraint data_source__data_cleaned_data_source_id_fk
                                                references data_source,
    "dataCleanedId"                 integer
                                    constraint data_source__data_cleaned_data_cleaned_id_fk
                                                references data_cleaned
);

create table ott_taxonomy
(
    uid                             integer primary key,
    parent_uid                      integer
                                    constraint ott_information_parent_id_fk
                                                references ott_taxonomy,
    name                            text,
    rank                            text,
    sourceinfo                      text,
    uniqname                        text,
    flags                           text
);
