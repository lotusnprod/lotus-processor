create table curation_type
(
    id   serial primary key,
    name text
);

create table structure_cleaned
(
    id                        serial primary key,
    traditional_name          text,
    iupac_name                text,
    inchikey                  text,
    inchikey2d                text,
    inchi                     text,
    inchi2d                   text,
    smiles                    text,
    smiles2d                  text,
    stereocenters_total       integer,
    stereocenters_unspecified integer,
    molecular_formula         text,
    exact_mass                numeric,
    xlogp                     numeric
);

create table structure_information_type
(
    id   serial primary key,
    name text
);

create table structure_information
(
    id                            serial primary key,
    structure_information_type_id integer
        constraint structure_information_structure_information_type_id_fk
            references structure_information_type,
    value_string                  text,
    value_int                     integer,
    value_real                    numeric
);

create table structure_type
(
    id   serial primary key,
    name text
);

create table structure_source
(
    id                serial primary key,
    value             text,
    structure_type_id integer
        constraint structure_source_structure_type_id_fk
            references structure_type
);

create table structure_cleaned__structure_information
(
    id                       serial primary key,
    structure_cleaned_id     integer
        constraint structure_cleaned__structure_information_structure_cleaned_id_f
            references structure_cleaned,
    structure_information_id integer
        constraint structure_cleaned__structure_information_structure_information_
            references structure_information
);

create table organism_cleaned
(
    id   serial primary key,
    name text
);

create table organism_synonym
(
    id                  serial primary key,
    name                text,
    organism_cleaned_id integer
        constraint organism_synonym_organism_cleaned_id_fk
            references organism_cleaned
);

create table organism_database
(
    id   serial primary key,
    name text
);

create table organism_information
(
    id                   serial primary key,
    organism_cleaned_id  integer
        constraint organism_information_organism_cleaned_id_fk
            references organism_cleaned,
    organism_database_id integer
        constraint organism_information_organism_database_id_fk
            references organism_database,
    taxon_id             text,
    ranks                text,
    taxonomy             text,
    rank                 text
);

create table reference_cleaned
(
    id    serial primary key,
    doi   text,
    pmcid text,
    pmid  text,
    title text
);

create table reference_database
(
    id   serial primary key,
    name text
);

create table reference_information
(
    id                    serial primary key,
    reference_cleaned_id  integer
        constraint reference_information_reference_cleaned_id_fk
            references reference_cleaned,
    reference_database_id integer
        constraint reference_information_reference_database_id_fk
            references reference_database,
    reference_id          text,
    data                  text
);

create table reference_type
(
    id   serial primary key,
    name text
);

create table reference_source
(
    id                serial primary key,
    value             text,
    reference_type_id integer
        constraint reference_source_reference_type_id_fk
            references reference_type
);

create table data_cleaned
(
    id                   serial primary key,
    structure_cleaned_id integer
        constraint data_cleaned_structure_cleaned_id_fk
            references structure_cleaned,
    organism_cleaned_id  integer
        constraint data_cleaned_organism_cleaned_id_fk
            references organism_cleaned,
    reference_cleaned_id integer
        constraint data_cleaned_reference_cleaned_id_fk
            references reference_cleaned,
    curation_type_id     integer
        constraint data_cleaned_curation_type_id_fk
            references curation_type
);

create table database_type
(
    id   serial primary key,
    name text
);

create table database_source
(
    id               serial primary key,
    name             text,
    database_type_id integer
        constraint database_source_database_type_id_fk
            references database_type
);

create table organism_type
(
    id   serial primary key,
    name text
);

create table organism_source
(
    id               serial primary key,
    value            text,
    organism_type_id integer
        constraint organism_source_organism_type_id_fk
            references organism_type
);

create table data_source
(
    id                  serial primary key,
    database_source_id  integer
        constraint data_source_database_source_id_fk
            references database_source,
    organism_source_id  integer
        constraint data_source_organism_source_id_fk
            references organism_source,
    structure_source_id integer
        constraint data_source_structure_source_id_fk
            references structure_source,
    reference_source_id integer
        constraint data_source_reference_source_id_fk
            references reference_source
);

create table data_source__data_cleaned
(
    id              serial primary key,
    data_source_id  integer
        constraint data_source__data_cleaned_data_source_id_fk
            references data_source,
    data_cleaned_id integer
        constraint data_source__data_cleaned_data_cleaned_id_fk
            references data_cleaned
);

create table ott_taxonomy
(
    uid        integer primary key,
    parent_uid integer
        constraint ott_information_parent_id_fk
            references ott_taxonomy,
    name       text,
    rank       text,
    sourceinfo text,
    uniqname   text,
    flags      text
);
