create table curation_type
(
    id   INTEGER not null
        constraint curation_type_pk
            primary key autoincrement,
    name TEXT
);

create
    unique index curation_type_id_uindex
    on curation_type (id);

create
    unique index curation_type_name_uindex
    on curation_type (name);

create table data_cleaned
(
    id                 INTEGER not null
        constraint data_cleaned_pk
            primary key autoincrement,
    structureCleanedId INTEGER
        references structure_cleaned,
    organismCleanedId  INTEGER
        references organism_cleaned,
    referenceCleanedId INTEGER
        references reference_cleaned,
    curationTypeId     INTEGER
        references curation_type
);

create
    unique index data_cleaned_id_uindex
    on data_cleaned (id);

create table data_source
(
    id                INTEGER not null
        constraint data_source_pk
            primary key autoincrement,
    databaseSourceId  INTEGER
        references database_source,
    organismSourceId  INTEGER
        references organism_source,
    structureSourceId INTEGER
        references structure_source,
    referenceSourceId INTEGER
        references reference_source
);

create
    unique index data_source_id_uindex
    on data_source (id);

create table data_source__data_cleaned
(
    id            INTEGER not null
        constraint data_source__data_cleaned_pk
            primary key autoincrement,
    dataSourceId  INTEGER
        references data_source,
    dataCleanedId INTEGER
        references data_cleaned
);

create
    unique index data_source__data_cleaned_id_uindex
    on data_source__data_cleaned (id);

create table database_type
(
    id   INTEGER not null
        constraint database_type_pk
            primary key autoincrement,
    name TEXT
);

create table database_source
(
    id             INTEGER not null
        constraint database_source_pk
            primary key autoincrement,
    name           TEXT,
    databaseTypeId INTEGER
        references database_type
);

create
    unique index database_source_id_uindex
    on database_source (id);

create
    unique index database_source_name_uindex
    on database_source (name);

create table organism_cleaned
(
    id   INTEGER not null
        constraint organism_cleaned_pk
            primary key autoincrement,
    name TEXT
);

create
    unique index organism_cleaned_id_uindex
    on organism_cleaned (id);

create table organism_synonym
(
    id                INTEGER not null
        constraint organism_synonym_pk
            primary key autoincrement,
    name              TEXT,
    organismCleanedId int
        references organism_cleaned
);

create table organism_database
(
    id   INTEGER not null
        constraint organism_database_pk
            primary key autoincrement,
    name TEXT
);

create
    unique index organism_database_id_uindex
    on organism_database (id);

create table organism_information
(
    id                 INTEGER not null
        constraint organism_information_pk
            primary key autoincrement,
    organismCleanedId  INTEGER
        references organism_cleaned,
    organismDatabaseId INTEGER
        references organism_database,
    taxonId            TEXT,
    ranks              TEXT,
    taxonomy           TEXT,
    rank               TEXT
);

create
    unique index organism_information_id_uindex
    on organism_information (id);

create
    unique index organism_synonym_id_uindex
    on organism_synonym (id);

create table organism_type
(
    id   INTEGER not null
        constraint organism_type_pk
            primary key autoincrement,
    name TEXT
);

create table organism_source
(
    id             INTEGER not null
        constraint organism_source_pk
            primary key autoincrement,
    value          TEXT,
    organismTypeId INTEGER
        references organism_type
);

create table reference_cleaned
(
    id    INTEGER not null
        constraint reference_cleaned_pk
            primary key autoincrement,
    doi   TEXT,
    pmcid TEXT,
    pmid  TEXT,
    title TEXT
);

create
    unique index reference_cleaned_id_uindex
    on reference_cleaned (id);


create table reference_database
(
    id   INTEGER not null
        constraint reference_database_pk
            primary key autoincrement,
    name TEXT
);

create
    unique index reference_database_id_uindex
    on reference_database (id);

create table reference_information
(
    id                  INTEGER not null
        constraint reference_information_pk
            primary key autoincrement,
    referenceCleanedId  INTEGER
        references reference_cleaned,
    referenceDatabaseId INTEGER
        references reference_database,
    referenceId         TEXT,
    data                TEXT
);

create
    unique index reference_information_id_uindex
    on reference_information (id);

create table reference_type
(
    id   INTEGER not null
        constraint reference_type_pk
            primary key autoincrement,
    name TEXT
);

create table reference_source
(
    id              INTEGER not null
        constraint reference_source_pk
            primary key autoincrement,
    value           TEXT,
    referenceTypeId INTEGER
        references reference_type
);

create table structure_cleaned
(
    id                       INTEGER not null
        constraint structure_cleaned_pk
            primary key autoincrement,
    traditionalName          TEXT,
    iupacName                TEXT,
    inchikey                 TEXT,
    inchikey2D               TEXT,
    inchi                    TEXT,
    inchi2D                  TEXT,
    smiles                   TEXT,
    smiles2D                 TEXT,
    stereocentersTotal       INTEGER,
    stereocentersUnspecified INTEGER,
    molecularFormula         TEXT,
    exactMass                REAL,
    xlogp                    REAL
);

create
    unique index structure_cleaned_id_uindex
    on structure_cleaned (id);

create table structure_cleaned__structure_information
(
    id                     INTEGER not null
        constraint structure_cleaned__structure_information_pk
            primary key autoincrement,
    structureCleanedId     INTEGER
        references structure_cleaned,
    structureInformationId INTEGER
        references structure_information
);

create
    unique index structure_cleaned__structure_information_id_uindex
    on structure_cleaned__structure_information (id);

create table structure_information
(
    id                         INTEGER not null
        constraint structure_information_pk
            primary key autoincrement,
    structureInformationTypeId INTEGER
        references structure_information_type,
    valueString                TEXT,
    valueInt                   INTEGER,
    valueReal                  REAL
);

create
    unique index structure_information_id_uindex
    on structure_information (id);

create table structure_information_type
(
    id   INTEGER not null
        constraint structure_information_type_pk
            primary key autoincrement,
    name TEXT
);

create table structure_source
(
    id              INTEGER not null
        constraint structure_source_pk
            primary key autoincrement,
    value           TEXT,
    structureTypeId INTEGER
        references structure_type
);

create table structure_type
(
    id   INTEGER not null
        constraint structure_type_pk
            primary key autoincrement,
    name TEXT
);