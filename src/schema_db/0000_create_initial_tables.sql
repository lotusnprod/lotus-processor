create table chemical_databases
(
	id INTEGER not null
		constraint chemical_databases_pk
			primary key autoincrement,
	name TEXT
);

create unique index chemical_databases_id_uindex
	on chemical_databases (id);

create table chemical_information
(
	id INTEGER not null
		constraint chemical_information_pk
			primary key autoincrement,
	cleanedStructureId INTEGER
		references structure_cleaned,
	chemicalDatabaseId INTEGER
		references chemical_databases,
	chemicalId TEXT,
	ranks TEXT,
	taxonomy TEXT,
	rank TEXT
);

create unique index chemical_information_id_uindex
	on chemical_information (id);
	
create table curation_states
(
	id INTEGER not null
		constraint curation_states_pk
			primary key autoincrement,
	name TEXT
);

create unique index curation_states_id_uindex
	on curation_states (id);

create unique index curation_states_name_uindex
	on curation_states (name);

create table data_processed
(
	id INTEGER not null
		constraint processed_data_pk
			primary key autoincrement,
	structureCleanedId INTEGER
		references structure_cleaned,
	organismCleanedId INTEGER
		references organism_cleaned,
	referenceCleanedId INTEGER
		references reference_cleaned
);

create unique index data_processed_id_uindex
	on data_processed (id);

create table data_source
(
	id INTEGER not null
		constraint source_data_pk
			primary key autoincrement,
	databaseSourceId INTEGER
		references database_source,
	organismSourceId INTEGER
		references organism_source,
	structureSourceId INTEGER
		references structure_source,
	referenceSourceId INTEGER
		references reference_source
);

create unique index data_source_id_uindex
	on data_source (id);

create table data_processed__data_source
(
	id INTEGER not null
		constraint processed_data__source_data_pk
			primary key autoincrement,
	dataSourceId INTEGER
		references data_source,
	dataProcessedId INTEGER
		references data_processed,
	curationStateId INTEGER
		references curation_states
);

create unique index data_processed__data_source_id_uindex
	on data_processed__data_source (id);
	
create table database_types
(
	id INTEGER not null
		constraint database_types_pk
			primary key autoincrement,
	name TEXT
);

create table database_source
(
	id INTEGER not null
		constraint source_database_pk
			primary key autoincrement,
	name TEXT,
	typeId INTEGER
		references database_types
);
	
create unique index database_source_id_uindex
	on database_source (id);

create unique index database_source_name_uindex
	on database_source (name);

create table organism_cleaned
(
	id INTEGER not null
		constraint organism_cleaned_pk
			primary key autoincrement,
	name TEXT
);

create unique index organism_cleaned_id_uindex
	on organism_cleaned (id);

create table organism_synonyms
(
	id INTEGER not null
		constraint organism_synonyms_pk
			primary key autoincrement,
	name TEXT,
	organismCleanedId int
		references organism_cleaned
);

create unique index organism_synonyms_id_uindex
	on organism_synonyms (id);

create table organism_types
(
	id INTEGER not null
		constraint organism_types_pk
			primary key autoincrement,
	name TEXT
);

create table organism_source
(
	id INTEGER not null
		constraint organism_source_pk
			primary key autoincrement,
	value TEXT,
	typeId INTEGER
		references organism_types
);

create table reference_cleaned
(
	id INTEGER not null
		constraint reference_cleaned_pk
			primary key autoincrement,
	doi TEXT,
	pmcid TEXT,
	pmid TEXT,
	title TEXT
);

create unique index reference_cleaned_id_uindex
	on reference_cleaned (id);

create table reference_types
(
	id INTEGER not null
		constraint reference_types_pk
			primary key autoincrement,
	name TEXT
);

create table reference_source
(
	id INTEGER not null
		constraint reference_source_pk
			primary key autoincrement,
	value TEXT,
	typeId INTEGER
		references reference_types
);

create table structure_cleaned
(
	id INTEGER not null
		constraint cleaned_compounds_pk
			primary key autoincrement,
	traditionalName TEXT,
	iupacName TEXT,
	inchikey TEXT,
	inchikey2D TEXT,
	inchi TEXT,
	inchi2D TEXT,
	smiles TEXT,
	smiles2D TEXT,
	stereocentersTotal INTEGER,
	stereocentersUnspecified INTEGER,
	molecularFormula TEXT,
	exactMass REAL,
	xlogp REAL
);


create unique index structure_cleaned_id_uindex
	on structure_cleaned (id);

create table structure_types
(
	id INTEGER not null
		constraint structure_types_pk
			primary key autoincrement,
	name TEXT
);

create table structure_source
(
	id INTEGER not null
		constraint structure_source_pk
			primary key autoincrement,
	value TEXT,
	typeId INTEGER
		references structure_types
);

create table reference_databases
(
	id INTEGER not null
		constraint reference_databases_pk
			primary key autoincrement,
	name TEXT
);

create unique index reference_databases_id_uindex
	on reference_databases (id);

create table reference_information
(
	id INTEGER not null
		constraint reference_information_pk
			primary key autoincrement,
	cleanedReferenceId INTEGER
		references reference_cleaned,
	referenceDatabaseId INTEGER
		references reference_databases,
	referenceId TEXT,
	data TEXT
);

create unique index reference_information_id_uindex
	on reference_information (id);
	
create table taxonomic_databases
(
	id INTEGER not null
		constraint taxonomic_databases_pk
			primary key autoincrement,
	name TEXT
);

create unique index taxonomic_databases_id_uindex
	on taxonomic_databases (id);

create table taxonomic_information
(
	id INTEGER not null
		constraint taxonomic_information_pk
			primary key autoincrement,
	cleanedOrganismId INTEGER
		references organism_cleaned,
	taxonomicDatabaseId INTEGER
		references taxonomic_databases,
	taxonomicId TEXT,
	ranks TEXT,
	taxonomy TEXT,
	rank TEXT
);

create unique index taxonomic_information_id_uindex
	on taxonomic_information (id);