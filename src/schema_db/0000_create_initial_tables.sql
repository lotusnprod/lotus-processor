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
		references structures_cleaned,
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
		references structures_cleaned,
	organismCleanedId INTEGER
		references organisms_cleaned,
	referenceCleanedId INTEGER
		references references_cleaned,
	curationStateId INTEGER
		references curation_states
);

create unique index data_processed_id_uindex
	on data_processed (id);

create table data_source
(
	id INTEGER not null
		constraint source_data_pk
			primary key autoincrement,
	databaseSourceId INTEGER
		references databases_source,
	organismSourceId INTEGER
		references organisms_source,
	structureSourceId INTEGER
		references structures_source,
	referenceSourceId INTEGER,
	constraint data_source_references_source_id_id_fk
		foreign key (id, referenceSourceId) references references_source (id, id)
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
		references data_processed
);

create unique index data_processed__data_source_id_uindex
	on data_processed__data_source (id);
	
create table databases_types
(
	id INTEGER not null
		constraint databases_types_pk
			primary key autoincrement,
	name TEXT
);

create table databases_source
(
	id INTEGER not null
		constraint source_database_pk
			primary key autoincrement,
	name TEXT,
	typeId INTEGER
		references databases_types
);
	
create unique index databases_source_id_uindex
	on databases_source (id);

create unique index databases_source_name_uindex
	on databases_source (name);

create table organisms_cleaned
(
	id INTEGER not null
		constraint cleaned_organisms_pk
			primary key autoincrement,
	name TEXT
);

create unique index organisms_cleaned_id_uindex
	on organisms_cleaned (id);

create table organisms_synonyms
(
	id INTEGER not null
		constraint organismSynonym_pk
			primary key autoincrement,
	name TEXT,
	organismCleanedId int
		references organisms_cleaned
);

create unique index organisms_synonyms_id_uindex
	on organisms_synonyms (id);

create table organisms_types
(
	id INTEGER not null
		constraint organisms_types_pk
			primary key autoincrement,
	name TEXT
);

create table organisms_source
(
	id INTEGER not null
		constraint organisms_source_pk
			primary key autoincrement,
	value TEXT,
	typeId INTEGER
		references organisms_types
);

create table references_cleaned
(
	id INTEGER not null
		constraint cleaned_references_pk
			primary key autoincrement,
	doi TEXT,
	pmcid TEXT,
	pmid TEXT,
	title TEXT
);

create unique index references_cleaned_id_uindex
	on references_cleaned (id);

create table references_types
(
	id INTEGER not null
		constraint references_types_pk
			primary key autoincrement,
	name TEXT
);

create table references_source
(
	id INTEGER not null
		constraint references_source_pk
			primary key autoincrement,
	value TEXT,
	typeId INTEGER
		references references_types
);

create table structures_cleaned
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


create unique index structures_cleaned_id_uindex
	on structures_cleaned (id);

create table structures_types
(
	id INTEGER not null
		constraint structure_types_pk
			primary key autoincrement,
	name TEXT
);

create table structures_source
(
	id INTEGER not null
		constraint structures_source_pk
			primary key autoincrement,
	value TEXT,
	typeId INTEGER
		references structures_types
);

create table references_databases
(
	id INTEGER not null
		constraint references_databases_pk
			primary key autoincrement,
	name TEXT
);

create unique index references_databases_id_uindex
	on references_databases (id);

create table references_information
(
	id INTEGER not null
		constraint references_information_pk
			primary key autoincrement,
	cleanedReferenceId INTEGER
		references references_cleaned,
	referenceDatabaseId INTEGER
		references references_databases,
	referenceId TEXT,
	data TEXT
);

create unique index references_information_id_uindex
	on references_information (id);
	
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
		references organisms_cleaned,
	taxonomicDatabaseId INTEGER
		references taxonomic_databases,
	taxonomicId TEXT,
	ranks TEXT,
	taxonomy TEXT,
	rank TEXT
);

create unique index taxonomic_information_id_uindex
	on taxonomic_information (id);