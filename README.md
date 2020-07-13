# OpenNPDB: an Open Natural Products Database

*OpenNPDB:* an Open Natural Products Database. Actually, this database (DB) consists of XXX'XXX structure - organism pairs, gathered and standardized among XX open DB’s. It represents the most exhaustive open DB for natural products (NP’s) dereplication. It encompasses XXX'XXX distinct structures in XX’XXX resolved organisms. Additionally, both chemical and biological taxonomy are given for each pair. Other basic chemical descriptors and in silico tandem mass spectrometry (MS/MS) spectrum are provided for each structure. Within the frame of current computational approaches to guide NP’s research, all these elements should allow a more complete understanding of organisms and their metabolites.


## Flowchart

```mermaid
graph TD

subgraph legend
style legend fill:#FFFFFF,stroke:#424242,stroke-width:2px
A([file])
B[[script]]
C[(database)]
end

subgraph 1_gathering
010[(external/db/...)] -- n times --> 020[[db/.../standardizing.R]] -- n times --> 030([interim/db/...]) 
end

subgraph 2_curating
subgraph 1_integrating
030([interim/db/...]) --> 040[[integrating.R]]
040[[integrating.R]] --> 100([0_original/organism/organism/*.tsv])
040[[integrating.R]] --> 210([0_original/structure/inchi])
040[[integrating.R]] --> 220([0_original/structure/smiles])
040[[integrating.R]] --> 230([0_original/structure/nominal])
040[[integrating.R]] --> 310([0_original/reference/full])
040[[integrating.R]] --> 320([0_original/reference/doi])
040[[integrating.R]] --> 330([0_original/reference/pubmed])
040[[integrating.R]] --> 340([0_original/reference/title])
040[[integrating.R]] --> 350([0_original/reference/unsplit])

040[[integrating.R]] --> 400([0_original/table])
end

subgraph 2_editing
subgraph organism
style organism fill:#E5F5E0
100([0_original/organism/organism/*.tsv]) -->
    101[[1_cleaningOriginal.R]] --> 
        102([2_cleaned/organism/original/*.json])
    101[[1_cleaningOriginal.R]] --> 
        104([original organisms with taxonomy]) -->
                                    112[[3_cleaningTranslatedOrganism.R]]
    102([2_cleaned/organism/original/*.json]) -->
        107[[2_translatingOrganism.R]] -->
            108([1_translated/organism/*.tsv]) -->
                109[[3_cleaningTranslatedOrganism.R]] --> 
                    111([translated organisms with taxonomy]) -->
                        112[[4_cleaningTaxonomy.R]] -->
                            120([2_cleaned/organism/final])
                109[[3_cleaningTranslatedOrganism.R]] --> 110([2_cleaned/organism/translated/*.json]) --> 109[[3_cleaningTranslatedOrganism.R]]
end

subgraph structure
style structure fill:#FEE6CE
210([0_original/structure/inchi]) -->
            240[[2_integrating.R]]
220([0_original/structure/smiles]) -->
    221[[1_translating/smiles.py]] -->
        222([1_translated/structure/smiles]) -->
            240[[2_integrating.R]]
230([0_original/structure/nominal]) -->
    231[[1_translating/names.R]] -->
        232([1_translated/structure/names]) -->
            240[[2_integrating.R]] -->
                250([1_translated/structure/unique]) -->
                    260[[3_CleaningAndEnriching/chemosanitizer.py]] -->
                        270([1_translated/structure/cleaned]) -->
                            |external| 281[[classyfire]] -->
                                |external| 291([structures enriched taxonomy]) -->
                                    298[[integrating enrichment]]
                        270([1_translated/structure/cleaned]) -->
                            |external| 282[[chemGPS]] -->
                                |external| 292([structures enriched chemGPS]) -->
                                    298[[integrating enrichment]]
                        270([1_translated/structure/cleaned]) -->
                            |external| 283[[in silico]] -->
                                |external| 293([structures enriched in silico spectra]) -->
                                    298[[integrating enrichment]] -->
                                        299([clean and enriched structures])

classDef NotDone color:red
class 281,291,282,292,283,293 NotDone

end

subgraph reference
style reference fill:#EFEDF5
320([0_original/reference/doi]) -->
    321[[1_translating/doi.R]] -->
        322([1_translated/reference/doi]) -->
            360[[2_integrating/integrating.R]]
330([0_original/reference/pubmed]) -->
    331[[1_translating/pubmed.R]] -->
        332([1_translated/reference/pubmed]) -->
            360[[2_integrating/integrating.R]]
340([0_original/reference/title]) -->
    341[[1_translating/title.R]] -->
        342([1_translated/reference/title]) -->
            360[[2_integrating/integrating.R]]
350([0_original/reference/unsplit]) -->
    351[[1_translating/unsplit.R]] -->
        352([1_translated/reference/unsplit]) -->
            360[[2_integrating/integrating.R]] -->
                370([1_translated/reference/reference]) -->
                    380[[3_cleaning.R]] -->
                        390([2_cleaned/reference/cleaned.tsv.zip])
end

subgraph 3_integrating
120([2_cleaned/organism/cleaned.tsv.zip]) -->
998[[integrating.R]]

299([clean and enriched structures]) -->
998[[integrating.R]]

390([2_cleaned/reference/cleaned.tsv.zip]) -->
998[[integrating.R]]

400([0_original/table]) --> 
    998[[integrating.R]] --> 999[(2_cleaned/table)]
end
end    
end
```

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

You need Git, [DVC](https://www.dvc.org).

You need to have access to the metabomaps server to be able to pull the data, if you don't you will have to pull all the DBs data
from scratch.

### Access to metabomaps

Add this to your ~/.ssh/config

```
Host metabomaps
  IdentityFile /home/<user_local>/.ssh/id_rsa_metabomaps
  User <user_remote>
  HostName metabomaps.nprod.net
  Port 10311
```

### Pull the repository

```console
git pull https://gitlab.unige.ch/Adriano.Rutz/opennaturalproductsdb.git

# If you need the data
dvc pull  # This will take a while
```

### Having the data in a different place

If you want to have the data in a different place (for example for running a test), you can set the environment variable **DATA_PATH**.

### If you want to build only a simple DB

make -C src/1_gathering/db -B alkamid 


### Use docker to build

Install docker on your machine, make sure it is on your path

then

````console
make docker-build
make docker-bash
````

This will bring you in a container that will already have all the dependencies installed so you can run your commands in it.


# TO UPDATE

```
Give examples
```

### Packages

#### Conda environement

A "loose" environment_loose.yml file is created and should allow to recreate a working env formthe project without beeing too restrictive on the versions to install.
Install it by running in the home directory. Else use the environment_notloose.yml for strict mirroring of the installed packages.

````
conda env create -f environment_loose.yml
````

If your environement is not directly sources by your default bash run the following lines

````
$ source ~/anaconda3/etc/profile.d/conda.sh

$ conda activate your_env
````

Your R working directory should be src

````
$ cd src
````

If you are using Visual Studio be sure to set your R path in the settings option to reflect your created conda environment.

## Minimal working example

A minimal working example containing XXX entries coming from various DB's is proposed.
Use this example to check if all steps are running correctly on your machine.

## To build a graph from the make
(Requires remake and gprof2dot https://github.com/jrfonseca/gprof2dot)

```
remake --profile -B curating
gprof2dot -f callgrind callgrind.out.50802 | dot -Tpng -o output_full.png

```



### Description

The minimal working example file contains following columns:

# TO UPDATE

```
Give an example
```

### Final output

If everything went well, the output of the minimal working example should be:

# TO UPDATE

```
Give an example 
```

## Explanations

Add additional notes 

## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

Should we include versioning?

## Authors

* **Adriano Rutz** - *Initial work* - [Adriano.Rutz](https://gitlab.unige.ch/Adriano.Rutz)
* **Pierre-Marie Allard** - *Investigator* - [Pierre-Marie.Allard](https://gitlab.unige.ch/Pierre-Marie.Allard)
* **Jonathan Bisson** - *To define* - [bjo](https://gitlab.unige.ch/bjo)

See also the list of [contributors](https://gitlab.unige.ch/Adriano.Rutz/opennaturalproductsdb/-/project_members) who participated in this project.

## License

This project is licensed under the GNU GPLv3 license - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc.
