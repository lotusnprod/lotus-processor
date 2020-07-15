# OpenNPDB: an Open Natural Products Database

*OpenNPDB:* an Open Natural Products Database. Actually, this database (DB) consists of XXX'XXX structure - organism pairs, gathered and standardized among XX open DB’s. It represents the most exhaustive open DB for natural products (NP’s) dereplication. It encompasses XXX'XXX distinct structures in XX’XXX resolved organisms. Additionally, both chemical and biological taxonomy are given for each pair. Other basic chemical descriptors and in silico tandem mass spectrometry (MS/MS) spectrum are provided for each structure. Within the frame of current computational approaches to guide NP’s research, all these elements should allow a more complete understanding of organisms and their metabolites.


## Flowchart

```mermaid
graph TD

subgraph legend
style legend fill:#FFFFFF,stroke:#424242,stroke-width:2px
A([file])
B[[script]]
end

subgraph 1_gathering
subgraph db
style db fill:#ffffcc
010([external/db/...]) -- x times --> 020[[db/.../standardizing.R]] -- x times --> 030([interim/db/...]) 
end

subgraph translation
style translation fill:#E5F5E0
010([external/db/...]) -- y times --> 040[[translation/common.R]]
011([external/translation/...]) -- z times --> 040[[translation/common.R]] --> 050([common/names])
010([external/db/...]) -- y2 times --> 060[[translation/tcm.R]]
011([external/translation/...]) -- z2 times --> 060[[translation/tcm.R]] --> 070([tcm/names])
end

end

050([common/names]) --> 105[[2_translating.kt]]
070([tcm/names]) --> 105[[2_translating.kt]]

subgraph 2_curating
subgraph 1_integrating
030([interim/db/...]) --> 080[[integrating.R]]
080[[integrating.R]] --> 100([organism/*.tsv])
080[[integrating.R]] --> 210([inchi.tsv.gz])
080[[integrating.R]] --> 220([smiles.tsv.gz])
080[[integrating.R]] --> 230([nominal.tsv.gz])
080[[integrating.R]] --> 310([full.tsv.gz])
080[[integrating.R]] --> 320([doi.tsv.gz])
080[[integrating.R]] --> 330([pubmed.tsv.gz])
080[[integrating.R]] --> 340([title.tsv.gz])
080[[integrating.R]] --> 350([unsplit.tsv.gz])
080[[integrating.R]] --> 400([0_original/table])
end

subgraph 2_editing
subgraph organism
style organism fill:#E5F5E0
100([organism/*.tsv]) -->
    101[[1_cleaningOriginal.R]] --> 
        102([original/*.json]) --> 
            103[[1_cleaningOriginal.R]] --> 
                104([original.tsv.gz]) -->
                    111[[4_cleaningTaxonomy.R]]
        102([original/*.json]) -->
            105[[2_translating.kt]] -->
                106([organism/*.tsv]) -->
                    107[[3_cleaningTranslated.R]] --> 
                        108([translated/*.json]) --> 
                            109[[3_cleaningTranslated.R]] -->  
                                110([translated.tsv.gz]) -->
                                    111[[4_cleaningTaxonomy.R]] -->
                                        120([cleaned.tsv.gz])
end

subgraph structure
style structure fill:#FEE6CE
subgraph str_1_translating
style str_1_translating fill:#FEE6CE
220([smiles.tsv.gz]) -->
    221[[smiles.py]] -->
        222([smiles.tsv.gz])
230([nominal.tsv.gz]) -->
    231[[names.R]] -->
        232([names.tsv.gz])
end

subgraph str_2_integrating
style str_2_integrating fill:#FEE6CE
210([inchi.tsv.gz]) --> 240[[integrating.R]] 
222([smiles.tsv.gz]) --> 240[[integrating.R]]
232([names.tsv.gz]) --> 240[[integrating.R]]
    --> 250([unique.tsv.gz]) 
end

subgraph str_3_sanitizing
style str_3_sanitizing fill:#FEE6CE
250([unique.tsv.gz]) --> 260[[chemosanitizer.py]] -->
    270([cleaned.tsv.gz])
end

subgraph str_4_enriching
style str_4_enriching fill:#FEE6CE
270([cleaned.tsv.gz]) -->
    |external| 281[[classyfire]] -->
        |external| 291([enriched taxonomy])
270([cleaned.tsv.gz]) -->
    |external| 282[[chemGPS]] -->
        |external| 292([enriched chemGPS])
270([cleaned.tsv.gz]) -->
    |external| 283[[in silico]] -->
        |external| 293([enriched in silico spectra])
end 

subgraph str_5_integrating
style str_5_integrating fill:#FEE6CE
291([enriched taxonomy]) --> 298[[integrating enrichment]]
292([enriched chemGPS]) --> 298[[integrating enrichment]]
293([enriched in silico spectra]) --> 298[[integrating enrichment]]
298[[integrating enrichment]] --> 299([clean and enriched structures])
end
classDef NotDone color:red
class 281,291,282,292,283,293,298 NotDone

end

subgraph reference
style reference fill:#EFEDF5
subgraph ref_1_translating
style ref_1_translating fill:#EFEDF5
320([doi.tsv.gz]) -->
    321[[doi.R]] -->
        322([doi.tsv.gz]) 
330([pubmed.tsv.gz]) -->
    331[[pubmed.R]] -->
        332([pubmed.tsv.gz])
340([title.tsv.gz]) -->
    341[[title.R]] -->
        342([title.tsv.gz])
350([unsplit.tsv.gz]) -->
    351[[unsplit.R]] -->
        352([unsplit.tsv.gz])
end
subgraph ref_2_integrating
style ref_2_integrating fill:#EFEDF5
310([full.tsv.gz]) --> 360[[2_integrating/integrating.R]]
322([doi.tsv.gz]) --> 360[[2_integrating/integrating.R]]
332([pubmed.tsv.gz]) --> 360[[2_integrating/integrating.R]]
342([title.tsv.gz]) --> 360[[2_integrating/integrating.R]]
352([unsplit.tsv.gz]) --> 360[[2_integrating/integrating.R]] -->
        370([reference.tsv.gz]) 
end
subgraph .
style . fill:#EFEDF5
370([reference.tsv.gz]) -->
    380[[3_cleaning.R]] -->
        390([cleaned.tsv.gz])
end
end

subgraph 3_integrating
120([cleaned.tsv.gz]) -->
998[[integrating.R]]

299([clean and enriched structures]) -->
998[[integrating.R]]

390([cleaned.tsv.gz]) -->
998[[integrating.R]]

400([0_original/table]) --> 
    998[[integrating.R]] --> 999([2_cleaned/table.tsv.gz])
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
* **Jonathan Bisson** - *Hacker in Chief* - [bjo](https://gitlab.unige.ch/bjo)

See also the list of [contributors](https://gitlab.unige.ch/Adriano.Rutz/opennaturalproductsdb/-/project_members) who participated in this project.

## License

This project is licensed under the GNU GPLv3 license - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc.
