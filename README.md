# OpenNPDB: an Open Natural Products Database

*OpenNPDB:* an Open Natural Products Database. Actually, this database (DB) consists of XXX'XXX structure - organism pairs, gathered and standardized among XX open DB’s. It represents the most exhaustive open DB for natural products (NP’s) dereplication. It encompasses XXX'XXX distinct structures in XX’XXX resolved organisms. Additionally, both chemical and biological taxonomy are given for each pair. Other basic chemical descriptors and in silico tandem mass spectrometry (MS/MS) spectrum are provided for each structure. Within the frame of current computational approaches to guide NP’s research, all these elements should allow a more complete understanding of organisms and their metabolites.

![Graphical abstract](data/processed/figures/graphical_abstract.png)

je coupe le son

## Flowchart

```mermaid
graph TD
000(adequate minimal input) -->
  100(organism) --> 
    110(scientific names recognition) -->
      120(subtraction of recognized scientific names and translation of vernacular and tcm names to scientific) -->
        130(scientific names recognition) -->
            140(taxonomies accros multiple taxonomy DBs) -->
              |could be discussed, human feeling|150(selection of the best taxonomy) -->
                |could be discussed, human feeling|160(comparison of all obtained taxonomies and upstream filling of taxa choosing best) -->
                  170(sanitized organism with clean taxonomy) -->
999(adequate minimal output)

000(adequate minimal input) --> 
  200(metabolite) --> 
    211(InChI) --> 
      220(InChI)
  200(metabolite) --> 
    212(SMILES) --> 
      220(InChI)
  200(metabolite) --> 
    213(name) --> 
      220(InChI) --> 
        230(ROMOL) --> 
          |think about 2D 3D| 240(sanitized ROMOL) -->
            251(InChI, SMILES, InChIKey) -->
              261(classyfire taxonomy)-->
                270(sanitized metabolite with taxonomy and metadata)
            251(InChI, SMILES, InChIKey) -->
              |still to do| 262(chem-GPS coordinates)-->
                270(sanitized metabolite with taxonomy and metadata)
            251(InChI, SMILES, InChIKey) -->
              |still to do| 263(in silico tandem ms spectra)-->
                270(sanitized metabolite with taxonomy and metadata)
          240(sanitized ROMOL) -->
            252(additional terms xlogP, MF, exact mass) -->
                270(sanitized metabolite with taxonomy and metadata) -->
999(adequate minimal output)

000(adequate minimal input) --> 
  300(reference) -->
    |splitting actual non-optimal| 311(DOIs) -->
      320(generated reference with metadata)
  300(reference) -->
    |splitting actual non-optimal| 312(PubmedIDs )-->
      |think about minimal required fields DOI, authors, date, title, journal?| 320
  300(reference) -->
    |splitting actual non-optimal| 313(text) -->
      320 -->
        |think about it|330(additional cleaning steps to define) -->
          340(sanitized reference with metadata) -->
999(adequate minimal output)
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

