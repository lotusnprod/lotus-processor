[![pipeline status](https://gitlab.unige.ch/Adriano.Rutz/opennaturalproductsdb/badges/master/pipeline.svg)](https://gitlab.unige.ch/Adriano.Rutz/opennaturalproductsdb/-/commits/master)

# LOTUS: a curated naturaL prOducTs occUrrences databaSe

## Overview

*LOTUS* actually, consists of XXX'XXX referenced structure - organism pairs, collected and curated among 34 open databases (DB).
It represents the most exhaustive open DB of natural products (NP).
It encompasses XXX'XXX distinct sanitized structures in XX’XXX resolved organisms.
Within the frame of current computational approaches to guide NP’s research, all these elements should allow a more complete understanding of organisms and their metabolites.

## Data Availability Statements

The data used to support the findings of this study have been deposited in the XXX Dataverse repository ([DOI or OTHER PERSISTENT IDENTIFIER]).
An snapshot of the repository at the time of publication is also available.

## Dataset list

![dataset](dataset.md)

## Computational requirements

see <https://social-science-data-editors.github.io/guidance/template-README.html>

### Software Requirements

### Memory and Runtime Requirements

## Description of programs

### Instructions

### Details

## List of tables and programs

## References

## Acknowledgments

- Hat tip to anyone whose code was used
- Inspiration
- etc.

## Flowchart

![flowchart](flowchart.md)

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

You need:

- Git
- [DVC](https://www.dvc.org)
- gnfinder (see how to install)
- gnverify (see how to install)

- If you want to be able to use the [PMID translation script](src/2_curating/2_editing/reference/subscripts/1_translating/pubmed.R) with less limitations and it to work correctly, you have to set an API key as described in the following [vignette](https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html)

- You need to have access to the metabomaps server to be able to pull the data, if you don't you will have to pull all the DBs data
from scratch.

### Access to metabomaps

Add this to your

```console
~/.ssh/config
```

```console
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

```console
make docker-build
make docker-bash
```

This will bring you in a container that will already have all the dependencies installed so you can run your commands in it.

### Packages

#### Conda environement

A "loose" environment.yml file is created and should allow to recreate a working env formthe project without beeing too restrictive on the versions to install. Install it by running in the home directory.

We will also create another environment for strict mirroring of the installed packages. (TO DO)

conda env create -f environment.yml

If your environement is not directly sources by your default bash run the following lines

```console
source ~/anaconda3/etc/profile.d/conda.sh

conda activate lotus_env
```

Your R working directory should be 'src'

```console
cd src
```

If you are using Visual Studio be sure to set your R path in the settings option to reflect your created conda environment.

## Minimal working example

A minimal working example containing XXX entries coming from various DB's is proposed.
Use this example to check if all steps are running correctly on your machine.

## Molconvert issue

At the moment, we use molconvert (commercial) for structure to chemical name conversion. Since we cannot disseminate it, you won't be able to proceed to the translation except if you modify following variables in src/paths.R accordingly:

```console
works_locally_only <- TRUE # FALSE
molconvertPath <- adapt_path_to # "~/../../Applications/MarvinSuite/bin/molconvert"
```

## To build a graph from the make

(Requires remake and gprof2dot <https://github.com/jrfonseca/gprof2dot>)

```console
remake --profile -B curating
gprof2dot -f callgrind callgrind.out.50802 | dot -Tpng -o output_full.png
```

### Description

The minimal working example file contains following columns:

## TO UPDATE 1

Give an example

### Final output

If everything went well, the output of the minimal working example should be:

## TO UPDATE 2

Give an example

## Explanations

Add additional notes

## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of conduct, and the process for submitting pull requests to us.

## Authors

- **Adriano Rutz** - _Initial work_ - [Adriano.Rutz](https://gitlab.unige.ch/Adriano.Rutz)
- **Pierre-Marie Allard** - _Investigator_ - [Pierre-Marie.Allard](https://gitlab.unige.ch/Pierre-Marie.Allard)
- **Jonathan Bisson** - _Hacker in Chief_ - [bjo](https://gitlab.unige.ch/bjo)

See also the list of [contributors](https://gitlab.unige.ch/Adriano.Rutz/opennaturalproductsdb/-/project_members) who participated in this project.

## License

This project is licensed under the GNU GPLv3 license - see the [LICENSE.md](LICENSE.md) file for details
