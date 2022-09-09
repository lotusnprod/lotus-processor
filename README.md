# The LOTUS Initiative

*LOTUS* is a comprehensive collection of documented structure-organism pairs.
Within the frame of current computational approaches in Natural Products research and related fields, 
these documented structure-organism pairs should allow a more complete understanding of organisms and their chemistry.

For additional details on this first project of the LOTUS Initiative please have a look at our preprint [https://lotusnprod.github.io/lotus-manuscript](https://lotusnprod.github.io/lotus-manuscript)

As the shorter READMEs are the best, we decided to give detailed documentation in our [Wiki](https://github.com/lotusnprod/lotus-processor/wiki). ❤️

If you are in a hurry and just want a quick test without reading the whole documentation, what you need is:

- R
- Python 3
- Java >= 17

### UNIX systems requirements

Please make sure to have [Make](https://www.gnu.org/software/make) installed.


#### OSX
If you observe this type of errors 
`objc[67570]: +[__NSPlaceholderDate initialize] may have been in progress in another thread when fork() was called.`, you might also need to append the following lines to your bash profile.

- export OBJC_DISABLE_INITIALIZE_FORK_SAFETY=YES
- export DISABLE_SPRING=TRUE

See https://stackoverflow.com/a/52230415 for details.


### Windows systems requirements

As steps are a bit longer, we made a dedicated [Wiki for Windows users](https://github.com/lotusnprod/lotus-processor/wiki/Windows-users).



## Test the processing workflow
```
git clone https://github.com/lotusnprod/lotus-processor.git
cd lotus-processor
conda env create --file environment.yml
conda activate lotus_env
Rscript -e 'remotes::install_github("ropensci/rcrossref")'
make MODE=test lotus-bloom
```

if everything worked smoothly you can then:

```
make MODE=test lotus-check
```

## Data Availability Statements

- The data used to support the findings of this study have been deposited on Zenodo [https://zenodo.org/communities/the-lotus-initiative](https://zenodo.org/communities/the-lotus-initiative).
A snapshot of the repository at the time of publication is also available under the same link.
- Outputs of the lotus-processor will be regularly archived at [https://zenodo.org/record/5665295](https://zenodo.org/record/5665295)


## Dataset list

All data sources used for this study are listed under [docs/dataset.csv](docs/dataset.csv)

## Contributing

Please read [contributing](CONTRIBUTING.md) for details on our code of conduct, and the process for submitting pull requests.

## Authors

- *Adriano Rutz* - _Initial work_ - [Adafede](https://github.com/Adafede)
- *Pierre-Marie Allard* - _Investigator_ - [oolonek](https://github.com/oolonek)
- *Jonathan Bisson* - _Hacker in Chief_ - [bjonnh](https://github.com/bjonnh)

See also the list of [contributors](https://github.com/lotusnprod/lotus-processor/graphs/contributors) who participated in this project.

## License

This project is licensed under the GNU GPLv3 license - see [license file](LICENSE.md) for details
